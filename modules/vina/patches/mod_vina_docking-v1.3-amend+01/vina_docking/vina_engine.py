''' AutoDock Vina Engine (caddmodules.workhorse.docking.vina_engine)
'''

import json
import logging
import re
import tempfile
from dataclasses import asdict, dataclass
from io import BytesIO
from pathlib import Path
from shlex import quote as _quote
import time
from typing import (Awaitable, Callable, Generator, List, Optional, Tuple,
                    TypeVar, Union)
from zipfile import ZIP_LZMA, ZipFile

import numpy as np
import pandas as pd
from mxf_files.files import PDBFile, PDBQTFile
from mxf_serum.target.local import async_call, call
from scipy.spatial.distance import cdist as spatial_cdist

from ._vina_engine import COLLISION_SCOPE
from ._vina_engine import detect as collision_detect
from ._vina_engine import restore_ligand as _restore_ligand
from .consts import (ADFR_AUTOGRID4_PATH, MGLTOOLS_AD4_PARMS_PATH,
                     MGLTOOLS_PREPARE_FLEX_RECEPTOR4, MGLTOOLS_PREPARE_GPF4,
                     MGLTOOLS_PREPARE_RECEPTOR4, MGLTOOLS_PROCESS_VINA_RESULT,
                     VINA_1_2_X_SEARCH_PATH)
from .exceptions import (ADFRAutoGridGeneralException,
                         GenerateEmptyFlexibleResidues,
                         MGLToolsPrepGPF4GeneralException,
                         MGLToolsPrepReceptor4GeneralException,
                         MGLToolsProcessVinaResultGeneralException,
                         RemoteProcedureCallFailed,
                         VinaAtomTypeNotSupport)
from .merge import MergeFlavor, pdb_merge


def quote(s: Union[str, Path]) -> str:
    if isinstance(s, Path):
        s = str(s)
    return _quote(s)


SerumCallReturn = TypeVar('SerumCallReturn')


GPF4_LIGAND_TYPES_PARTTEN = re.compile(r'(<?ligand_types)(.*?)(=?# ligand atom types)')

VINA_DOCKING_LOG = 'docking.log'
VINA_ERROR_LOG = 'docking.err'
VINA_SUCCESS_SIGNAL = 'dock_success'
VIN_COMMANDLINE_PIPE = f'> {VINA_DOCKING_LOG} 2> {VINA_ERROR_LOG} && touch {VINA_SUCCESS_SIGNAL}'


class Vina112:
    ''' engine for lagency AutoDock Vina 1.1.2
    '''
    pass


@dataclass
class Vina123Parms:
    center_x: float  # X coordinate of the center (Angstrom)
    center_y: float  # Y coordinate of the center (Angstrom)
    center_z: float  # Z coordinate of the center (Angstrom)
    size_x: float  # size in the X dimension (Angstrom)
    size_y: float  # size in the Y dimension (Angstrom)
    size_z: float  # size in the Z dimension (Angstrom)
    scoring: str = 'vina'   # scoring function (ad4, vina or vinardo)
    maps: str = ''  # affinity maps for the autodock4.2 (ad4) or vina scoring function, stem of map eg. 1iep_receptor for 1iep_receptor.*.map

    # these weight is written here just for completeness, not using for now
    #
    # weight_gauss1: float = -0.035579  # gauss_1 weight
    # weight_gauss2: float = -0.005156  # gauss_2 weight
    # weight_repulsion: float = 0.84024500000000002  # repulsion weight
    # weight_hydrophobic: float = -0.035069000000000003  # hydrophobic weight
    # weight_hydrogen: float = -0.58743900000000004  # Hydrogen bond weight
    # weight_rot: float = 0.058459999999999998  # N_rot weight
    # weight_vinardo_gauss1: float = -0.044999999999999998  # Vinardo gauss_1 weight
    # weight_vinardo_repulsion: float = 0.80000000000000004  # Vinardo repulsion weight
    # weight_vinardo_hydrophobic: float = -0.035000000000000003  # Vinardo hydrophobic weight
    # weight_vinardo_hydrogen: float = -0.59999999999999998  # Vinardo Hydrogen bond weight
    # weight_vinardo_rot: float = 0.058459999999999998  # Vinardo N_rot weight
    # weight_ad4_vdw: float = 0.16619999999999999  # ad4_vdw weight
    # weight_ad4_hb: float = 0.12089999999999999  # ad4_hb weight
    # weight_ad4_elec: float = 0.1406  # ad4_elec weight
    # weight_ad4_dsolv: float = 0.13220000000000001  # ad4_dsolv weight
    # weight_ad4_rot: float = 0.29830000000000001  # ad4_rot weight
    # weight_glue: float = 50  # macrocycle glue weight

    cpu: int = 0  # the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1
    seed: int = 0  # explicit random seed
    exhaustiveness: int = 8  # exhaustiveness of the global search (roughly proportional to time): 1+
    max_evals: int = 0  # number of evaluations in each MC run (if zero, which is the default, the number of MC steps is based on heuristics
    num_modes: int = 9  # maximum number of binding modes to generate
    min_rmsd: float = 1.  # minimum RMSD between output poses
    energy_range: float = 3.  # maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)
    spacing: float = 0.375  # grid spacing (Angstrom)
    verbosity: int = 1  # verbosity (0=no output, 1=normal, 2=verbose)

    def __str__(self) -> str:
        ''' convert to command line arguements
        '''
        return ''.join([f' --{k} {quote(str(v))}' for k, v in asdict(self).items() if v])


class Vina123:
    ''' engine for Modern AutoDock Vina 1.2.3
    '''
    #######################
    # Prepare protein
    #######################

    @classmethod
    def prepare_receptor(cls: 'Vina123', receptor: Path, flex: Optional[List[str]] = None, *args, logger: logging.Logger = logging.getLogger('vina_engine.vina123.prepare_receptor')) -> Tuple[str, Union[None, str]]:
        """ prepare ligand via MGLTools

        Args:
            receptor (Path): Receptor pdb file path
            flex (Optional[List[str]], optional): list of flex residue in form of `<chain id>:<residue name><residue id>` eg. `A:ARG1`. Defaults to None.
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.prepare_receptor').

        Raises:
            MGLToolsPrepLigand4GeneralException: prepare_receptor4 failed

        Returns:
            Tuple[str, Union[None, str]]: rigid part and flexible part of receptor in pdbqt format
        """
        with tempfile.TemporaryDirectory() as workdir:
            receptor_path = Path(workdir).joinpath(receptor.name)
            receptor_path.symlink_to(receptor)
            rigid_path = receptor_path.with_name('rigid.pdbqt')
            flex_path = receptor_path.with_name('flex.pdbqt')
            flex = ','.join(flex) if isinstance(flex, list) else None

            if flex:
                # prepare_flexreceptor4 need a pdbqt format file
                if receptor_path.suffix != '.pdbqt':
                    pdbqt_receptor, _ = cls.prepare_receptor(receptor, logger=logger)
                    receptor_path = receptor_path.with_suffix('.pdbqt')
                    receptor_path.write_text(pdbqt_receptor)

                cmd = f'{MGLTOOLS_PREPARE_FLEX_RECEPTOR4} -r {quote(receptor_path)} -s {flex} -g {quote(rigid_path)} -x {quote(flex_path)}'
                proc = call(cmd, cwd=workdir, logger=logger)
                if rigid_path.exists() and flex_path.exists():
                    if len(flex_path.read_text().strip()) == 0:
                        raise GenerateEmptyFlexibleResidues('prepare_receptor4 generate empty flexible pdbqt file')
                    logger.info(f'rigid and flexible pdbqt file generated for receptor ({receptor})')
                else:
                    logger.error(f'pdbqt file fail to generated for receptor ({receptor})')
                    logger.debug(f'STDOUT: {proc["stdout"]}')
                    logger.debug(f'STDERR: {proc["stderr"]}')
                    raise MGLToolsPrepReceptor4GeneralException('prepare_receptor4 script in MGLTools fail to generate pdbqt file (General Error)')

                return rigid_path.read_text(), flex_path.read_text()

            else:
                cmd = f'{MGLTOOLS_PREPARE_RECEPTOR4} -r {quote(receptor_path)} -o {quote(rigid_path)} -A hydrogens -U nphs_lps_waters_nonstdres_deleteAltB'
                proc = call(cmd, cwd=workdir, logger=logger)
                if rigid_path.exists():
                    logger.info(f'pdbqt file generated for receptor ({receptor})')
                else:
                    logger.error(f'pdbqt file fail to generated for receptor ({receptor})')
                    logger.debug(f'STDOUT: {proc["stdout"]}')
                    logger.debug(f'STDERR: {proc["stderr"]}')
                    raise MGLToolsPrepReceptor4GeneralException('prepare_receptor4 script in MGLTools fail to generate pdbqt file (General Error)')

                return rigid_path.read_text(), None

    #######################
    # Proceeding docking
    #######################

    @staticmethod
    async def gen_affinity_map(ligand: Path, receptor: Path, center_xyz: Tuple[float, float, float], *args, with_water: bool, with_flex: bool, logger: logging.Logger = logging.getLogger('vina_engine.vina123.gen_affinity_map')) -> bytes:
        """ Generate the affinity map for Autodock4 scoring function.

        Args:
            ligand (Path): Ligand pdbqt file path
            receptor (Path): Receptor pdbqt file path
            center_xyz (Tuple[float, float, float]): The center coordinate of searching area
            with_water (bool): wether docking with water
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.gen_affinity_map').

        Raises:
            MGLToolsPrepGPF4GeneralException: prepare_gpf4 failed
            ADFRAutoGridGeneralException: autogrid4 failed
            VinaExmapleScriptGeneralException: mapwater.py failed

        Returns:
            bytes: affinity map zip file content
        """
        with tempfile.TemporaryDirectory() as workdir:
            ligand_path = Path(workdir).joinpath('ligand.pdbqt')
            receptor_path = Path(workdir).joinpath('receptor.pdbqt')
            receptor_gpf = receptor_path.with_suffix('.gpf')
            ligand_path.symlink_to(ligand)
            receptor_path.symlink_to(receptor)

            vec_center = f'{center_xyz[0]:.3f},{center_xyz[1]:.3f},{center_xyz[2]:.3f}'

            if with_water:
                raise NotImplementedError

            else:
                # generate gpf file
                logger.debug('start to generate gpf file')
                cmd = f'{MGLTOOLS_PREPARE_GPF4} -l {quote(ligand_path)} -r {quote(receptor_path)} -p parameter_file={quote(MGLTOOLS_AD4_PARMS_PATH)} -p gridcenter={vec_center} -o {quote(receptor_gpf)}'
                if with_flex:
                    # flexible residue will counting inside of ligand
                    # expanding the atom type to full from AD4.1_bound.dat
                    cmd += ' -p ligand_types=H,HD,HS,C,A,N,NA,NS,OA,OS,F,Mg,MG,P,SA,S,Cl,CL,Ca,CA,Mn,MN,Fe,FE,Zn,ZN,Br,BR,I,Z'
                proc = await async_call(cmd, cwd=workdir, logger=logger)
                if receptor_gpf.exists():
                    logger.info(f'gpf file generated for receptor ({receptor}) and ligand ({ligand})')
                else:
                    logger.error(f'gpf file fail to generated for receptor ({receptor}) and ligand ({ligand})')
                    logger.debug(f'STDOUT: {proc["stdout"]}')
                    logger.debug(f'STDERR: {proc["stderr"]}')
                    raise MGLToolsPrepGPF4GeneralException

                # run autogrid4
                cmd = f'{ADFR_AUTOGRID4_PATH} -p {quote(receptor_gpf)} && touch autogrid4_success'
                proc = await async_call(cmd, cwd=workdir, logger=logger)
                if Path(workdir).joinpath('autogrid4_success').exists():
                    logger.info(f'autogrid completed for receptor ({receptor}) and ligand ({ligand})')
                else:
                    logger.error(f'autogrid fail for receptor ({receptor}) and ligand ({ligand})')
                    logger.info('usually it caused by mapping atoms failure')
                    logger.debug(f'STDOUT: {proc["stdout"]}')
                    logger.debug(f'STDERR: {proc["stderr"]}')
                    raise ADFRAutoGridGeneralException

            # grab *.map* files
            buf = BytesIO()
            with ZipFile(buf, mode='w', compression=ZIP_LZMA) as tarball:
                for file in Path(workdir).glob('*.map*'):
                    tarball.write(str(file), arcname=str(file.name))

            # usually it will not be to big, but still leave some foot print
            logger.info(f'grid map collected, {buf.tell()} bytes.')
            buf.seek(0)
            b_content = buf.read()
            buf.close()
            return b_content

    @classmethod
    async def docking_simple(cls, async_call: Callable[[str, str], Awaitable[SerumCallReturn]], receptor_pdbqt: Path, ligand_pdbqt: Path, ref_receptor_pdb: Path, ref_ligand_pdb: Path, workdir: Path, docking_args: Vina123Parms, logger: logging.Logger = logging.getLogger('vina_engine.vina123.dock_simple')) -> Tuple[SerumCallReturn, pd.DataFrame]:
        """ simple docking subroutine

        Args:
            async_call (Callable[[str, str], Awaitable[SerumCallReturn]]): mxf_serum.target.Runner.async_call or other compatible callback
            receptor_pdbqt (Path): path of prepared receptor in pdbqt format
            ligand_pdbqt (Path): path of prepared ligand in pdbqt format
            ref_receptor_pdb (Path): path of original receptor in pdb format
            ref_ligand_pdb (Path): path of prepared ligand in pdb format
            workdir (Path): working directory to store docking result
            docking_args (Vina123Parms): docking parameters
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.dock_simple').

        Returns:
            Tuple[SerumCallReturn, pd.DataFrame]: async_call returned `proc`, along with parsed docking result, docked complex is store at column `_com_fn`
        """
        receptor = workdir.joinpath('receptor.pdbqt')
        receptor.symlink_to(receptor_pdbqt)
        ligand = workdir.joinpath('ligand.pdbqt')
        ligand.symlink_to(ligand_pdbqt)

        if docking_args.scoring in ['ad4', ]:
            zipped_map = await cls.gen_affinity_map(ligand, receptor, (docking_args.center_x, docking_args.center_y, docking_args.center_z), with_water=False, with_flex=False, logger=logger)
            zipped_map_buf = BytesIO(zipped_map)
            with ZipFile(zipped_map_buf) as zipped_map_fp:
                zipped_map_fp.extractall(workdir)
            zipped_map_buf.close()
            del zipped_map
            docking_args.maps = 'receptor'
            cmd = f'PATH={quote(VINA_1_2_X_SEARCH_PATH)}:$PATH vina --ligand ligand.pdbqt {docking_args} {VIN_COMMANDLINE_PIPE}'
        else:
            cmd = f'PATH={quote(VINA_1_2_X_SEARCH_PATH)}:$PATH vina --receptor receptor.pdbqt --ligand ligand.pdbqt {docking_args} {VIN_COMMANDLINE_PIPE}'

        try:
            proc = await async_call(cmd, str(workdir))
        except RuntimeError:
            raise RemoteProcedureCallFailed('Calculation submit failed (Sugon/Changchun/...)')
        if workdir.joinpath(VINA_SUCCESS_SIGNAL).exists():
            logger.info(f'docking success for {receptor_pdbqt} receptor and {ligand_pdbqt} ligand')
            dock_df = await cls.split_docked(workdir.joinpath('ligand_out.pdbqt'), receptor, logger=logger)
            dock_df['_lig_fn'] = dock_df['_lig_fn'].apply(lambda x: str(workdir.joinpath(x)))
            dock_df['_rec_fn'] = str(receptor)

            # create complex pdb file
            unpatched_ligands_pdbqt = [Path(fn) for fn in dock_df['_lig_fn'].to_list()]
            docked_complex_pdb = []
            # 1. restore ligand hydrogen and connect
            patched_ligands = cls.restore_ligand(ref_ligand_pdb, ligand, unpatched_ligands_pdbqt, logger=logger.getChild('restore_ligand'))

            # 2. detect collision
            collision_detected = []
            for ligands_path, patched_pdb_content in zip(unpatched_ligands_pdbqt, patched_ligands):
                buf = pdb_merge(ref_receptor_pdb.read_text(), patched_pdb_content, MergeFlavor.MACRO_LIGAND_DOCKING, logger=logger)

                complex_path = ligands_path.with_suffix('').with_suffix('.complex.pdb')
                complex_path.write_text(buf)
                docked_complex_pdb.append(str(complex_path))
                del buf

                if (
                    len(collision_detect([patched_pdb_content], COLLISION_SCOPE.inter)) > 0
                ) or (
                    len(collision_detect([ref_receptor_pdb.read_text(), patched_pdb_content], COLLISION_SCOPE.intra)) > 0
                ):
                    logger.warning('Collision detected, this dock mode is going to be considered as invalid')
                    collision_detected.append(True)
                else:
                    collision_detected.append(False)

            dock_df['_com_fn'] = docked_complex_pdb
            dock_df['_collision'] = collision_detected

            interactions = cls.prepare_interaction(ref_receptor_pdb.read_text(), patched_ligands, workdir, logger=logger.getChild('prepare_interaction'))
            dock_df['_interaction'] = [str(i) for i in interactions]

            return proc, dock_df
        else:
            logger.error(f'docking fail for {receptor_pdbqt} receptor and {ligand_pdbqt} ligand')

            if workdir.joinpath(VINA_ERROR_LOG).exists():
                logger.debug('error log spot, start to parse')
                cls.vina_error_parse(proc['stderr'], workdir.joinpath(VINA_ERROR_LOG).read_text())

            if proc['code'] == 0:
                proc['code'] = 1
            raise RuntimeError('unclassified vina error')

    @classmethod
    async def docking_flex(cls, async_call: Callable[[str, str], Awaitable[SerumCallReturn]], receptor_rigid_pdbqt: Path, receptor_flex_pdbqt: Path, ligand_pdbqt: Path, ref_receptor_pdb: Path, ref_ligand_pdb: Path, workdir: Path, docking_args: Vina123Parms, logger: logging.Logger = logging.getLogger('vina_engine.vina123.dock_simple')) -> Tuple[SerumCallReturn, pd.DataFrame]:
        """ flexible docking subroutine

        Args:
            async_call (Callable[[str, str], Awaitable[SerumCallReturn]]): mxf_serum.target.Runner.async_call or other compatible callback
            receptor_rigid_pdbqt (Path): path of prepared rigid receptor in pdbqt format
            receptor_flex_pdbqt (Path): path of prepared flexible receptor in pdbqt format
            ligand_pdbqt (Path): path of prepared ligand in pdbqt format
            ref_receptor_pdb (Path): path of original receptor in pdb format
            ref_ligand_pdb (Path): path of prepared ligand in pdb format
            workdir (Path): working directory to store docking result
            docking_args (Vina123Parms): docking parameters
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.dock_simple').

        Returns:
            Tuple[SerumCallReturn, pd.DataFrame]: async_call returned `proc`, along with parsed docking result, docked complex is store at column `_com_fn`
        """
        receptor_rigid = workdir.joinpath('receptor.rigid.pdbqt')
        receptor_rigid.symlink_to(receptor_rigid_pdbqt)
        receptor_flex = workdir.joinpath('receptor.flex.pdbqt')
        receptor_flex.symlink_to(receptor_flex_pdbqt)
        ligand = workdir.joinpath('ligand.pdbqt')
        ligand.symlink_to(ligand_pdbqt)

        if docking_args.scoring in ['ad4', ]:
            zipped_map = await cls.gen_affinity_map(ligand, receptor_rigid, (docking_args.center_x, docking_args.center_y, docking_args.center_z), with_water=False, with_flex=True, logger=logger)
            zipped_map_buf = BytesIO(zipped_map)
            with ZipFile(zipped_map_buf) as zipped_map_fp:
                zipped_map_fp.extractall(workdir)
            zipped_map_buf.close()
            del zipped_map
            docking_args.maps = 'receptor'
            cmd = f'PATH={quote(VINA_1_2_X_SEARCH_PATH)}:$PATH vina --flex receptor.flex.pdbqt --ligand ligand.pdbqt {docking_args} {VIN_COMMANDLINE_PIPE}'
        else:
            cmd = f'PATH={quote(VINA_1_2_X_SEARCH_PATH)}:$PATH vina --receptor receptor.rigid.pdbqt --flex receptor.flex.pdbqt --ligand ligand.pdbqt {docking_args} {VIN_COMMANDLINE_PIPE}'

        try:
            proc = await async_call(cmd, str(workdir))
        except RuntimeError:
            raise RemoteProcedureCallFailed('Calculation submit failed (Sugon/Changchun/...)')
        if workdir.joinpath(VINA_SUCCESS_SIGNAL).exists():
            logger.info(f'docking success for {receptor_rigid_pdbqt}/{receptor_flex_pdbqt} receptor and {ligand_pdbqt} ligand')
            dock_df = await cls.split_docked(workdir.joinpath('ligand_out.pdbqt'), receptor_rigid, logger=logger)
            dock_df['_lig_fn'] = dock_df['_lig_fn'].apply(lambda x: str(workdir.joinpath(x)))
            dock_df['_rec_fn'] = str(receptor_rigid)

            # create complex pdb file
            unpatched_ligands_pdbqt = [Path(fn) for fn in dock_df['_lig_fn'].to_list()]
            docked_complex_pdb = []
            # 1. restore ligand hydrogen and connect
            patched_ligands = cls.restore_ligand(ref_ligand_pdb, ligand, unpatched_ligands_pdbqt, logger=logger.getChild('restore_ligand'))
            # 2. resotre flexible receptor atoms
            patched_receptors = cls.restore_receptor(ref_receptor_pdb, receptor_flex, unpatched_ligands_pdbqt)

            # 3. detect collision
            collision_detected = []
            for ligands_path, patched_ligand_content, patched_receptor_content in zip(unpatched_ligands_pdbqt, patched_ligands, patched_receptors):
                buf = pdb_merge(patched_receptor_content, patched_ligand_content, MergeFlavor.MACRO_LIGAND_DOCKING, logger=logger)

                complex_path = ligands_path.with_suffix('').with_suffix('.complex.pdb')
                complex_path.write_text(buf)
                docked_complex_pdb.append(str(complex_path))
                del buf

                if (
                    len(collision_detect([patched_ligand_content], COLLISION_SCOPE.inter)) > 0
                ) or (
                    len(collision_detect([patched_receptor_content, patched_ligand_content], COLLISION_SCOPE.intra)) > 0
                ):
                    logger.warning('Collision detected, this dock mode is going to be considered as invalid')
                    collision_detected.append(True)
                else:
                    collision_detected.append(False)

            dock_df['_com_fn'] = docked_complex_pdb
            dock_df['_collision'] = collision_detected

            interactions = [
                cls.prepare_interaction(single_patched_receptor, [single_patched_ligand], workdir, logger=logger.getChild('prepare_interaction'))
                for single_patched_receptor, single_patched_ligand in zip(patched_receptors, patched_ligands)
            ]
            dock_df['_interaction'] = [str(i[0]) if i else '' for i in interactions]

            return proc, dock_df
        else:
            # filter known issue

            # 1. OOM might happen when flexible residue is too much
            # TODO: this should use proc['code'], and will implement in mxf_serum later
            if 'oom-kill event' in proc['stderr'].decode():
                logger.warning('OOM event is detected')
                raise RuntimeError('too many flexible residue, system memory is insufficient.')

            if workdir.joinpath(VINA_ERROR_LOG).exists():
                logger.debug('error log spot, start to parse')
                cls.vina_error_parse(proc['stderr'], workdir.joinpath(VINA_ERROR_LOG).read_text())

            logger.error(f'docking fail for {receptor_rigid_pdbqt}/{receptor_flex_pdbqt} receptor and {ligand_pdbqt} ligand')
            logger.debug(f'stdout: {proc["stdout"]}')
            logger.debug(f'stderr: {proc["stderr"]}')
            raise RuntimeError('unclassified vina error')

    @staticmethod
    def vina_error_parse(stderr: str, vina_err: str, *, logger: logging.Logger = logging.getLogger('vina_engine.vina123.vina_error_parse')) -> None:
        """ parse vina stderr and raise proper exception

        Args:
            stderr (str): content of stderr
            vina_err (str): content of stderr from vina
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.vina_error_parse').

        Raises:
            VinaAtomTypeNotSupport: given atom type is not supported

        Returns:
            None
        """
        # judge atom type error case
        result = re.search(r'Atom type ([^\s]*) is not a valid AutoDock type', vina_err)
        logger.debug(f'try to judge atom type error... {result}')
        if result is not None:
            invalid_atom_name = result.groups()[0]
            raise VinaAtomTypeNotSupport(f'Atom type {invalid_atom_name} is not a valid AutoDock type')

        logger.debug(f'unclassified stderr: {stderr}')
        return None

    @staticmethod
    async def split_docked(docked_lig: Path, receptor_: Path, logger: logging.Logger = logging.getLogger('vina_engine.vina123.split_docked')) -> pd.DataFrame:
        """split docked ligand pdbqt file, store at the same place of `docked_lig`

        Args:
            docked_lig (Path): docked lignad file in pdbqt format
            receptor_ (Path): receptor file, used to analysis contract
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.split_docked').

        Returns:
            pd.DataFrame: parsed docking result with column `_lig_fn` point to splitted docked ligand
        """
        with tempfile.TemporaryDirectory() as workdir:
            ligand = Path(workdir).joinpath(docked_lig.name)
            ligand.symlink_to(docked_lig)
            receptor = Path(workdir).joinpath(receptor_.name)
            receptor.symlink_to(receptor_)

            cmd = f'{MGLTOOLS_PROCESS_VINA_RESULT} -f {ligand.name} -r {receptor.name} -p'
            proc = await async_call(cmd, workdir, logger=logger)
            if proc['code'] != 0:
                logger.warning(f'process_vina.py in MGLTools fail with code {proc["code"]}')
                logger.debug(f'STDERR: {proc["stderr"]}')
                logger.debug(f'STDOUT: {proc["stdout"]}')

            col_dock_mode = []
            col_affinty = []
            col_rmsd_lb = []
            col_rmsd_ub = []
            col_fn = []
            for splitted_pdbqt in Path(workdir).glob('*_model*.pdbqt'):
                content = splitted_pdbqt.read_text()
                dock_mode_pattern = re.search(r'USER\s+AD>\s+(\d+)\s+of\s+\d+\s+MODELS', content)
                if dock_mode_pattern is None:
                    raise MGLToolsProcessVinaResultGeneralException()
                else:
                    dock_mode = dock_mode_pattern.groups()
                if len(dock_mode) == 0:
                    logger.error('fail to extract dock mode')
                    logger.debug(f'edge case of splitted docked pdbqt file:\n{content}')
                dock_mode = dock_mode[0]
                dock_status = re.search(r'USER\s+AD>\s+(-?\d+\.\d+),\s+((-?\d+\.\d+)|0),\s+((-?\d+\.\d+)|0)', content).groups()
                if len(dock_status) < 5:
                    logger.error('fail to extract dock status')
                    logger.debug(f'edge case of splitted docked pdbqt file:\n{content}')
                affinty, rmsd_lb, _, rmsd_ub, _ = dock_status

                report_affinity = float(affinty)
                if abs(report_affinity) > 100:
                    logger.warning(f'reporting affinity ({report_affinity}) is over ±100, which is wired, discarded.')
                    continue

                col_dock_mode.append(int(dock_mode))
                col_affinty.append(float(affinty))
                col_rmsd_lb.append(float(rmsd_lb))
                col_rmsd_ub.append(float(rmsd_ub))
                col_fn.append(splitted_pdbqt.name)
                docked_lig.with_name(splitted_pdbqt.name).write_text(content)

            return pd.DataFrame({
                'mode': col_dock_mode,
                'affinity (kcal/mol)': col_affinty,
                'dist from best mode rmsd l.b.': col_rmsd_lb,
                'dist from best mode rmsd u.b.': col_rmsd_ub,
                '_lig_fn': col_fn,
            })

    @staticmethod
    def restore_receptor(pdbfn: Path, flexpdbqt: Path, patchfn: List[Path], logger: logging.Logger = logging.getLogger('vina_engine.vina123.restore_ligand')) -> Generator[str, None, None]:
        origin_pdb = PDBFile(pdbfn)
        origin_atoms = list(origin_pdb.atoms)
        origin_coords = origin_pdb.coords
        prepared_pdbqt = PDBQTFile(flexpdbqt)
        prepared_atoms = list(prepared_pdbqt.atoms)
        rearrange = prepared_pdbqt.coords

        # it reflect from prepared flexible part of macros to rigid one
        reflect = np.argwhere(np.isclose(spatial_cdist(rearrange, origin_coords), 0))
        if reflect.shape[0] < rearrange.shape[0]:
            logger.warning(f'#{rearrange.shape[0]} flexible sidechain atoms, but #{reflect.shape[0]} reflection found. some hydrogen may lost.')

        flex_resi = set([(atom.nameChain, atom.idxResi) for atom in [prepared_atoms[i] for i in reflect[:, 0]]])
        # find atoms is flexible residue but not record in pdbqt file as `flex_vague`
        flex_vague = set([idx for idx, atom in enumerate(origin_atoms) if (atom.nameChain, atom.idxResi) in flex_resi])
        flex_vague = set([idx for idx in flex_vague if origin_atoms[idx].nameAtom not in ['O', 'N', 'C', 'CA']])  # mainchain is known to be rigid
        flex_vague = flex_vague - set(reflect[:, 1])
        flex_atom_number = rearrange.shape[0]

        patched_content = []
        for fn in patchfn:
            patch_pdbqt = PDBQTFile(fn)
            patch_atoms = list(patch_pdbqt.atoms)[-flex_atom_number:]
            origin_atoms = list(origin_pdb.atoms)
            for floatted, sattled in reflect:
                origin_atoms[sattled].fX = patch_atoms[floatted].fX
                origin_atoms[sattled].fY = patch_atoms[floatted].fY
                origin_atoms[sattled].fZ = patch_atoms[floatted].fZ
            buf = '\n'.join([atom.to_pdb() for atom in origin_atoms])
            buf += '\n'
            buf += '\n'.join([conn.to_pdb() for conn in origin_pdb.connects])
            buf += '\nEND\n'
            patched_content.append(buf)

        return patched_content

    @staticmethod
    def restore_ligand(pdbfn: Path, pdbqtfn: Path, patchfn: List[Path], *, molstar_compatible: bool = True, molstar_unknown_as: str = '***', logger: logging.Logger = logging.getLogger('vina_engine.vina123.restore_ligand')) -> List[str]:
        """a workaround to revert the PDBQT atom sequence rearrangement

        subroutine prepare_ligand4.py inside MGLTools would distrube the atom order and
        discard non-polar hydrogen this function is intend to recover original atom order
        and pulling its connection back.

        Args:
            pdbfn (Path): full path to original PDB file
            pdbqtfn (Path): full path to docked PDBQT file
            patchfn (List[Path]): full path to patched docked pdb/pdbqt files
            molstar_compatible (bool, optional): if a molstar compatible pdb file is needed. Defaults to True.
            molstar_unknown_as (str, optional): cast molstar unknow residue to this name. Defaults to '***'.
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.restore_ligand').

        Raises:
            NotImplementedError: provided file to patch is not pdb nor pdbqt file

        Returns:
            List[str]: converted pdb block
        """
        time_start_restore = time.time()
        try:
            patched_pdbs: List[str] = _restore_ligand(
                pdbfn.read_text(),
                pdbqtfn.read_text(),
                [f.read_text() for f in patchfn],
            )
        except Exception as e:
            logger.error('fail to restore ligand', exc_info=e)
            raise e
        else:
            logger.info(f'ligand restored in {time.time() - time_start_restore:.3f} sec(s)')

        # map the residue name to molstar compatible name
        if molstar_compatible:
            if molstar_unknown_as[:3] != molstar_unknown_as:
                logger.warning(f'residue name `{molstar_unknown_as}` is requested, which is over 3 char')
            molstar_unknown_as = f'{molstar_unknown_as[:3]:>3}'
            logger.info(f'residue name `{molstar_unknown_as}` will be used to replace `UNK` `UNL` `UNX` and `N`')

            for patched_pdb_idx, patched_pdb in enumerate(patched_pdbs):
                renamed_pdb = []
                for line in patched_pdb.splitlines():
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        if line[17:20] in ['UNK', 'UNL', 'UNX', 'N']:
                            line = line[:17] + molstar_unknown_as + line[20:]
                    renamed_pdb.append(line)
                patched_pdbs[patched_pdb_idx] = '\n'.join(renamed_pdb)

        return patched_pdbs

    #######################
    # Post docking analysis
    #######################

    @staticmethod
    def prepare_interaction(receptor: str, ligands: List[str], cache_file_root: Path, *, logger: logging.Logger = logging.getLogger('vina_engine.vina123.prepare_interaction')) -> List[Path]:
        """ prepare interaction cache file for further plotting generation

        Args:
            receptor (str): pdb file content of receptor
            ligands (List[str]): pdb file content of ligands
            cache_file_root (Path): root path for cache file storage
            logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.prepare_interaction').

        Raises:
            RuntimeError: _description_

        Returns:
            List[Path]: cache files one by one to the `ligands`
        """
        wd = str(cache_file_root.resolve())
        cache_path = []
        for ligand in ligands:
            with tempfile.NamedTemporaryFile('w+', dir=wd, suffix='.2d_interaction.cache', delete=False) as cache_file:
                cache_file.write(json.dumps({
                    'receptor': receptor,
                    'ligand': ligand,
                }))
                cache_path.append(Path(cache_file.name).resolve())
        return cache_path
