''' Pocket detection (caddmodules.workhorse.docking.pocket)
'''

import itertools
import logging
import re
import tempfile
from pathlib import Path
from shlex import quote
from typing import List, Optional, Tuple

from mxf_files.files import PDBFile
from mxf_files.selection import BoundingShape, Selection
from mxf_serum.target.local import call
from typing_extensions import Protocol

from .consts import (ADFR_AUTOSITE_PATH, FPOCKET_PATH,
                     MGLTOOLS_PREPARE_RECEPTOR4)


class PocketInfo:
    pass


class PocketDetectorProt(Protocol):
    '''protocol of pocket detector
    '''
    def __init__(self, pdb: Path, ligand_identification: Optional[str] = None, logger: logging.Logger = logging.getLogger('docking.pocket')) -> None:
        ...

    def search(self, number: int = 1, bounding_shape: BoundingShape = BoundingShape.AABB, vdw_volume: Tuple[float, float] = (0., 1.e9)) -> List[Tuple[Selection, PocketInfo]]:
        """ search for the pocket

        Args:
            number (int, optional): maximal pocket to identify. Defaults to 1.
            bounding_shape (BoundingShape, optional): bounding shape, used to optimize coordinate system. Defaults to BoundingShape.AABB.
            vdw_volume (Tuple[float, float], optional): strict the pocket volume in this range in unit A^3. Defaults to (0., 1.e9).

        Raises:
            NotImplementedError: ...

        Returns:
            List[Tuple[Selection, PocketInfo]]: list of pockets, might be empty
        """
        raise NotImplementedError


class FPocket:
    def __init__(self, pdb: Path, ligand_identification: Optional[str] = None, logger: logging.Logger = logging.getLogger('docking.pocket')) -> None:
        self.logger = logger.getChild(self.__class__.__name__)
        self.whole_pdb = pdb

    def search(self, number: int = 1, bounding_shape: BoundingShape = BoundingShape.AABB, vdw_volume: Tuple[float, float] = (0., 1.e9)) -> List[Tuple[Selection, PocketInfo]]:
        with tempfile.TemporaryDirectory() as workdir:
            pdb_path = Path(workdir).joinpath(f'target{self.whole_pdb.suffix}')
            pdb_path.symlink_to(self.whole_pdb)
            output_dir = Path(workdir).joinpath('target_out').joinpath('pockets')

            cmd = f'{FPOCKET_PATH} -f {quote(pdb_path.name)} -w pdb'
            proc = call(cmd, cwd=workdir, logger=self.logger)
            if output_dir.exists():
                self.logger.info(f'FPocket search complete for {self.whole_pdb}')
            else:
                self.logger.error(f'egde case for FPocket to seach {self.whole_pdb}')
                self.logger.debug(f'STDERR: {proc["stderr"]}')
                self.logger.debug(f'STDOUT: {proc["stdout"]}')

            # pocket will in form of
            #   `target_out/pockets/pocket<id>_atm.pdb` to show neighbor atoms
            #   `target_out/pockets/pocket<id>_vert.pqr` to represent pocket
            result = []
            for idx in itertools.count(1):
                pqr_path = output_dir.joinpath(f'pocket{idx}_vert.pqr')
                if not pqr_path.exists():
                    self.logger.info(f'#{idx - 1} pocket(s) extracted')
                    break

                # check the pocket volume
                volume = re.search(r'Real volume[^\d]*\s+(\d+\.\d+)', pqr_path.read_text())
                if volume is None:
                    self.logger.warning(f'cannot find volume of pocket {idx}')
                    continue
                volume = float(volume.groups()[0])
                vdw_volume = sorted(vdw_volume)
                if not (vdw_volume[0] <= volume <= vdw_volume[1]):
                    self.logger.info(f'volume of pocket {idx} ({volume}) out of range of {vdw_volume}, discarded.')
                    continue

                pdbfile = PDBFile.cast(pqr_path.read_text())
                bounding = Selection.from_luck('system', pdbfile=pdbfile)
                result.append((bounding, PocketInfo()))
            else:
                self.logger.info(f'#{idx} pocket(s) extracted')
            return result[:number]


class AutoSite:
    def __init__(self, pdb: Path, ligand_identification: Optional[str] = None, logger: logging.Logger = logging.getLogger('docking.pocket')) -> None:
        self.logger = logger.getChild(self.__class__.__name__)
        self.whole_pdb = pdb

    def search(self, number: int = 1, bounding_shape: BoundingShape = BoundingShape.AABB, vdw_volume: Tuple[float, float] = (0., 1.e9)) -> List[Tuple[Selection, PocketInfo]]:

        with tempfile.TemporaryDirectory() as workdir:
            pdbqt_path = Path(workdir).joinpath('target.pdbqt')
            output_dir = Path(workdir).joinpath('target')

            if self.whole_pdb.suffix == '.pdb':
                self.logger.info('conversion from pdb to pdbqt is performed via MGLTools/prepare_receptor4.py')
                cmd = f'{MGLTOOLS_PREPARE_RECEPTOR4} -r {quote(str(self.whole_pdb))} -A hydrogens -U nphs_lps_waters_nonstdres_deleteAltB -o {quote(str(pdbqt_path))}'
                proc = call(cmd, cwd=workdir, logger=self.logger)
                if pdbqt_path.exists():
                    self.logger.info(f'conversion from {self.whole_pdb} to {pdbqt_path} success.')
                else:
                    self.logger.info(f'conversion from {self.whole_pdb} to {pdbqt_path} failed.')
                    self.logger.debug(f'STDERR: {proc["stderr"]}')
                    self.logger.debug(f'STDOUT: {proc["stdout"]}')
                    raise RuntimeError('prepare_receptor4.py fail to conversion receptor from pdb to pdbqt')
            elif self.whole_pdb.suffix == '.pdbqt':
                pdbqt_path.symlink_to(self.whole_pdb)
            else:
                raise NotImplementedError(f'request a {self.whole_pdb.suffix} file, while only .pdb/.pdbqt file is support')

            cmd = f'{ADFR_AUTOSITE_PATH} -r {quote(str(pdbqt_path))} -n 50'
            proc = call(cmd, cwd=workdir, logger=self.logger)
            if output_dir.exists() and output_dir.joinpath('target_summary.csv').exists():
                self.logger.info(f'autosite for {self.whole_pdb} success.')
            else:
                self.logger.info(f'autosite for {self.whole_pdb} failed.')
                self.logger.debug(f'STDERR: {proc["stderr"]}')
                self.logger.debug(f'STDOUT: {proc["stdout"]}')

            # pocket will in form of
            #   `target/target_fp_<id:03d>.pdb` to show neighbor atoms
            #   `target/target_cl_<id:03d>.pdb` to show neighbor atoms
            result = []
            for idx in itertools.count(1):
                cluster_path = output_dir.joinpath(f'target_cl_{idx:03d}.pdb')
                if not cluster_path.exists():
                    self.logger.info(f'#{idx - 1} pocket(s) extracted')
                    break

                pdbfile = PDBFile.cast(cluster_path.read_text())
                bounding = Selection.from_luck('system', pdbfile=pdbfile)

                # autosite will not report pocket volume, so here we use OBB volume as it upper limit
                # TODO: maybe we should use better criteria
                obb_dimension = bounding.bounding.OBB.dimension
                volume_upper_limit = obb_dimension[0] * obb_dimension[1] * obb_dimension[2]
                vdw_volume = sorted(vdw_volume)
                if not (vdw_volume[0] <= volume_upper_limit <= vdw_volume[1]):
                    self.logger.info(f'volume of pocket {idx} (<={volume_upper_limit}) out of range of {vdw_volume}, discarded.')
                    break

                result.append((bounding, PocketInfo()))
            else:
                self.logger.info(f'#{idx} pocket(s) extracted')
            return result[:number]


def gen_dock_site_smart(receptor: Path, site_num: int, ligand_identification: Optional[str] = None, *args, engine: str, pocket_volume_lb: float = 0, logger: logging.Logger = logging.getLogger('vina_engine.vina123.gen_dock_site_smart')) -> List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """ automatic dock site pick

    Args:
        receptor (Path): receptor path (in format of pdb)
        site_num (int): maxima site
        engine (str): picker engine
        ligand_identification (Optional[str], optional): ligand identification to specify pocket site. Defaults to None.
        pocket_volume_lb (float, optional): _description_. Defaults to 0.
        logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('vina_engine.vina123.gen_dock_site_smart').

    Raises:
        NotImplementedError: pocker finder not implemented

    Returns:
        List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]: list of AABB box about center and dimension
    """
    if engine.lower() == 'fpocket':
        from .pocket import FPocket
        picker = FPocket(receptor, ligand_identification=ligand_identification, logger=logger)
    elif engine.lower() == 'autosite':
        from .pocket import AutoSite
        picker = AutoSite(receptor, ligand_identification=ligand_identification, logger=logger)
    else:
        raise NotImplementedError(f'request pocket finder {engine} not implemented.')
    try:
        sites = picker.search(site_num, vdw_volume=(pocket_volume_lb, 1.e9))
    except RuntimeError as e:
        logger.error(f'fail to generate site, {e}', exc_info=e)
        raise ValueError(str(e))
    boxes = [(sele.bounding.AABB.center, sele.bounding.AABB.dimension) for sele, _ in sites]
    return boxes
