'''Node for Docking (Macros & Ligands) (ID: LS-02-000-0003)
'''
import itertools
from pathlib import Path
from typing import Generator, Iterable, List, Tuple

import pandas as pd
from mxf_files.files import MXFFile, PDBFile
from mxf_files.selection import MaxFlowSetType, Selection
from mxf_serum import CanRunRemote, Serum
from mxf_serum.cds.structure_data_pipeline import SimsFeed, SimsStructureType
from mxf_serum.dp.template_method import StaticBufferedAsyncSimsPattern
from mxf_serum.ormp.base import ParameterDelegation
from typing_extensions import TypedDict

from .safe_fn import safe_fn
from .ui import MacrosLigandsDockingParameter


class TriageResult(TypedDict):
    receptor: List[Path]    # path to the triaged receptor pdb
    ligand: List[Path]      # path to the triaged ligand pdb
    errno: int              # error number
                            # 0 - safe and sound
                            # 1 - no .mxf file found


class LigandInfo(TypedDict):
    name: str        # ligand name
    inlet_pdb: Path  # provided pdb file path
    pdbqt: Path      # prepared pdbqt file
    volume: float    # molecule volume
    smiles: str      # SMILES represent
    _deprecated: str


class ReceptorInfo(TypedDict):
    name: str
    inlet_mxf: Path
    inlet_pdb: Path
    rigid_pdbqt: Path
    flex_pdbqt: Path
    center: Tuple[float, float, float]
    dimension: Tuple[float, float, float]
    site_idx: int
    _deprecated: str


class MacrosLigandsDocking(Serum, CanRunRemote, ParameterDelegation, StaticBufferedAsyncSimsPattern):
    parameter_definition = MacrosLigandsDockingParameter

    channel_muxer = 'squash'
    concurrent = 50
    task_chunk_size = 500

    _triage_dir = '.splited_structure'
    _triage_counter = 0

    _ligand_dir = '.ligands'
    _ligand_counter = 0

    _receptor_dir = '.receptors'
    _receptor_counter = 0

    @property
    def triage_dir(self) -> Path:
        self._triage_counter += 1

        if isinstance(self._triage_dir, Path):
            return self._triage_dir.joinpath(str(self._triage_counter))
        self._triage_dir = self.mxfin.node_path.joinpath(self._triage_dir)
        return self._triage_dir.joinpath(str(self._triage_counter))

    def triage_mols(self, inlet_feed: SimsFeed) -> TriageResult:
        from .mxffile_op import split_conformers
        self.logger.info(f'split incoming {inlet_feed.name} (mxf files: {inlet_feed.mxf})')

        if len(inlet_feed.mxf) == 0:
            self.logger.error(f'fail to process {inlet_feed.name} as it doesnt having any mxf files')
            return TriageResult(receptor=[], ligand=[], errno=1)

        if len(inlet_feed.mxf) > 1:
            self.logger.warning(f'structure {inlet_feed.name} having multiple mxf file, only first one is used.')
        inlet_mxf = inlet_feed.mxf[0]

        raw_structure_store_at = self.triage_dir
        raw_structure_store_at.mkdir(parents=True, exist_ok=True)

        result = TriageResult(receptor=[], ligand=[], errno=0)
        split_receptor, split_ligand = split_conformers(inlet_mxf.decode(), store_at=raw_structure_store_at, logger=self.logger.getChild('split_conformers'))
        processing_idx = 0
        for split_mxf_path in split_receptor:
            converted_pdb = PDBFile.cast(MXFFile(split_mxf_path))
            converted_pdb_path = raw_structure_store_at.joinpath(f'{processing_idx:05d}.pdb')
            converted_pdb_path.write_text(str(converted_pdb))
            result['receptor'].append(converted_pdb_path)
            processing_idx += 1
        for split_mxf_path in split_ligand:
            converted_pdb = PDBFile.cast(MXFFile(split_mxf_path))
            converted_pdb_path = raw_structure_store_at.joinpath(f'{processing_idx:05d}.pdb')
            converted_pdb_path.write_text(str(converted_pdb))
            result['ligand'].append(converted_pdb_path)
            processing_idx += 1

        self.logger.info(f'#{len(split_receptor)} receptor(s), #{len(split_ligand)} ligand(s) are triaged for {inlet_feed.name}.')
        return result

    def prepare_ligands(self, ligand: LigandInfo) -> LigandInfo:
        from .organic_compound import (get_smiles_from_file, prepare_ligand,
                                       vdw_volume)

        # store SMILES of ligand
        try:
            ligand_smi = get_smiles_from_file(ligand['inlet_pdb'])[0]
        except Exception as e:
            self.logger.warning(f'fail to generate SMILES for `{ligand["inlet_pdb"]}`', exc_info=e)
            ligand_smi = 'Null'

        try:
            prep_pdbqt = prepare_ligand(ligand['inlet_pdb'], logger=self.logger.getChild('prep_ligand'))
        except ValueError:
            # possible empty ligand, deprecate this ligand
            result = ligand.copy()
            result['_deprecated'] = 'Invalid structure file, fail to convert as PDBQT file.'
            return result

        prep_pdbqt_path = self.mxfin.node_path.joinpath(self._ligand_dir, f'{self._ligand_counter:09d}.pdbqt')
        prep_pdbqt_path.parent.mkdir(exist_ok=True)
        prep_pdbqt_path.write_text(prep_pdbqt)
        self._ligand_counter += 1
        result = ligand.copy()
        result['pdbqt'] = prep_pdbqt_path
        result['volume'] = vdw_volume(ligand['inlet_pdb'])
        result['smiles'] = ligand_smi
        return result

    def prepare_docking_site(self, receptor: ReceptorInfo, max_ligand: float) -> List[ReceptorInfo]:
        from .pocket import gen_dock_site_smart

        if self.parameter.domain_type == 1:
            # use automatic docking site
            if self.parameter.auto_domain_lb == 0:
                max_ligand_volume = max_ligand
            else:
                max_ligand_volume = self.parameter.auto_domain_lb

            aabbs = gen_dock_site_smart(
                receptor['inlet_pdb'],
                self.parameter.auto_domain_num,
                engine=self.parameter.auto_domain_method,
                pocket_volume_lb=max_ligand_volume,
                logger=self.logger.getChild('gen_docksite'),
            )

            if len(aabbs) == 0:
                return list()

            buf = []
            for site_idx, (center, dimension) in enumerate(aabbs):
                payload = receptor.copy()
                payload['center'] = center
                payload['dimension'] = dimension
                payload['site_idx'] = site_idx + 1
                buf.append(payload)
                self.logger.info(f'recpetor {payload["inlet_pdb"]} with docking site at {center} dimension {dimension} prepared')
            return buf

        elif self.parameter.domain_type == 2:
            # use box select UI
            try:
                mxf_sets_lists = self.parameter.box_select([receptor['name'], ], MXFFile(receptor['inlet_mxf']))
            except (ValueError, NotImplementedError):
                # invalid or no group select would break here
                payload = receptor.copy()
                payload['_deprecated'] = 'No docking site set for this receptor'
                return [payload, ]
            # unrolling sets
            mxf_set_lists = itertools.chain(*[sets_list.sets for sets_list in mxf_sets_lists])
            # filter bounding box type set
            mxf_bb_sets_lists = filter(lambda x: x.type == MaxFlowSetType.CuboidBounding, mxf_set_lists)
            mxf_bb_sets_list: List[Selection] = itertools.chain(*[bb_sets.atoms for bb_sets in mxf_bb_sets_lists])
            buf = []
            for site_idx, mxf_bb in enumerate(mxf_bb_sets_list):
                payload = receptor.copy()
                payload['center'] = mxf_bb.bounding.AABB.center
                payload['dimension'] = mxf_bb.bounding.AABB.dimension
                payload['site_idx'] = site_idx + 1
                buf.append(payload)
                self.logger.info(f'recpetor {receptor["name"]} with docking site at {buf[-1]["center"]} dimension {buf[-1]["dimension"]} prepared')
            return buf
        else:
            raise NotImplementedError(f'Docking site method {self.parameter.domain_type} is not implemented')

    def prepare_receptor(self, receptor: ReceptorInfo) -> ReceptorInfo:
        from .exceptions import GeneralVinaException
        from .vina_engine import Vina123

        try:
            if self.parameter.flex_docking:
                prep_pdbqt_rigid, prep_pdbqt_flex = Vina123.prepare_receptor(
                    receptor['inlet_pdb'],
                    flex=[self.parameter.flex_res, ],
                    logger=self.logger.getChild('prep_receptor'),
                )
            else:
                prep_pdbqt_rigid, prep_pdbqt_flex = Vina123.prepare_receptor(
                    receptor['inlet_pdb'],
                    flex=None,
                    logger=self.logger.getChild('prep_receptor'),
                )

        except GeneralVinaException as e:
            self.logger.error(f'fail to prepare {receptor["name"]}, due to {e}')
            raise

        prep_pdbqt_rigid_path = self.mxfin.node_path.joinpath(self._receptor_dir, f'{self._receptor_counter:06d}.rigid.pdbqt')
        prep_pdbqt_flex_path = self.mxfin.node_path.joinpath(self._receptor_dir, f'{self._receptor_counter:06d}.flex.pdbqt')
        prep_pdbqt_rigid_path.parent.mkdir(parents=True, exist_ok=True)
        prep_pdbqt_rigid_path.write_text(prep_pdbqt_rigid)
        receptor['rigid_pdbqt'] = prep_pdbqt_rigid_path
        if prep_pdbqt_flex:
            prep_pdbqt_flex_path.write_text(prep_pdbqt_flex)
            receptor['flex_pdbqt'] = prep_pdbqt_flex_path
        else:
            prep_pdbqt_flex_path = None

        return receptor

    def on_deliver(self, feed_gen: Iterable[Tuple[int, SimsFeed]]) -> Generator[Tuple[int, Tuple[ReceptorInfo, LigandInfo, Path, Path]], None, None]:
        from .exceptions import GeneralVinaException

        self.logger.info('Start to analyze receptors and ligands structure file')

        pending_receptors = []
        pending_ligands = []
        for feed_idx, feed in feed_gen:
            triage = self.triage_mols(feed)
            if triage['errno'] != 0:
                self.logger.warning(f'fail to triage {feed.name}')
                continue

            for pdb_path in triage['receptor']:
                pending_receptors.append(
                    ReceptorInfo(
                        name=feed.name,
                        inlet_mxf=feed.mxf[0].decode(),
                        inlet_pdb=pdb_path,
                    )
                )

            for pdb_path in triage['ligand']:
                pending_ligands.append(
                    LigandInfo(
                        name=feed.name,
                        inlet_pdb=pdb_path
                    )
                )

        prepared_ligands = []
        for pending in pending_ligands:
            try:
                prepared = self.prepare_ligands(pending)
            except Exception as e:
                self.logger.warning(f'fail to prepare ligand {pending["name"]} due to {e}')
            prepared_ligands.append(prepared)

        max_ligand_volume = max([i['volume'] for i in prepared_ligands])
        self.logger.info(f'max ligand volume is {max_ligand_volume:.3f}')

        prepared_receptors = []
        for pending in pending_receptors:
            sites = self.prepare_docking_site(pending, max_ligand_volume)
            if self.parameter.domain_type == 1:
                sites = sites[:self.parameter.auto_domain_num]

            site_added = 0
            for site in sites:
                if site.get('_deprecated', False):
                    self.logger.info('a receptor with no docking site defined has been detected')
                    prepared_receptors.append(site)
                else:
                    try:
                        prepared = self.prepare_receptor(site)
                    except GeneralVinaException as e:
                        site['_deprecated'] = str(e)
                        self.logger.warning('a receptor was failed to prepare, deprecated.')
                        prepared_receptors.append(site)
                    else:
                        prepared_receptors.append(prepared)
                site_added += 1

            # patch up site number showing schema
            # only one site -> not showing site number
            # multiple site -> showing site number
            if site_added == 1:
                prepared_receptors[-1]['site_idx'] = -1

        self.logger.info(f'SUMMARY: #{len(prepared_ligands)} ligands and #{len(prepared_receptors)} receptors are prepared')
        self.task_total_count = len(prepared_ligands) * len(prepared_receptors)

        intermediate_files_folder = self.mxfin.node_path.joinpath('intermediate_files')
        results_files_folder = self.mxfin.node_path.joinpath('result_files')
        job_name_pool = set()
        for job_idx, (receptor, ligand) in enumerate(itertools.product(prepared_receptors, prepared_ligands)):
            partition_idx = job_idx // self.task_chunk_size
            intermediate_dir = intermediate_files_folder.joinpath(f'partition_{partition_idx:07d}', )
            result_dir = results_files_folder.joinpath(f'partition_{partition_idx:07d}')
            intermediate_dir.mkdir(exist_ok=True, parents=True)
            result_dir.mkdir(exist_ok=True, parents=True)

            # locate precise and non-collision folder
            if receptor['site_idx'] > 0:
                site_info_str = f'_Site{receptor["site_idx"]}'
            else:
                site_info_str = ''
            for loop in range(100):
                job_name = f'{safe_fn(receptor["name"], "Receptor")}_{safe_fn(ligand["name"], "Ligand")}{site_info_str}'
                if job_name not in job_name_pool:
                    break
            else:
                raise RuntimeError(f'try to avoid name collsion for {loop + 1} times, name collsion cannot be resolved.')

            intermediate_dir = intermediate_dir.joinpath(job_name)
            result_dir = result_dir.joinpath(job_name)
            # This dir should not collision
            intermediate_dir.mkdir()
            result_dir.mkdir()

            yield job_idx, (receptor, ligand, intermediate_dir, result_dir)

    def on_fail(self, receptor: ReceptorInfo, ligand: LigandInfo, workdir: Path, cleaned_dir: Path, *args, exception: Exception, **kwargs):
        return pd.DataFrame([
            {
                'structure_name': receptor['name'] + ' - ' + ligand['name'],
                'structure_path': '',
                'receptor_name': receptor['name'],
                'ligand_name': ligand['name'],
                'status': self.NodeStatus.FAIL.value,
                'error': str(exception)
            }
        ])

    async def partial_async_call(self, cmd: str, cwd: Path):
        return await self.remote_subprocess.async_call(
            cmd, cwd,
            tmpdir=True,
            identification=self.mxfin.cast_user_identity,
            logger=self.logger,
            job_conf={'cpu': 32, 'mem': 20480, 'app': 'Vina123'},
        )

    async def atomic_run(self, receptor: ReceptorInfo, ligand: LigandInfo, workdir: Path, cleaned_dir: Path, *args, process_idx: int, **kwargs):
        from .vina_engine import Vina123, Vina123Parms

        if receptor.get('_deprecated', False):
            # this receptor needs to be give up
            raise self.ModerateException(receptor['_deprecated'])

        if ligand.get('_deprecated', False):
            # this ligand needs to be give up
            raise self.ModerateException(ligand['_deprecated'])

        docking_name = workdir.name
        docking_parms = Vina123Parms(
            scoring=self.parameter.scoring_func,
            center_x=receptor['center'][0],
            center_y=receptor['center'][1],
            center_z=receptor['center'][2],
            size_x=receptor['dimension'][0],
            size_y=receptor['dimension'][1],
            size_z=receptor['dimension'][2],
            exhaustiveness=self.parameter.exhaustiveness,
            max_evals=self.parameter.max_evals,
            num_modes=self.parameter.n_poses,
            min_rmsd=self.parameter.min_rmsd,
            energy_range=self.parameter.energy_range,
        )

        try:
            if not self.parameter.flex_docking:
                proc, result = await Vina123.docking_simple(
                    self.partial_async_call,
                    receptor['rigid_pdbqt'],
                    ligand['pdbqt'],
                    receptor['inlet_pdb'],
                    ligand['inlet_pdb'],
                    workdir,
                    docking_parms,
                    logger=self.logger,
                )
            if self.parameter.flex_docking:
                proc, result = await Vina123.docking_flex(
                    self.partial_async_call,
                    receptor['rigid_pdbqt'],
                    receptor['flex_pdbqt'],
                    ligand['pdbqt'],
                    receptor['inlet_pdb'],
                    ligand['inlet_pdb'],
                    workdir,
                    docking_parms,
                    logger=self.logger,
                )
        except Exception as e:
            self.logger.error(f'fail to dock {receptor} & {ligand}', exc_info=e)
            return pd.DataFrame([{
                'structure_name': docking_name,
                'structure_path': str(workdir),
                'receptor_name': receptor['name'],
                'ligand_name': ligand['name'],
                'error': str(e),
                'docking site (center)': str(receptor['center']),
                'docking site (dimension)': str(receptor['dimension']),
                'status': self.NodeStatus.FAIL.value,
            }])
        else:
            self.logger.debug(f'Diagnosis of call: {proc}')
            self.logger.debug(f'Brief of result {result}')

        if len(result) == 0:
            return pd.DataFrame([{
                'structure_name': docking_name,
                'structure_path': str(workdir),
                'receptor_name': receptor['name'],
                'ligand_name': ligand['name'],
                'error': 'No docking mode is found.',
                'docking site (center)': str(receptor['center']),
                'docking site (dimension)': str(receptor['dimension']),
                'status': self.NodeStatus.FAIL.value,
            }])

        display_col = list(filter(lambda x: not x.startswith('_'), result.columns))
        docked_pdb = result['_com_fn'].to_list()
        docked_interaction = result['_interaction'].to_list()
        has_collision = result['_collision'].to_list()
        result['structure_name'] = [f'{docking_name}_{idx}' for idx, _ in enumerate(docked_pdb)]
        docked_mxf = []
        for docked_pdb_fn in docked_pdb:
            mxffp = MXFFile.cast(PDBFile(docked_pdb_fn), aux_mxffile=MXFFile(receptor['inlet_mxf']))
            mxffp.fn = str(Path(docked_pdb_fn).with_suffix('.mxf'))
            mxffp.flush()
            docked_mxf.append(str(Path(mxffp.fn).resolve()))
        result['structure_path'] = docked_mxf
        result['receptor_name'] = receptor['name']
        result['ligand_name'] = ligand['name']
        result['SMILES'] = ligand['smiles']
        result['docking site (center)'] = str(receptor['center'])
        result['docking site (dimension)'] = str(receptor['dimension'])
        result['multimedia'] = [
            [
                dict(
                    name=f'{docking_name}_{com_idx}',
                    type='PRO_NA_3D',
                    path=com_fn,
                ),
                dict(
                    name=f'2D interaction of {docking_name} Mode {com_idx}',
                    type='image',
                    type_species='INTERACTION_2D',
                    path=docked_interaction[com_idx],
                ),
                dict(
                    name=f'{docking_name} Mode {com_idx}',
                    type=SimsStructureType.PDB_FILE.value,
                    path=com_fn,
                ),
                dict(
                    name=f'{docking_name} Mode {com_idx}',
                    type=SimsStructureType.MDPDB_FILE.value,
                    path=com_fn,
                ),
                dict(
                    name=f'{docking_name} Mode {com_idx}',
                    type=SimsStructureType.MXF_FILE.value,
                    path=docked_mxf[com_idx],
                ),
            ]
            for com_idx, com_fn in enumerate(docked_pdb)
        ]
        result['status'] = [self.NodeStatus.FAIL.value if i else self.NodeStatus.SUCCESS.value for i in has_collision]
        result['error'] = ['Collision detected' if i else '' for i in has_collision]
        result = result[[
            'structure_name',
            'receptor_name',
            'ligand_name',
            'SMILES',
            'structure_path',
            *display_col,
            'docking site (center)',
            'docking site (dimension)',
            'multimedia',
            'status',
            'error',
        ]]

        # construct result directory
        output_filter = [True if i == self.NodeStatus.SUCCESS.value else False for i in result['status']]
        failed_dir = cleaned_dir.joinpath('fail')
        for com_idx, (com_fn, need_out) in enumerate(zip(docked_pdb, output_filter)):
            if need_out:
                user_friendly_fn = cleaned_dir.joinpath(f'{docking_name}_{com_idx}.pdb')
                try:
                    user_friendly_fn.symlink_to(com_fn)
                except Exception as e:
                    self.logger.warning(f'fail to symlink user friendly name: {e}')
                    continue

                # swap file registered for molstar display
                for media_idx, media in enumerate(result.iloc[com_idx]['multimedia']):
                    if media['type'] != 'PRO_NA_3D':
                        continue

                    result.iloc[com_idx]['multimedia'][media_idx]['path'] = str(user_friendly_fn)

            else:
                failed_dir.mkdir(exist_ok=True)
                failed_dir.joinpath(f'{docking_name}_{com_idx}.pdb').symlink_to(com_fn)

        return result
