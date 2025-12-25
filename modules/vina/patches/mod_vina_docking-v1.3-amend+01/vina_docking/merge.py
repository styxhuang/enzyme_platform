'''model merge module (caddmodules.workhorse.modelling.merge)
'''
import logging
from pathlib import Path
from enum import Enum
from warnings import warn
import tempfile

from mxf_files.files import PDBFile

try:
    import ost.io
    from ost.mol import INCLUDE_ALL
except (ImportError, ModuleNotFoundError):
    warn('openstructure failed to import, some functionality may compromised')
    _OPENSTRUCTURE_ENABLE = False
else:
    _OPENSTRUCTURE_ENABLE = True


CHAIN_CHARS = set([chr(i) for i in range(65, 91)])  # A-Z


class MergeFlavor(Enum):
    DEFAULT = SPARE = 0
    SHIFT = 1
    MACRO_LIGAND_DOCKING = 2


def pdb_merge(fix_pdb: str, float_pdb: str, flavor: MergeFlavor = MergeFlavor.DEFAULT, *, no_hetatms: bool = False, logger: logging.Logger = logging.getLogger('caddmodules.workhorse.modelling.merge.pdb_merge')) -> str:
    if flavor.value == MergeFlavor.SPARE.value:
        if _OPENSTRUCTURE_ENABLE:
            with tempfile.TemporaryDirectory() as workdir:
                fix = Path(workdir).joinpath('fix.pdb')
                flt = Path(workdir).joinpath('flt.pdb')
                fix.write_text(fix_pdb)
                flt.write_text(float_pdb)
                return pdb_merge_spare_via_ost(fix, flt, no_hetatms=no_hetatms, logger=logger.getChild('via_ost'))
        else:
            raise NotImplementedError(f'Missing package to process `{flavor}` flavor merge')
    elif flavor.value == MergeFlavor.SHIFT.value:
        if _OPENSTRUCTURE_ENABLE:
            with tempfile.TemporaryDirectory() as workdir:
                fix = Path(workdir).joinpath('fix.pdb')
                flt = Path(workdir).joinpath('flt.pdb')
                fix.write_text(fix_pdb)
                flt.write_text(float_pdb)
                return pdb_merge_shift_via_ost(fix, flt, no_hetatms=no_hetatms, logger=logger.getChild('via_ost'))
        else:
            raise NotImplementedError(f'Missing package to process `{flavor}` flavor merge')

    elif flavor.value == MergeFlavor.MACRO_LIGAND_DOCKING.value:
        return pdb_merge_macros_ligand_via_native(fix_pdb, float_pdb, logger=logger)

    else:
        raise NotImplementedError(f'`{flavor}` flavor is not implemented')


def pdb_merge_spare_via_ost(fix_pdb: Path, float_pdb: Path, *, no_hetatms: bool = False, logger: logging.Logger = logging.getLogger('caddmodules.workhorse.modelling.merge.pdb_merge_spare_via_ost')) -> str:
    """ Merge pdb file using spare chain id via OpenStructure

    Args:
        fix_pdb (Path): path to pdb file that chain name should be fixed
        float_pdb (Path): path to pdb file that chain name can be floatting
        no_hetatms (bool, optional): keep HETATM record. Defaults to False.
        logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('caddmodules.workhorse.modelling.merge.pdb_merge_spare_via_ost').

    Returns:
        str: merged pdb content
    """
    pdb_mol1 = ost.io.LoadPDB(str(fix_pdb), no_hetatms=no_hetatms, fault_tolerant=True)
    pdb_mol2 = ost.io.LoadPDB(str(float_pdb), no_hetatms=no_hetatms, fault_tolerant=True)

    mol1_chains = set([chainhandler.name for chainhandler in pdb_mol1.chains])
    mol2_chains = set([chainhandler.name for chainhandler in pdb_mol2.chains])
    overlap_chains = mol1_chains & mol2_chains
    logger.debug(f'chain(s) from fixed pdb are {mol1_chains}; from float pdb are {mol2_chains}')
    if len(overlap_chains) == 0:
        logger.debug('no overlap about chain name, merge directly')
    elif len((mol1_chains | mol2_chains) ^ CHAIN_CHARS) < len(overlap_chains):
        logger.warning(f'insufficient spare chain id for relocation of overlapped chain {overlap_chains}, merge directly')
    else:
        logger.info(f'relocating chain {overlap_chains} from protein2')
        spare_chains = list((mol1_chains | mol2_chains) ^ CHAIN_CHARS)
        if ' ' in spare_chains:
            # we dont want to relocate as empty
            spare_chains.pop(spare_chains.index(' '))
        spare_chains = sorted(spare_chains)
        mol2_editor = pdb_mol2.EditXCS()
        for overlap_chain, spare_chain in zip(overlap_chains, spare_chains):
            logger.info(f'chain {overlap_chain} is relocated to {spare_chain}')
            mol2_editor.RenameChain(pdb_mol2.FindChain(overlap_chain), spare_chain)

    mol1_view = pdb_mol1.CreateFullView()
    for chain in pdb_mol2.chains:
        mol1_view.AddChain(chain, INCLUDE_ALL)

    return ost.io.EntityToPDBStr(mol1_view)


def pdb_merge_shift_via_ost(fix_pdb: Path, float_pdb: Path, *, no_hetatms: bool = False, logger: logging.Logger = logging.getLogger('caddmodules.workhorse.modelling.merge.pdb_merge_spare_via_ost')) -> str:
    """ Merge pdb file using next available chain id via OpenStructure

    Args:
        fix_pdb (Path): path to pdb file that chain name should be fixed
        float_pdb (Path): path to pdb file that chain name can be floatting
        no_hetatms (bool, optional): keep HETATM record. Defaults to False.
        logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('caddmodules.workhorse.modelling.merge.pdb_merge_spare_via_ost').

    Returns:
        str: merged pdb content
    """
    pdb_mol1 = ost.io.LoadPDB(str(fix_pdb), no_hetatms=no_hetatms, fault_tolerant=True)
    pdb_mol2 = ost.io.LoadPDB(str(float_pdb), no_hetatms=no_hetatms, fault_tolerant=True)

    mol1_chains = set([chainhandler.name for chainhandler in pdb_mol1.chains])
    mol2_chains = set([chainhandler.name for chainhandler in pdb_mol2.chains])
    overlap_chains = mol1_chains & mol2_chains
    logger.debug(f'chain(s) from fixed pdb are {mol1_chains}; from float pdb are {mol2_chains}')
    if len(overlap_chains) == 0:
        logger.debug('no overlap about chain name, merge directly')
    elif len((mol1_chains | mol2_chains) ^ CHAIN_CHARS) < len(overlap_chains):
        logger.warning(f'insufficient spare chain id for relocation of overlapped chain {overlap_chains}, merge directly')
    else:
        logger.info(f'relocating chain {overlap_chains} from protein2')
        spare_chains = list((mol1_chains | mol2_chains) ^ CHAIN_CHARS)
        if ' ' in spare_chains:
            # we dont want to relocate as empty
            spare_chains.pop(spare_chains.index(' '))
        spare_chains = sorted(spare_chains)  # FIXME
        mol2_editor = pdb_mol2.EditXCS()
        for overlap_chain, spare_chain in zip(overlap_chains, spare_chains):
            logger.info(f'chain {overlap_chain} is relocated to {spare_chain}')
            mol2_editor.RenameChain(pdb_mol2.FindChain(overlap_chain), spare_chain)

    mol1_view = pdb_mol1.CreateFullView()
    for chain in pdb_mol2.chains:
        mol1_view.AddChain(chain, INCLUDE_ALL)

    return ost.io.EntityToPDBStr(mol1_view)


def pdb_merge_macros_ligand_via_native(macro_pdb_content: str, ligand_pdb_content: str, *, logger: logging.Logger = logging.getLogger('caddmodules.workhorse.modelling.merge.pdb_merge_macros_ligand_via_native')) -> str:
    macros = PDBFile.cast(macro_pdb_content)
    ligand = PDBFile.cast(ligand_pdb_content)

    macros_chains = set([atom.nameChain for atom in macros.atoms])
    ligand_chain = set([atom.nameChain for atom in ligand.atoms])

    if len(ligand_chain) > 1:
        logger.warning(f'ligand has multiple chain ({ligand_chain}), it will be squeezed.')
    overlap_chains = macros_chains & ligand_chain
    if len(overlap_chains) == 0:
        logger.debug('no overlap about chain name, merge directly')
    elif len((macros_chains | ligand_chain) ^ CHAIN_CHARS) < len(overlap_chains):
        logger.warning(f'insufficient spare chain id for relocation of overlapped chain {overlap_chains}, merge directly')
    else:
        logger.info(f'relocating chain {overlap_chains} from ligand')
        spare_chains = list((macros_chains | ligand_chain) ^ CHAIN_CHARS)
        if ' ' in spare_chains:
            # we dont want to relocate as empty
            spare_chains.pop(spare_chains.index(' '))
        ligand_chain = {sorted(spare_chains)[0]}
    ligand_chain = list(ligand_chain)[0]

    buf = ''
    atom_idx = 1
    rec_connect_map = dict()
    lig_connect_map = dict()
    for rec_atom in macros.atoms:
        rec_connect_map[rec_atom.idxAtom] = atom_idx
        rec_atom.idxAtom = atom_idx
        buf += f'{rec_atom.to_pdb()}\n'
        atom_idx += 1
    buf += 'TER\n'
    for lig_atom in ligand.atoms:
        lig_connect_map[lig_atom.idxAtom] = atom_idx
        lig_atom.idxAtom = atom_idx
        lig_atom.nameChain = ligand_chain
        buf += f'{lig_atom.to_pdb()}\n'
        atom_idx += 1

    # 4. updating connect record
    for rec_conn in macros.connects:
        rec_conn.atom1 = rec_connect_map.get(rec_conn.atom1, '')
        rec_conn.atom2 = rec_connect_map.get(rec_conn.atom2, '')
        rec_conn.atom3 = rec_connect_map.get(rec_conn.atom3, '')
        rec_conn.atom4 = rec_connect_map.get(rec_conn.atom4, '')
        rec_conn.atom5 = rec_connect_map.get(rec_conn.atom5, '')
        buf += f'{rec_conn.to_pdb()}\n'
    for lig_conn in ligand.connects:
        lig_conn.atom1 = lig_connect_map.get(lig_conn.atom1, '')
        lig_conn.atom2 = lig_connect_map.get(lig_conn.atom2, '')
        lig_conn.atom3 = lig_connect_map.get(lig_conn.atom3, '')
        lig_conn.atom4 = lig_connect_map.get(lig_conn.atom4, '')
        lig_conn.atom5 = lig_connect_map.get(lig_conn.atom5, '')
        buf += f'{lig_conn.to_pdb()}\n'

    buf += 'END\n'
    return buf
