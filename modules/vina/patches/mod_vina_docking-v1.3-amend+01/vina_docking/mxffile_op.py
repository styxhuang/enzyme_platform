import logging
from pathlib import Path
from typing import Generator, List, Tuple, overload

from mxf_files.files import MXFFile, PDBFile
from openbabel import openbabel

from .exceptions import InputStructureNotGoodException


@overload
def split_conformers(target: Path, store_at: None, logger: logging.Logger = logging.getLogger('mxffile_op.split_conformers')) -> Tuple[List[bytes], List[bytes]]:
    ...


@overload
def split_conformers(target: Path, store_at: Path, logger: logging.Logger = logging.getLogger('mxffile_op.split_conformers')) -> Tuple[List[Path], List[Path]]:
    ...


def split_conformers(target, store_at, logger=logging.getLogger('mxffile_op.split_conformers')):
    """ spliter of mxf files about single frame or multiple frame into receptors or ligands

    Args:
        target (Path): mxf file path
        store_at (Path, None): where to store the file, if None, data will be in memory
        logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('mxffile_op.split_conformers').

    Raises:
        InputStructureNotGoodException: target path is not a mxf suffix path

    Returns:
        Tuple[List[bytes], List[bytes]] | Tuple[List[Path], List[Path]]: mxf file content if `store_at` is None, else mxf file path of receptors and ligands
    """
    if target.suffix != '.mxf':
        logger.error(f'provided conformers is not in format of .mxf ({target})')
        raise InputStructureNotGoodException(f'file format not good, need .mxf, get {target.suffix}')

    try:
        mxffile = MXFFile(target)
    except Exception as e:
        logger.error('mxf file is malformed or mxf toolkit version mismatch.', exc_info=e)
        raise e

    # first triage single frame or multiple frame, spliting multiple frame into list of single frame
    if hasattr(mxffile.fp, 'StructureType') and (mxffile.fp.StructureType == 'MultiStructure'):
        mxf_buf = []
        for submol_idx, submol in enumerate(mxffile.fp.Structures):
            plain_mxf = MXFFile()
            plain_mxf.fp = submol
            if isinstance(store_at, Path):
                plain_mxf.fn = store_at.joinpath(Path(mxffile.fn).name).with_suffix(f'.{submol_idx}.mxf')
                plain_mxf.flush()
            else:
                raise NotImplementedError
            mxf_buf.append(plain_mxf)
    else:
        mxf_buf = [mxffile, ]
    mxffiles: Generator[MXFFile, None, None] = (i for i in mxf_buf)

    # then check molecule weight to determine receptor (MW > 3000) or ligand
    receptors = []
    ligands = []
    for single_mxf in mxffiles:
        pdb_content = str(PDBFile.cast(single_mxf))

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('pdb', 'pdb')
        conv.ReadString(mol, pdb_content)
        if mol.GetMolWt() > 3000:
            receptors.append(Path(single_mxf.fn))
        else:
            ligands.append(Path(single_mxf.fn))

    return receptors, ligands
