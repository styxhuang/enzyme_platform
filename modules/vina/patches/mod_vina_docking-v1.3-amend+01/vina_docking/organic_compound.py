''' Organic compound gadget module
'''
import logging
import tempfile
from enum import Enum
from pathlib import Path
from typing import Callable, List
from warnings import warn

try:
    from openbabel import openbabel
except (ImportError, ModuleNotFoundError):
    _OPENBABEL_ENABLE = False
else:
    _OPENBABEL_ENABLE = True

# derivatived from 10.1021/J100785A001
VDW_RADII_VOLUME_BONDI = dict(
    H=7.24,
    C=20.58,
    N=15.6,
    O=14.71,
    F=13.31,
    Cl=22.45,
    Br=26.52,
    I=32.52,
    P=24.43,
    S=24.43,
    As=26.52,
    B=40.48,
    Si=38.79,
    Se=28.73,
    Te=36.62,
)

# derivatived from Cordero et al., Dalton Trans. 2832-2838, 2008
COVALENT_RADII = dict(
    H=0.31,
    He=0.28,
    Li=1.28,
    Be=0.96,
    B=0.84,
    C_sp3=0.76,
    C_sp2=0.73,
    C_sp=0.69,
    N=0.71,
    O=0.66,
    F=0.57,
    Ne=0.58,
    Na=1.66,
    Mg=1.41,
    Al=1.21,
    Si=1.11,
    P=1.07,
    S=1.05,
    Cl=1.02,
    Ar=1.06,
    K=2.03,
    Ca=1.76,
    Sc=1.7,
    Ti=1.6,
    V=1.53,
    Cr=1.39,
    Mn_ls=1.39,
    Mn_hs=1.61,
    Fe_ls=1.32,
    Fe_hs=1.52,
    Co_ls=1.26,
    Co_hs=1.5,
    Ni=1.24,
    Cu=1.32,
    Zn=1.22,
    Ga=1.22,
    Ge=1.2,
    As=1.19,
    Se=1.2,
    Br=1.2,
    Kr=1.16,
    Rb=2.2,
    Sr=1.95,
    Y=1.9,
    Zr=1.75,
    Nb=1.64,
    Mo=1.54,
    Tc=1.47,
    Ru=1.46,
    Rh=1.42,
    Pd=1.39,
    Ag=1.45,
    Cd=1.44,
    In=1.42,
    Sn=1.39,
    Sb=1.39,
    Te=1.38,
    I=1.39,
    Xe=1.4,
    Cs=2.44,
    Ba=2.15,
    La=2.07,
    Ce=2.04,
    Pr=2.03,
    Nd=2.01,
    Pm=1.99,
    Sm=1.98,
    Eu=1.98,
    Gd=1.96,
    Tb=1.94,
    Dy=1.92,
    Ho=1.92,
    Er=1.89,
    Tm=1.9,
    Yb=1.87,
    Lu=1.87,
    Hf=1.75,
    Ta=1.7,
    W=1.62,
    Re=1.51,
    Os=1.44,
    Ir=1.41,
    Pt=1.36,
    Au=1.36,
    Hg=1.32,
    Tl=1.45,
    Pb=1.46,
    Bi=1.48,
    Po=1.4,
    At=1.5,
    Rn=1.5,
    Fr=2.6,
    Ra=2.21,
    Ac=2.15,
    Th=2.06,
    Pa=2.00,
    U=1.96,
    Np=1.9,
    Pu=1.87,
    Am=1.8,
    Cm=1.69,
)


class vdwVolumeEstimator(Enum):
    VABC = 'vabc'


def vdw_volume(mol_fn: Path, flavor: vdwVolumeEstimator = vdwVolumeEstimator.VABC) -> float:
    if flavor == vdwVolumeEstimator.VABC:
        return _vdw_volume_vabc(mol_fn)


def _vdw_volume_vabc(fn: Path) -> float:
    """ implement of VABC issued at 10.1021/jo034808o

    Args:
        fn (Path): file path to a molecule file

    Returns:
        float: vdw volume
    """
    if not _OPENBABEL_ENABLE:
        raise NotImplementedError('implement of VABC depend on openbabel')

    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(fn.suffix.replace('.', ''), 'xyz')
    conv.ReadFile(mol, str(fn))
    n_aromatic_ring = 0
    n_nonaromatic_ring = 0
    for ring in mol.GetSSSR():
        if ring.IsAromatic():
            n_aromatic_ring += 1
        else:
            n_nonaromatic_ring += 1
    volume = 0
    n_atoms = 0
    for atom in openbabel.OBMolAtomIter(mol):
        symbol = openbabel.GetSymbol(atom.GetAtomicNum())
        atom_volume = VDW_RADII_VOLUME_BONDI.get(symbol, 0)
        if atom_volume == 0:
            warn(f'element {symbol} dont have Bondi vdw volume, fallback to 0', UserWarning)
        volume += VDW_RADII_VOLUME_BONDI.get(symbol, 0)
        n_atoms += 1

    volume = volume - 5.92 * (n_atoms - 1 + n_aromatic_ring + n_nonaromatic_ring) - 14.7 * n_aromatic_ring - 3.8 * n_nonaromatic_ring
    return volume


def _get_smiles_via_openbabel(file_path: Path) -> List[str]:
    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(file_path.suffix.lstrip('.'), 'smi')
    if conv.ReadString(mol, file_path.read_text()):
        return conv.WriteString(mol).splitlines()
    else:
        return list()


def get_smiles_from_file(file_path: Path, on_fail: Callable[[Path], str] = lambda x: 'NULL') -> List[str]:
    """ generate SMILES for given file

    Args:
        file_path (Path): file path
        on_fail (Callable[[Path], str], optional): string that on fail. Defaults to lambdax:'NULL'.

    Raises:
        NotImplementedError: openbabel is not installed

    Returns:
        List[str]: list of SMILES
    """
    if _OPENBABEL_ENABLE:
        result = _get_smiles_via_openbabel(file_path)
        if result:
            return result
        else:
            return [on_fail(file_path)]

    raise NotImplementedError('Openbabel is not installed, SMILES cannot be converted')


def prepare_ligand(ligand: Path, *args, logger: logging.Logger = logging.getLogger('organic_compound.prepare_ligand')) -> str:
    """ prepare ligand via MGLTools and wet.py in Vina example script

    Args:
        ligand (Path): Ligand pdb file path
        hydrated (bool): wether to make ligand wet
        logger (logging.Logger, optional): logging facility. Defaults to logging.getLogger('organic_compound.prepare_ligand').

    Raises:
        ValueError: conversion failed

    Returns:
        str: content of prepared pdbqt file
    """
    with tempfile.TemporaryDirectory() as workdir:
        ligand_path = Path(workdir).joinpath('ligand').with_suffix(ligand.suffix)
        ligand_path.symlink_to(ligand)

        obconv = openbabel.OBConversion()
        obconv.SetInAndOutFormats('pdb', 'pdbqt')
        obmol = openbabel.OBMol()
        obconv.ReadFile(obmol, str(ligand_path))
        ligand_pdbqt_content = obconv.WriteString(obmol)

        ligand_pdbqt = ligand_path.with_suffix('.pdbqt')
        if ligand_pdbqt_content:
            ligand_pdbqt.write_text(ligand_pdbqt_content)
            logger.info(f'pdbqt file generated for ligand ({ligand})')
        else:
            logger.error(f'pdbqt file fail to generated for and ligand ({ligand})')
            raise ValueError('Openbabel fail to convert pdb to pdbqt')

        result = ligand_pdbqt.read_text()
        return result
