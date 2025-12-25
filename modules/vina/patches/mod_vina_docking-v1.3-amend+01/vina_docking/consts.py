from pathlib import Path
from shlex import quote

_EXT_ROOT = Path(__file__).parent.joinpath('ext')

ADFR_AUTOSITE_PATH = f'PATH={quote(str(_EXT_ROOT / "ADFRsuite" / "bin"))}:$PATH autosite'
ADFR_AUTOGRID4_PATH = f'PATH={quote(str(_EXT_ROOT / "ADFRsuite" / "bin"))}:$PATH autogrid4'

_MGLTOOLS_PYTHONSH_PATH = str(_EXT_ROOT / 'MGLTools' / 'bin' / 'pythonsh')
_MGLTOOLS_ADT = _EXT_ROOT / 'MGLTools' / 'MGLToolsPckgs' / 'AutoDockTools'
_MGLTOOLS_ADT_UTILS24 = _MGLTOOLS_ADT / 'Utilities24'
MGLTOOLS_PREPARE_RECEPTOR4 = f'{quote(str(_MGLTOOLS_PYTHONSH_PATH))} {quote(str(_MGLTOOLS_ADT_UTILS24 / "prepare_receptor4.py"))}'
MGLTOOLS_PREPARE_GPF4 = f'{quote(str(_MGLTOOLS_PYTHONSH_PATH))} {quote(str(_MGLTOOLS_ADT_UTILS24 / "prepare_gpf4.py"))}'
MGLTOOLS_PREPARE_FLEX_RECEPTOR4 = f'{quote(str(_MGLTOOLS_PYTHONSH_PATH))} {quote(str(_MGLTOOLS_ADT_UTILS24 / "prepare_flexreceptor4.py"))}'
MGLTOOLS_PDBQT_TO_PDB = f'{quote(str(_MGLTOOLS_PYTHONSH_PATH))} {quote(str(_MGLTOOLS_ADT_UTILS24 / "pdbqt_to_pdb.py"))}'
MGLTOOLS_PROCESS_VINA_RESULT = f'{quote(str(_MGLTOOLS_PYTHONSH_PATH))} {quote(str(_MGLTOOLS_ADT_UTILS24 / "process_VinaResult.py"))}'
MGLTOOLS_AD4_PARMS_PATH = str(_MGLTOOLS_ADT / 'AD4_parameters.dat')

VINA_1_2_X_SEARCH_PATH = quote(str(_EXT_ROOT / 'vina' / '1.2.3'))

FPOCKET_PATH = quote(str(_EXT_ROOT / 'fpocket-4.0.3' / 'fpocket'))
