class GeneralVinaException(Exception):
    pass


class InputStructureNotGoodException(GeneralVinaException):
    msg = 'Input structure is not good (General Error)'


class MGLToolsPrepGPF4GeneralException(GeneralVinaException):
    msg = 'prepare_gpf4 script in MGLTools fail to generate gpf file (General Error)'


class MGLToolsPrepLigand4GeneralException(GeneralVinaException):
    msg = 'prepare_ligand4 script in MGLTools fail to generate pdbqt file (General Error)'


class MGLToolsPrepReceptor4GeneralException(GeneralVinaException):
    msg = 'prepare_receptor4 script in MGLTools fail to generate pdbqt file (General Error)'


class MGLToolsProcessVinaResultGeneralException(GeneralVinaException):
    msg = 'process_VinaResult script failed, ligand is likely to be unreasonable'


class GenerateEmptyFlexibleResidues(MGLToolsPrepReceptor4GeneralException):
    msg = 'prepare_receptor4 generate empty flexible pdbqt file'


class ADFRAutoGridGeneralException(GeneralVinaException):
    msg = 'autogrid failed (General Error)'


class VinaExmapleScriptGeneralException(GeneralVinaException):
    msg = 'AutoDock vina example script failed (General Error)'


class VinaAtomTypeNotSupport(GeneralVinaException):
    msg = 'AutoDock vina cannot recognize specific atom type'


class VagueLigandBond(GeneralVinaException):
    msg = 'No bonding infomation found in provided ligand, forget to preprocess it?'


class RemoteProcedureCallFailed(GeneralVinaException):
    msg = 'Calculation submit failed (Sugon/Changchun/...'
