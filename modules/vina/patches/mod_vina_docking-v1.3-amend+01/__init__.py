'''entry point of application to be compatible with legacy back-end
'''

_ONLINE = False

try:
    if not _ONLINE:
        from utiles import MXF
except (ImportError, ModuleNotFoundError):
    _ONLINE = False
else:
    _ONLINE = True

try:
    if not _ONLINE:
        from new_pkg.pipelines.utiles import MXF
except (ImportError, ModuleNotFoundError, ValueError):
    _ONLINE = False
else:
    _ONLINE = True

try:
    if not _ONLINE:
        from ..utiles import MXF
except (ImportError, ModuleNotFoundError, ValueError):
    _ONLINE = False
else:
    _ONLINE = True

if _ONLINE:
    class MOD_CADD_VINA_DOCKING_1_3:
        @staticmethod
        def macros_ligand_docking(tag_id, *ID):
            from .vina_docking import MacrosLigandsDocking
            with MacrosLigandsDocking(MXF(tag_id)) as node:
                return node.run()

else:
    import warnings
    warnings.warn('fail to import utils in pipelines')

__all__ = ['MOD_CADD_VINA_DOCKING_1_3']
