from unittest import TestCase


class TestImport(TestCase):
    def test_import_vina_engine(self):
        from vina_docking._vina_engine import COLLISION_SCOPE  # noqa
        from vina_docking._vina_engine import detect as collision_detect  # noqa
        from vina_docking._vina_engine import restore_ligand as _restore_ligand  # noqa
