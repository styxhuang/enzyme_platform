import tempfile
from pathlib import Path
from unittest import TestCase

ASSETS = Path(__file__).parent.joinpath('assets')


class TestFlexDock(TestCase):
    def test_regression_of_bug45817(self):
        assets = ASSETS.joinpath('bug45817')

        pdbqt = assets / 'ligand.pdbqt'
        docked = assets / 'ligand_out_model1.pdbqt'

        from vina_docking.vina_engine import Vina123
        with tempfile.NamedTemporaryFile('w+', suffix='.pdb') as fp:
            for line in pdbqt.read_text().splitlines():
                if line.startswith('ATOM'):
                    fp.write(line)
                    fp.write('\n')
            fp.flush()

            ret = Vina123.restore_ligand(Path(fp.name), pdbqt, [docked])

        # there should be 29 atoms
        self.assertEqual(ret[0].count('ATOM'), 29)
