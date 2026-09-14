import copy
from pathlib import Path
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from reactingflow import validate_topology, MechanismError


class TopologyTests(unittest.TestCase):
    def setUp(self):
        self.data = dict(schema_version=1, kind='neutral_gas_topology',
            units={'molar_mass':'kg/mol'}, elements={'H':.001, 'O':.016},
            species=[dict(name='H2',composition={'H':2}), dict(name='O2',composition={'O':2}),
                     dict(name='H2O',composition={'H':2,'O':1})],
            reactions=[dict(id='r1',reactants={'H2':1,'O2':.5},products={'H2O':1},reversible=True)])

    def test_conservation_and_order(self):
        m = validate_topology(self.data)
        self.assertEqual(m.species, ('H2','O2','H2O'))
        self.assertAlmostEqual(m.molar_masses[2], .018)
        self.assertEqual(m.stoichiometry[0][1], -.5)
        self.assertAlmostEqual(sum(float(v)*w for v,w in zip(m.stoichiometry[0],m.molar_masses)), 0)
        self.data['reactions'] = []
        self.assertEqual(validate_topology(self.data).stoichiometry, ())

    def test_reject_invalid(self):
        changes = [lambda d: d['reactions'][0]['products'].update(H2O=2),
                   lambda d: d['reactions'][0]['reactants'].update(unknown=1),
                   lambda d: d['species'].append(d['species'][0]),
                   lambda d: d['reactions'].append(d['reactions'][0]),
                   lambda d: d['elements'].update(H=float('nan')),
                   lambda d: d['units'].update(molar_mass='g/mol'),
                   lambda d: d['reactions'][0].update(rate_constant=1),
                   lambda d: d['reactions'][0].update(reversible='true'),
                   lambda d: d['species'][0]['composition'].update(H=True),
                   lambda d: d['reactions'][0]['reactants'].update(H2=-1),
                   lambda d: d.update(schema_version=True)]
        for i, change in enumerate(changes):
            with self.subTest(i=i):
                d=copy.deepcopy(self.data)
                change(d)
                with self.assertRaises(MechanismError): validate_topology(d)

    def test_no_fixed_species_limit(self):
        self.data['species'] += [dict(name=f'inert{i}',composition={'H':1}) for i in range(100)]
        self.assertEqual(len(validate_topology(self.data).species), 103)

    def test_standalone_cli(self):
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run([sys.executable, str(root/'tools/validate_mechanism.py'),
                                 str(root/'examples/topology_demo.yaml')], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('5 species, 2 reactions', result.stdout)


if __name__ == '__main__':
    unittest.main()
