from pathlib import Path
import sys
import tempfile
import unittest
import math

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from reactingflow.thermo import NASA, IdealGas, R
from reactingflow.mechanism import MechanismError
from reactingflow.importer import import_cantera
from reactingflow.units import ReferenceScales


class ThermoTests(unittest.TestCase):
    def setUp(self):
        self.nasa = NASA('NASA7',(200,6000),((3.5,0,0,0,0,-1000,2),))
        self.gas = IdealGas(('A',),(.02,),(self.nasa,))

    def test_constant_cp_and_inverse(self):
        for t in (200,300,1000,5999,6000):
            p=self.gas.properties(t,[1])
            self.assertAlmostEqual(p['cp'],3.5*R/.02)
            self.assertAlmostEqual(p['e'],(2.5*t-1000)*R/.02)
            self.assertAlmostEqual(self.gas.temperature_from_energy(p['e'],[1]),t,places=5)
            self.assertAlmostEqual(p['h']-t*p['s'],p['g'])

    def test_nasa9_equivalent(self):
        nine=NASA('NASA9',(200,6000),((0,0,3.5,0,0,0,0,-1000,2),))
        for a,b in zip(nine.molar(1200),self.nasa.molar(1200)):
            self.assertAlmostEqual(a,b)

    def test_rejections(self):
        for y in ([.9],[-1],[float('nan')],[True], [1,0]):
            with self.assertRaises(MechanismError): self.gas.properties(300,y)
        for t in (0,199,6001,float('nan')):
            with self.assertRaises(MechanismError): self.nasa.molar(t)
        with self.assertRaises(MechanismError): self.gas.temperature_from_energy(1e20,[1])
        with self.assertRaises(MechanismError): NASA('NASA7',(200,100),((1,)*7,))

    def test_mixing_entropy(self):
        gas=IdealGas(('A','B'),(.02,.02),(self.nasa,self.nasa))
        p=gas.properties(500,[.5,.5])
        pure=self.gas.properties(500,[1])
        self.assertAlmostEqual(p['s']-pure['s'],R/.02*math.log(2))

    def test_reference_scales(self):
        scales=ReferenceScales(2,10,.5,300)
        self.assertEqual(scales.to_si(3,'pressure'),600)
        self.assertEqual(scales.to_si(2,'time'),.1)
        self.assertEqual(scales.from_si(-200,'energy'),-2)
        with self.assertRaises(MechanismError): ReferenceScales(0,10,.5,300)
        with self.assertRaises(MechanismError): scales.scale('unknown')


class CanteraReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import cantera as ct
        except ImportError:
            raise unittest.SkipTest('Install Cantera 3.2.0 for reference comparisons')
        cls.ct=ct

    def compare(self, reference, own):
        for t in (350,800,1500,2800):
            y=[1/len(own.names)]*len(own.names)
            reference.TPY=t,230000,y
            p=own.properties(t,y,230000)
            for key,value in [('cp',reference.cp_mass),('cv',reference.cv_mass),
                ('h',reference.enthalpy_mass),('e',reference.int_energy_mass),
                ('s',reference.entropy_mass),('g',reference.gibbs_mass),
                ('density',reference.density),('sound_speed',reference.sound_speed)]:
                self.assertAlmostEqual(p[key]/value,1,places=10,msg=key)
            self.assertAlmostEqual(own.temperature_from_energy(p['e'],y),t,places=5)

    def test_hydrogen_nasa7_and_rate_preservation(self):
        ct=self.ct
        path=Path(ct.__file__).parent/'data/h2o2.yaml'
        imported=import_cantera(path)
        self.compare(ct.Solution(str(path)),imported.gas)
        restored=ct.Solution(yaml=imported.canonical_yaml)
        self.assertEqual(restored.n_reactions,len(imported.topology.reaction_ids))
        self.assertEqual(restored.n_reactions,len(imported.reaction_types))
        self.assertEqual(len(imported.source_sha256),64)
        self.assertFalse(any(line.startswith('date:') for line in imported.canonical_yaml.splitlines()))

    def test_neutral_nasa9(self):
        ct=self.ct
        species=[s for s in ct.Species.list_from_file('airNASA9.yaml') if s.name in ('N2','O2','NO')]
        reference=ct.Solution(thermo='ideal-gas',species=species)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'neutral.yaml'
            reference.write_yaml(path)
            self.compare(reference,import_cantera(path).gas)


if __name__ == '__main__':
    unittest.main()
