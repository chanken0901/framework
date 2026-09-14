from pathlib import Path
import sys
import unittest
import math
import subprocess
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from reactingflow.importer import import_cantera
from reactingflow.kinetics import compile_kinetics, Kinetics, Reaction, Rate, Arrhenius
from reactingflow.thermo import IdealGas, NASA
from reactingflow.mechanism import MechanismError


class AnalyticRatesTests(unittest.TestCase):
    def test_bimolecular_mass_action_without_cantera(self):
        nasa=NASA('NASA7',(200,3000),((3.5,0,0,0,0,0,0),))
        gas=IdealGas(('A','B'),(.01,.02),(nasa,nasa))
        reaction=Reaction(Rate('Arrhenius',high=Arrhenius(math.log(4),0,0)),(2,0),(0,1),(2,0),False)
        result=Kinetics(gas,(reaction,)).evaluate(1000,1,[.5,.5])
        self.assertAlmostEqual(result.forward[0],10000)
        self.assertEqual(result.reverse[0],0)
        self.assertAlmostEqual(result.mass_production[0],-200)
        self.assertAlmostEqual(sum(result.mass_production),0)


class ReferenceRatesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import cantera as ct
            import numpy as np
        except ImportError:
            raise unittest.SkipTest('Cantera 3.2.0 required for kinetics reference tests')
        cls.ct,cls.np=ct,np
        cls.path=Path(ct.__file__).parent/'data/h2o2.yaml'
        cls.imported=import_cantera(cls.path)

    def check_state(self,ref,own,t,p,y=None):
        if y is None: y=[1/ref.n_species]*ref.n_species
        ref.TPY=t,p,y
        result=own.evaluate(t,float(ref.density),list(ref.Y))
        for actual,expected in [(result.forward,ref.forward_rates_of_progress*1000),
            (result.reverse,ref.reverse_rates_of_progress*1000),
            (result.net,ref.net_rates_of_progress*1000),
            (result.molar_production,ref.net_production_rates*1000),
            (result.mass_production,ref.net_production_rates*ref.molecular_weights)]:
            self.np.testing.assert_allclose(actual,expected,rtol=3e-10,atol=1e-12)
        scale=max(1,sum(abs(v) for v in result.mass_production))
        self.assertLess(abs(sum(result.mass_production)),1e-13*scale)
        for element in ref.element_names:
            terms=[float(ref.n_atoms(s,element))*v for s,v in zip(ref.species_names,result.molar_production)]
            self.assertLess(abs(math.fsum(terms)),1e-12*max(1,sum(abs(v) for v in terms)))
        self.assertAlmostEqual(result.heat_release/float(ref.heat_release_rate),1,places=9)
        return result

    def test_hydrogen_complete(self):
        ref=self.ct.Solution(str(self.path))
        for t in (400,900,1500,2800):
            for p in (1e3,1e5,1e7):
                with self.subTest(t=t,p=p): self.check_state(ref,self.imported.kinetics,t,p)

    def test_gri30_complete(self):
        path=Path(self.ct.__file__).parent/'data/gri30.yaml'
        own=import_cantera(path)
        ref=self.ct.Solution(str(path))
        self.assertEqual(len(own.reaction_types),325)
        for t in (500,1200,2500):
            for p in (1e4,1e5,1e7): self.check_state(ref,own.kinetics,t,p)

    def test_evaluation_has_no_cantera_calls(self):
        with patch.object(self.ct,'Solution',side_effect=AssertionError('runtime Cantera call')):
            result=self.imported.kinetics.evaluate(1000,1,[.1]*10)
        self.assertEqual(len(result.net),29)

    def test_cli_rates(self):
        root=Path(__file__).resolve().parents[1]
        result=subprocess.run([sys.executable,str(root/'tools/inspect_thermo.py'),str(self.path),
            '--temperature','1000','--mass-fractions','{"H2":0.1,"O2":0.9}','--rates'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('"kinetics_evaluated": true',result.stdout)

    def test_pure_species(self):
        ref=self.ct.Solution(str(self.path))
        y=[0.]*ref.n_species
        y[ref.species_index('H2')]=1
        result=self.imported.kinetics.evaluate(1000,.1,y)
        self.np.testing.assert_allclose(result.molar_production,
            self.check_state(ref,self.imported.kinetics,1000,.1*self.imported.gas.properties(1000,y)['gas_constant']*1000,y).molar_production,
            rtol=1e-10,atol=1e-15)

    def test_zero_effective_collider(self):
        ct=self.ct
        rate=ct.LindemannRate(low=ct.ArrheniusRate(1e12,0,0),high=ct.ArrheniusRate(1e8,0,0))
        reaction=ct.Reaction(equation='H + O2 (+M) <=> HO2 (+M)',rate=rate)
        reaction.third_body.default_efficiency=0
        reaction.third_body.efficiencies={'AR':1}
        ref=ct.Solution(thermo='ideal-gas',kinetics='gas',species=ct.Solution(str(self.path)).species(),reactions=[reaction])
        ref.TPX=1000,1e5,'H:1,O2:1,HO2:1'
        own=compile_kinetics(self.imported.gas,ref)
        result=own.evaluate(ref.T,ref.density,list(ref.Y))
        self.assertEqual(result.forward,(0.,))
        self.assertEqual(result.reverse,(0.,))
        self.assertEqual(result.heat_release,0.)

    def synthetic(self,rate,falloff=False,orders=None):
        ct=self.ct
        equation='H + O2 (+M) <=> HO2 (+M)' if falloff else 'H + O2 => HO2'
        reaction=ct.Reaction(equation=equation,rate=rate)
        if falloff:
            reaction.third_body.efficiencies={'H2O':6,'AR':.7}
        if orders:
            reaction.allow_nonreactant_orders=True
            reaction.orders=orders
        ref=ct.Solution(thermo='ideal-gas',kinetics='gas',species=ct.Solution(str(self.path)).species(),reactions=[reaction])
        return ref,compile_kinetics(self.imported.gas,ref)

    def test_falloff(self):
        ct=self.ct
        for cls,params in [(ct.LindemannRate,()),(ct.TroeRate,(.5,100,1000)),
                           (ct.TroeRate,(.5,100,1000,5000)),(ct.SriRate,(1.2,100,1500)),
                           (ct.SriRate,(1.2,100,1500,.9,.1))]:
            rate=cls(low=ct.ArrheniusRate(1e12,.1,2e6),high=ct.ArrheniusRate(1e8,-.2,1e6),falloff_coeffs=params)
            ref,own=self.synthetic(rate,True)
            for t in (500,1000,2500):
                for p in (1e-2,1e5,1e12):
                    with self.subTest(cls=cls.__name__,t=t,p=p): self.check_state(ref,own,t,p)

    def test_plog_duplicate_pressure_and_clamping(self):
        ct=self.ct
        rates=[(1e4,ct.ArrheniusRate(1e8,0,1e6)),(1e4,ct.ArrheniusRate(-1e7,0,1e6)),
               (1e6,ct.ArrheniusRate(2e9,.1,2e6))]
        ref,own=self.synthetic(ct.PlogRate(rates))
        for t in (500,1000,2500):
            for p in (1e2,1e4,1e5,1e6,1e8): self.check_state(ref,own,t,p)

    def test_chebyshev(self):
        ct=self.ct
        ref,own=self.synthetic(ct.ChebyshevRate(temperature_range=(300,3000),pressure_range=(1e3,1e7),
                                               data=[[8,.2,-.02],[.5,.1,.04],[.03,.01,.005]]))
        for t in (300,700,1800,3000):
            for p in (1e3,1e5,1e7): self.check_state(ref,own,t,p)
        with self.assertRaises(MechanismError): self.check_state(ref,own,1000,1e8)

    def test_irreversible_custom_orders(self):
        ref,own=self.synthetic(self.ct.ArrheniusRate(1e8,.1,2e6),orders={'H':.5,'O2':1.2,'H2O':.1})
        self.check_state(ref,own,1000,1e5)

    def test_zero_species_and_invalid_state(self):
        ref=self.ct.Solution(str(self.path))
        y=[0.]*ref.n_species
        y[ref.species_index('H2')]=.1; y[ref.species_index('O2')]=.9
        self.check_state(ref,self.imported.kinetics,1000,1e5,y)
        with self.assertRaises(MechanismError): self.imported.kinetics.evaluate(1000,-1,y)
        with self.assertRaises(MechanismError): self.imported.kinetics.evaluate(1000,1,[math.nan]*len(y))

    def test_equilibrium_detailed_balance(self):
        ref=self.ct.Solution(str(self.path))
        ref.TPX=1500,1e5,'H2:2,O2:1,N2:4'
        ref.equilibrate('TP')
        result=self.imported.kinetics.evaluate(ref.T,ref.density,list(ref.Y))
        scale=max(result.forward)
        self.assertLess(max(abs(v) for v in result.net)/scale,1e-10)

    def test_unsupported_rejected(self):
        ct=self.ct
        rate=ct.TroeRate(low=ct.ArrheniusRate(1e12,0,0),high=ct.ArrheniusRate(1e8,0,0),falloff_coeffs=[.5,100,1000])
        rate.chemically_activated=True
        with self.assertRaises(MechanismError): self.synthetic(rate,True)


if __name__=='__main__':
    unittest.main()
