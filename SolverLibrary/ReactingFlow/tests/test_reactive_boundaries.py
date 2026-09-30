"""Reacting inlet, analytic reaction/advection, and ZND-initialized CFD checks."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import cantera as ct
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from export_mechanism import export,import_cantera


class ReactiveBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.environ.get('RF_FORTRAN_BUILD'): raise unittest.SkipTest('Set RF_FORTRAN_BUILD')
        cls.exe=Path(os.environ['RF_FORTRAN_BUILD'])/('rf_flow1d.exe' if os.name=='nt' else 'rf_flow1d')

    def run_case(self,tmp,data,controls,yl,yr,ok=True):
        inp=tmp/'case.in';out=tmp/'flow.csv'
        if out.exists(): out.unlink()
        inp.write_text('&flow1d '+controls+' /\n'+' '.join(map(str,yl))+'\n'+' '.join(map(str,yr))+'\n')
        run=subprocess.run([str(self.exe),str(data),str(inp),str(out)],capture_output=True,text=True,timeout=600)
        if not ok:
            self.assertNotEqual(run.returncode,0)
            if out.exists(): self.assertNotIn('# SUCCESS',out.read_text())
            return
        self.assertEqual(run.returncode,0,run.stdout+run.stderr)
        text=out.read_text();self.assertTrue(text.rstrip().endswith('# SUCCESS'))
        rows=[s for s in text.splitlines() if not s.startswith('#')]
        values=np.loadtxt(rows[1:],delimiter=',',ndmin=2)
        diag=dict(s[2:].split('=',1) for s in text.splitlines() if s.startswith('# ') and '=' in s)
        for key in ('mass_error','momentum_error','energy_error','element_error'):
            self.assertLess(float(diag[key]),1.e-8)
        return values

    def test_inlet_both_sides_and_rejection(self):
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);data=tmp/'m.rf'
            source=Path(ct.__file__).parent/'data/h2o2.yaml';export(import_cantera(source),data)
            gas=ct.Solution(str(source));gas.TPX=1100,101325,'N2:1';y=gas.Y
            for side,sgn in [('left',1),('right',-1)]:
                for speed in (100.,2000.):
                    u=sgn*speed
                    controls=f"nx=8,end_time=1.e-7,left_velocity={u},right_velocity={u}," \
                        f"{side}_bc='reacting_inlet',{side}_boundary_temperature=1100," \
                        f"{side}_boundary_pressure=101325,{side}_boundary_velocity={u}," \
                        f"{side}_boundary_y="+','.join(map(str,y))
                    for extra in ('',",transport_model='constant',viscosity=1.e-5,thermal_conductivity=.03"):
                        values=self.run_case(tmp,data,controls+extra,y,y)
                        np.testing.assert_allclose(values[:,4],u,rtol=1.e-10)
                        np.testing.assert_allclose(values[:,5],1100,rtol=1.e-10)
                    self.run_case(tmp,data,controls.replace(f'{side}_boundary_velocity={u}',f'{side}_boundary_velocity={-u}'),y,y,False)
                    self.run_case(tmp,data,controls.replace(f'{side}_boundary_temperature=1100,',''),y,y,False)

    def test_reaction_advection_convergence(self):
        species=[]
        for name in ('A','B'):
            s=ct.Species(name,{'H':2})
            s.thermo=ct.NasaPoly2(200,4000,101325,[1000,3.5,0,0,0,0,0,0,3.5,0,0,0,0,0,0])
            species.append(s)
        reaction=ct.Reaction(equation='A => B',rate=ct.ArrheniusRate(500.,0.,0.))
        gas=ct.Solution(thermo='ideal-gas',kinetics='gas',species=species,reactions=[reaction])
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);source=tmp/'analytic.yaml';gas.write_yaml(source)
            data=tmp/'m.rf';export(import_cantera(source),data)
            errors=[]
            for nx in (40,80):
                controls=f"nx={nx},end_time=.001,max_dt=.00002,write_every=100000,chemistry=.true.," \
                    "reconstruction='muscl',left_velocity=100,right_velocity=100,left_bc='reacting_inlet'," \
                    "left_boundary_temperature=1100,left_boundary_pressure=101325,left_boundary_velocity=100," \
                    "left_boundary_y=1,0"
                rows=self.run_case(tmp,data,controls,[1.,0.],[.2,.8]);last=rows[rows[:,0]==rows[-1,0]]
                # Cell-average analytic solution, including continuous fresh inlet supply.
                x=last[:,2,None]+np.linspace(-.5+.5/64,.5-.5/64,64)[None,:]/nx
                exact=np.where(x<.1,np.exp(-500*x/100),np.where(x<.6,1.,.2)*np.exp(-.5)).mean(1)
                errors.append(np.mean(abs(last[:,7]-exact)))
                np.testing.assert_allclose(last[:,5],1100,rtol=1.e-7)
                np.testing.assert_allclose(last[:,6],101325,rtol=1.e-7)
            self.assertLess(errors[1],.8*errors[0])
            self.assertLess(errors[1],.025)

    def test_initial_profile_validation(self):
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);source=Path(ct.__file__).parent/'data/h2o2.yaml'
            mechanism=import_cantera(source);data=tmp/'m.rf';export(mechanism,data)
            gas=ct.Solution(str(source));gas.TPX=1100,101325,'N2:1';y=gas.Y
            header='RF_FLOW_PROFILE_V1\n'+mechanism.canonical_sha256+'\n10 8 1\n'
            row='1100 101325 100 '+' '.join(map(str,y))+'\n'
            original=header+row*8;path=tmp/'profile.rf';path.write_text(original)
            controls="nx=8,end_time=1.e-8,initial_profile='profile.rf'"
            rows=self.run_case(tmp,data,controls,y,y)
            np.testing.assert_allclose(rows[:,4],100,rtol=1.e-10)
            for damaged in (original.replace(mechanism.canonical_sha256,'wrong'),
                            original.replace('10 8 1','10 9 1'),original.replace('10 8 1','10 8 2'),
                            header+row*7,original+'EXTRA\n',original.replace('1100 101325','1100 -1',1)):
                path.write_text(damaged);self.run_case(tmp,data,controls,y,y,False)

    def test_znd_stationary_wave_grid_refinement(self):
        self.check_znd_wave('muscl')

    def test_znd_first_order_grid_time_chemistry(self):
        self.check_znd_wave('first_order', sensitivity=True)

    def check_znd_wave(self, reconstruction, sensitivity=False):
        # Independent Cantera/Radau steady ZND reference in shock-fixed downstream coordinate.
        source=Path(ct.__file__).parent/'data/h2o2.yaml'
        gas=ct.Solution(str(source));gas.TPX=300,101325,'H2:2,O2:1,N2:3.76'
        y0=gas.Y.copy();rho0=gas.density;h0=gas.enthalpy_mass;speed=5*gas.sound_speed;flux=rho0*speed
        def jump(ratio):
            gas.DPY=rho0*ratio,101325+rho0*speed**2*(1-1/ratio),y0
            return gas.enthalpy_mass+.5*(speed/ratio)**2-h0-.5*speed**2
        jump(brentq(jump,1.00001,rho0*speed**2/101325))
        initial=np.r_[gas.P,gas.density,0.,y0]
        def rhs(time,state):
            gas.DPY=state[1],state[0],state[3:]
            u=flux/state[1];dy=gas.net_production_rates*gas.molecular_weights/state[1]
            sigma=np.dot(gas.mean_molecular_weight/gas.molecular_weights-
                         gas.partial_molar_enthalpies/gas.molecular_weights/(gas.cp_mass*gas.T),dy)
            drho=-state[1]*sigma/(1-u*u/gas.sound_speed**2)
            return np.r_[u*u*drho,drho,u,dy]
        def jac(time,state):
            base=rhs(time,state);matrix=np.empty((len(state),len(state)))
            for j in range(len(state)):
                shifted=state.copy();step=1.e-7*max(abs(state[j]),1.e-5)
                shifted[j]+=step;matrix[:,j]=(rhs(time,shifted)-base)/step
            return matrix
        ref=solve_ivp(rhs,[0,1.e-5],initial,method='Radau',rtol=1.e-10,atol=1.e-15,max_step=1.e-8,
                      dense_output=True,jac=jac)
        self.assertTrue(ref.success,ref.message)
        domain=.002;shock=.0005
        def profile(x):
            if x<shock: return np.r_[300.,101325.,speed,y0]
            tau=brentq(lambda t: ref.sol(t)[2]-(x-shock),0,ref.t[-1],xtol=1.e-16)
            state=ref.sol(tau);gas.DPY=state[1],state[0],state[3:]
            return np.r_[gas.T,gas.P,flux/state[1],gas.Y]
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);data=tmp/'m.rf';mechanism=import_cantera(source);export(mechanism,data)
            errors=[];solutions=[]
            cases=[(128,2.e-9,1.e-9,1.e-16,1.e-8),(256,2.e-9,1.e-9,1.e-16,1.e-8)]
            if sensitivity:
                cases.extend([(128,1.e-9,1.e-9,1.e-16,1.e-8),(128,2.e-9,1.e-10,1.e-17,1.e-9)])
            for nx,max_dt,rtol,atoly,atolt in cases:
                states=np.array([profile((i+.5)*domain/nx) for i in range(nx)])
                (tmp/'initial.rf').write_text('RF_FLOW_PROFILE_V1\n'+mechanism.canonical_sha256+
                    f'\n{len(y0)} {nx} {domain}\n'+'\n'.join(' '.join(map(str,s)) for s in states)+'\n')
                right=profile(domain)
                controls=f"nx={nx},length={domain},initial_profile='initial.rf',end_time=2.e-7,max_dt={max_dt}," \
                    f"chemistry_rtol={rtol},chemistry_atol_species={atoly},chemistry_atol_temperature={atolt}," \
                    f"write_every=100000,chemistry=.true.,reconstruction='{reconstruction}'," \
                    "left_bc='reacting_inlet',left_boundary_temperature=300,left_boundary_pressure=101325," \
                    f"left_boundary_velocity={speed},left_boundary_y="+','.join(map(str,y0))+','+ \
                    "right_bc='characteristic',"+f"right_boundary_temperature={right[0]},right_boundary_pressure={right[1]}," \
                    f"right_boundary_velocity={right[2]},right_boundary_y="+','.join(map(str,right[3:]))
                rows=self.run_case(tmp,data,controls,y0,y0);last=rows[rows[:,0]==rows[-1,0]]
                errors.append(np.mean(abs(last[:,6]-states[:,1]))/max(states[:,1]))
                solutions.append(last[:,6].copy())
                print(f'ZND {reconstruction}: nx={nx}, max_dt={max_dt}, rtol={rtol}, pressure error={errors[-1]}',flush=True)
                self.assertTrue(np.all(last[:,6]>0))
                front=(last[:-1,2]+last[1:,2])[np.argmax(abs(np.diff(last[:,6])))]/2
                self.assertLessEqual(abs(front-shock),2*domain/nx)
            self.assertLess(errors[1],errors[0],str(errors))
            self.assertLess(errors[1],.03,str(errors))
            if sensitivity:
                for label,index in [('time',2),('chemistry',3)]:
                    difference=np.mean(abs(solutions[index]-solutions[0]))/max(solutions[0])
                    print(f'ZND {reconstruction} {label} sensitivity={difference}',flush=True)
                    self.assertLess(difference,1.e-3)
