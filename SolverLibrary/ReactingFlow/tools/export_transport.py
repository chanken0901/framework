"""Offline Cantera transport tabulation; the consuming solver is entirely Fortran."""
import argparse
import hashlib
from pathlib import Path
import numpy as np
import cantera as ct
from export_mechanism import import_cantera


def export_transport(source, output, temperatures, pressure=101325., phase=None):
    """Create a NEW directory containing portable property and binary tables."""
    source, output = Path(source), Path(output)
    t = np.asarray(temperatures, dtype=float)
    if t.ndim != 1 or len(t) < 2 or not np.all(np.isfinite(t)) or np.any(t <= 0) or np.any(np.diff(t) <= 0):
        raise ValueError('At least two finite positive ascending temperatures are required')
    if not np.isfinite(pressure) or pressure <= 0:
        raise ValueError('Reference pressure must be positive and finite')
    own = import_cantera(source, phase)
    gas = ct.Solution(str(source), phase) if phase else ct.Solution(str(source))
    if gas.thermo_model != 'ideal-gas':
        raise ValueError('Only neutral ideal gas transport is supported')
    gas.transport_model = 'mixture-averaged'
    if t[0] < gas.min_temp or t[-1] > gas.max_temp:
        raise ValueError('Temperature grid outside common thermodynamic validity interval')
    if list(own.gas.names) != gas.species_names:
        raise ValueError('Species order differs from mechanism exporter')
    ns = gas.n_species
    f = lambda v: format(float(v), '.17g')
    comments = [f'# Cantera={ct.__version__} model=mixture-averaged',
                f'# source_sha256={hashlib.sha256(source.read_bytes()).hexdigest()}',
                f'# canonical_sha256={own.canonical_sha256}',
                '# SI: temperature K, pressure Pa, molar mass kg/mol, viscosity Pa.s, conductivity W/(m.K), D m2/s']
    species = [name+' '+f(mass) for name,mass in zip(gas.species_names,own.gas.molar_masses)]
    props = ['RF_TRANSPORT_TABLE_V1',*comments,f'{ns} {len(t)}',*species]
    binary = ['RF_BINARY_TABLE_V1',*comments,str(ns),f'{len(t)} {f(pressure)}',*species]
    for temp in t:
        gas.TPX = temp,pressure,np.ones(ns)
        mu = gas.species_viscosities.copy()
        dij = gas.binary_diff_coeffs.copy()
        k = np.empty(ns)
        for i in range(ns):
            pure = np.zeros(ns); pure[i] = 1
            gas.TPX = temp,pressure,pure
            k[i] = gas.thermal_conductivity
        if not all(np.all(np.isfinite(v)) and np.all(v>0) for v in (mu,k,dij)):
            raise ValueError('Nonpositive/nonfinite Cantera transport value')
        if not np.allclose(dij,dij.T,rtol=1.e-12,atol=0):
            raise ValueError('Binary diffusion matrix is not symmetric')
        props.append(f(temp));binary.append(f(temp))
        props.extend(f'{name} {f(a)} {f(b)}' for name,a,b in zip(gas.species_names,mu,k))
        binary.extend(f'{gas.species_names[i]} {gas.species_names[j]} {f(dij[i,j])}'
                      for i in range(ns) for j in range(i+1,ns))
    # Validate and evaluate everything before creating output. Never overwrite.
    output.mkdir(parents=False,exist_ok=False)
    (output/'properties.rf').write_text('\n'.join(props)+'\n',encoding='ascii')
    (output/'binary.rf').write_text('\n'.join(binary)+'\n',encoding='ascii')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mechanism',type=Path)
    parser.add_argument('output_directory',type=Path)
    parser.add_argument('--phase')
    parser.add_argument('--tmin',type=float,default=300.)
    parser.add_argument('--tmax',type=float,default=3500.)
    parser.add_argument('--points',type=int,default=201)
    parser.add_argument('--pressure',type=float,default=101325.)
    args=parser.parse_args()
    if args.points<2 or not np.isfinite(args.tmin) or not np.isfinite(args.tmax) or not 0<args.tmin<args.tmax:
        parser.error('Require 0<tmin<tmax and points>=2')
    export_transport(args.mechanism,args.output_directory,np.geomspace(args.tmin,args.tmax,args.points),args.pressure,args.phase)
    print('[OK] Created properties.rf and binary.rf (offline data; no runtime Cantera dependency)')
