#!/usr/bin/env python3
"""Inspect Cantera-format mechanism with independent ReactingFlow thermodynamics."""
import argparse
from pathlib import Path
import sys
import json
from dataclasses import asdict
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from reactingflow.importer import import_cantera


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--phase')
    parser.add_argument('--rates',action='store_true',help='Also evaluate instantaneous reaction/source rates (no time integration)')
    parser.add_argument('--temperature',type=float,default=1000)
    parser.add_argument('--pressure',type=float,default=101325)
    parser.add_argument('--mass-fractions',required=True,help='JSON mapping, e.g. {"H2":0.1,"O2":0.9}')
    args=parser.parse_args()
    try:
        imported=import_cantera(args.input,args.phase)
        composition=json.loads(args.mass_fractions)
        if not isinstance(composition,dict) or set(composition)-set(imported.gas.names):
            raise ValueError('Mass fractions contain unknown species')
        y=[composition.get(s,0) for s in imported.gas.names]
        props=imported.gas.properties(args.temperature,y,args.pressure)
        report=dict(species=imported.gas.names, reaction_types=dict(Counter(imported.reaction_types)),
                    source_sha256=imported.source_sha256,canonical_sha256=imported.canonical_sha256,
                    temperature_bounds=imported.gas.temperature_bounds,properties_si=props,
                    recovered_temperature=imported.gas.temperature_from_energy(props['e'],y),
                    kinetics_evaluated=False)
        if args.rates:
            rates=imported.kinetics.evaluate(args.temperature,props['density'],y)
            report['kinetics_evaluated']=True
            report['rates_si']=asdict(rates)
        print(json.dumps(report,indent=2,allow_nan=False))
    except ImportError:
        parser.exit(1,'Install input adapter: python -m pip install cantera==3.2.0\n')
    except Exception as exc:
        parser.exit(1,f'[ERROR] {exc}\n')


if __name__=='__main__':
    main()
