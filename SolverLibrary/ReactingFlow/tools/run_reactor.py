"""Development-only 0D reactor CLI. Does not modify a CFD case or manifest."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from reactingflow.importer import import_cantera
from reactingflow.mechanism import MechanismError
from reactingflow.reactor import integrate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mechanism', type=Path)
    parser.add_argument('--phase')
    parser.add_argument('--mode', choices=['constant_volume', 'constant_pressure'], default='constant_volume')
    parser.add_argument('--temperature', type=float, required=True, help='K')
    parser.add_argument('--pressure', type=float, default=101325., help='Pa')
    parser.add_argument('--mass-fractions', required=True, help='JSON map; must sum to one')
    parser.add_argument('--end-time', type=float, required=True, help='s')
    parser.add_argument('--rtol', type=float, default=1.e-7)
    parser.add_argument('--atol-species', type=float, default=1.e-14)
    parser.add_argument('--atol-temperature', type=float, default=1.e-6)
    parser.add_argument('--max-step', type=float)
    parser.add_argument('--max-steps', type=int, default=100000)
    parser.add_argument('--ignition-rise', type=float, default=400.)
    parser.add_argument('--output', type=Path, help='New JSON file (existing files are never overwritten)')
    args = parser.parse_args()
    try:
        if args.output and args.output.exists():
            raise MechanismError(f'Output already exists: {args.output}')
        mechanism = import_cantera(args.mechanism, args.phase)
        values = json.loads(args.mass_fractions)
        if not isinstance(values, dict) or set(values)-set(mechanism.gas.names):
            raise MechanismError('Mass fractions must map known species names to numbers')
        result = integrate(mechanism, temperature=args.temperature, pressure=args.pressure,
            mass_fractions=[values.get(s, 0.) for s in mechanism.gas.names], end_time=args.end_time,
            mode=args.mode, rtol=args.rtol, atol_species=args.atol_species,
            atol_temperature=args.atol_temperature, max_step=args.max_step,
            max_steps=args.max_steps, ignition_temperature_rise=args.ignition_rise)
        report = dict(source_sha256=mechanism.source_sha256,
                      canonical_sha256=mechanism.canonical_sha256,
                      settings=vars(args).copy(), result=asdict(result))
        report['settings']['mechanism'] = str(args.mechanism.resolve())
        report['settings']['output'] = str(args.output.resolve()) if args.output else None
        text = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n'
        if args.output:
            with args.output.open('x', encoding='utf-8') as stream:
                stream.write(text)
            print(f'[OK] {args.output.resolve()} ({result.diagnostics["steps"]} accepted steps)')
        else:
            print(text, end='')
        return 0
    except (MechanismError, ValueError, OSError, ImportError) as exc:
        print(f'[ERROR] {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
