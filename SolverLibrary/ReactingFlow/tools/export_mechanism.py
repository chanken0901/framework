"""Input conversion only: Cantera YAML -> portable SI data for Fortran."""
import argparse
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'reference/python'))
from reactingflow.importer import import_cantera

KINDS = ['Arrhenius', 'three-body', 'Lindemann', 'Troe', 'SRI', 'PLOG', 'Chebyshev']


def export(mechanism, path):
    lines = ['RFMECH1', mechanism.source_sha256, mechanism.canonical_sha256]
    def row(values):
        lines.append(' '.join(format(v, '.17g') if isinstance(v, float) else str(v) for v in values))
    def arr(a):
        row([a.log_a if math.isfinite(a.log_a) else -1.e300, a.b, a.ea_over_r, a.sign] if a else [-1.e300, 0., 0., 1])
    gas = mechanism.gas
    row([len(gas.names), len(mechanism.topology.elements), len(mechanism.kinetics.reactions)])
    for name, mass, nasa, atoms in zip(gas.names, gas.molar_masses, gas.thermo, mechanism.topology.compositions):
        if len(name) > 128 or any(c.isspace() or c in "'\",/" for c in name):
            raise ValueError(f'Unsupported portable species name: {name}')
        row([name, mass, int(nasa.model[4:]), len(nasa.coefficients), nasa.reference_pressure])
        row(nasa.bounds)
        for a in nasa.coefficients: row(list(a)+[0.]*(9-len(a)))
        row(atoms)
    for reaction in mechanism.kinetics.reactions:
        r = reaction.rate
        row([KINDS.index(r.kind)+1, int(reaction.reversible), len(r.parameters), len(r.plog),
             len(r.coefficients), len(r.coefficients[0]) if r.coefficients else 0])
        arr(r.high); arr(r.low)
        row(reaction.reactants); row(reaction.products); row(reaction.orders)
        row(r.efficiencies or [0.]*len(gas.names))
        if r.parameters: row(r.parameters)
        for p, terms in r.plog:
            row([p, len(terms)])
            for a in terms: arr(a)
        if r.coefficients:
            row(r.bounds)
            for a in r.coefficients: row(a)
    # Refuse overwrite. No chemistry is evaluated by this converter.
    with Path(path).open('x', encoding='ascii', newline='\n') as out:
        out.write('\n'.join(lines)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mechanism', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--phase')
    args = parser.parse_args()
    export(import_cantera(args.mechanism, args.phase), args.output)
