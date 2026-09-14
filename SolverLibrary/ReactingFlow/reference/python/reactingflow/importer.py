"""Offline Cantera INPUT adapter. Not linked to the Fortran solver.

Supported reactions are compiled to independent SI rate evaluators at import.
"""
from dataclasses import dataclass
import hashlib
from pathlib import Path
from .mechanism import MechanismError, validate_topology
from .thermo import NASA, IdealGas
from .kinetics import Kinetics, compile_kinetics


@dataclass(frozen=True)
class ImportedMechanism:
    gas: IdealGas
    topology: object
    reaction_types: tuple[str, ...]
    canonical_yaml: str
    source_sha256: str
    canonical_sha256: str
    kinetics: Kinetics


def import_cantera(path, phase=None):
    import cantera as ct
    path = Path(path).resolve()
    if not path.is_file():
        raise MechanismError(f'Mechanism file not found: {path}')
    source = path.read_bytes()
    solution = ct.Solution(str(path), phase) if phase else ct.Solution(str(path))
    if solution.thermo_model != 'ideal-gas':
        raise MechanismError('Only ideal-gas phases are supported')
    species, thermo = [], []
    for s in solution.species():
        if s.charge or any(int(n)!=n for n in s.composition.values()):
            raise MechanismError('Only neutral species with integer atom counts are supported')
        species.append(dict(name=s.name, composition={e:int(n) for e,n in s.composition.items()}))
        data = s.thermo.input_data
        if data['model'] not in ('NASA7','NASA9'):
            raise MechanismError(f'{s.name}: unsupported thermo {data["model"]}')
        thermo.append(NASA(data['model'], tuple(data['temperature-ranges']),
                          tuple(tuple(row) for row in data['data']), s.thermo.reference_pressure))
    reactions = [dict(id=f'r{i+1}', reactants=dict(r.reactants), products=dict(r.products),
                      reversible=r.reversible) for i,r in enumerate(solution.reactions())]
    topology = validate_topology(dict(schema_version=1,kind='neutral_gas_topology',
        units={'molar_mass':'kg/mol'},
        elements={e:float(m)/1000 for e,m in zip(solution.element_names,solution.atomic_weights)},
        species=species, reactions=reactions))
    # Expanded canonical input preserves rate/transport data and external references
    # semantically. Original source bytes and canonical contents have separate hashes.
    # Exclude generated wall-clock metadata from the reproducibility hash.
    canonical = '\n'.join(line for line in solution.write_yaml(header=False).splitlines()
                          if not line.startswith('date:')) + '\n'
    gas = IdealGas(topology.species, topology.molar_masses, tuple(thermo))
    kinetics = compile_kinetics(gas, solution)
    return ImportedMechanism(gas,
        topology, tuple(r.reaction_type for r in solution.reactions()), canonical,
        hashlib.sha256(source).hexdigest(), hashlib.sha256(canonical.encode()).hexdigest(), kinetics)
