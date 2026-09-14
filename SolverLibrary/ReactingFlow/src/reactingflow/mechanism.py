"""Validated neutral-gas stoichiometry, not a reaction-rate/thermo evaluator.

This internal topology contract is deliberately not advertised as Cantera YAML.
Unknown fields are rejected so unsupported kinetics cannot silently disappear.
"""
from dataclasses import dataclass
from fractions import Fraction
import math


class MechanismError(ValueError):
    pass


def mapping(value, keys, label):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise MechanismError(f'{label}: required keys are {sorted(keys)} (no extras)')
    return value


def name(value, label):
    if not isinstance(value, str) or not value or value.strip() != value:
        raise MechanismError(f'{label}: expected a nonempty trimmed string')
    return value


def positive(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MechanismError(f'{label}: expected a finite positive number')
    if not math.isfinite(value) or value <= 0:
        raise MechanismError(f'{label}: expected a finite positive number')
    return value


@dataclass(frozen=True)
class Mechanism:
    elements: tuple[str, ...]
    species: tuple[str, ...]
    molar_masses: tuple[float, ...]  # kg/mol; derived, never independently supplied
    compositions: tuple[tuple[int, ...], ...]  # species x elements
    reaction_ids: tuple[str, ...]
    reactants: tuple[tuple[Fraction, ...], ...]  # reactions x species
    products: tuple[tuple[Fraction, ...], ...]
    reversible: tuple[bool, ...]

    @property
    def stoichiometry(self):
        """Net products minus reactants, exact rational coefficients."""
        return tuple(tuple(p-r for p, r in zip(products, reactants))
                     for products, reactants in zip(self.products, self.reactants))


def validate_topology(document):
    mapping(document, ('schema_version', 'kind', 'units', 'elements', 'species', 'reactions'), 'topology')
    if type(document['schema_version']) is not int or document['schema_version'] != 1:
        raise MechanismError('schema_version must be integer 1')
    if document['kind'] != 'neutral_gas_topology':
        raise MechanismError('only neutral_gas_topology is implemented')
    if document['units'] != {'molar_mass': 'kg/mol'}:
        raise MechanismError('molar mass unit must explicitly be kg/mol')
    elements = document['elements']
    if not isinstance(elements, dict) or not elements:
        raise MechanismError('elements must be a nonempty atomic-mass mapping')
    for symbol, mass in elements.items():
        name(symbol, 'element')
        positive(mass, symbol)
    species, compositions, masses = [], [], []
    if not isinstance(document['species'], list) or not document['species']:
        raise MechanismError('species must be a nonempty list')
    for entry in document['species']:
        mapping(entry, ('name', 'composition'), 'species')
        symbol = name(entry['name'], 'species name')
        if symbol in species:
            raise MechanismError(f'duplicate species: {symbol}')
        composition = entry['composition']
        if not isinstance(composition, dict) or not composition or set(composition)-set(elements):
            raise MechanismError(f'{symbol}: invalid/unknown elements')
        if any(type(n) is not int or n <= 0 for n in composition.values()):
            raise MechanismError(f'{symbol}: atom counts must be positive integers')
        row = tuple(composition.get(e, 0) for e in elements)
        mass = sum(n*m for n, m in zip(row, elements.values()))
        positive(mass, symbol+' molar mass')
        species.append(symbol)
        compositions.append(row)
        masses.append(mass)
    if not isinstance(document['reactions'], list):
        raise MechanismError('reactions must be a list (empty is nonreacting)')
    ids, reactants, products, reversible = [], [], [], []
    for reaction in document['reactions']:
        mapping(reaction, ('id', 'reactants', 'products', 'reversible'), 'reaction')
        rid = name(reaction['id'], 'reaction id')
        if rid in ids:
            raise MechanismError(f'duplicate reaction id: {rid}')
        if type(reaction['reversible']) is not bool:
            raise MechanismError(f'{rid}: reversible must be boolean')
        sides = []
        for key in ('reactants', 'products'):
            side = reaction[key]
            if not isinstance(side, dict) or not side or set(side)-set(species):
                raise MechanismError(f'{rid}: empty side or unknown species in {key}')
            for s, coefficient in side.items():
                positive(coefficient, f'{rid}/{s}')
            sides.append(tuple(Fraction(str(side.get(s, 0))) for s in species))
        net = tuple(p-r for p, r in zip(sides[1], sides[0]))
        if not any(net):
            raise MechanismError(f'{rid}: reaction has no net change')
        for i, element in enumerate(elements):
            if sum(n*row[i] for n, row in zip(net, compositions)) != 0:
                raise MechanismError(f'{rid}: element {element} is not conserved')
        ids.append(rid)
        reactants.append(sides[0])
        products.append(sides[1])
        reversible.append(reaction['reversible'])
    return Mechanism(tuple(elements), tuple(species), tuple(masses), tuple(compositions),
                     tuple(ids), tuple(reactants), tuple(products), tuple(reversible))
