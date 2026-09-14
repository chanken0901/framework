#!/usr/bin/env python3
"""Validate an internal topology YAML, without simulating chemical reactions."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from reactingflow import validate_topology, MechanismError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    args = parser.parse_args()
    try:
        import yaml
        class UniqueLoader(yaml.SafeLoader):
            pass
        def unique_mapping(loader, node):
            result = {}
            for key_node, value_node in node.value:
                key = loader.construct_object(key_node)
                if not isinstance(key, str) or key in result:
                    raise MechanismError(f'duplicate/non-string YAML key: {key!r}')
                result[key] = loader.construct_object(value_node)
            return result
        UniqueLoader.add_constructor('tag:yaml.org,2002:map', unique_mapping)
        with args.input.open(encoding='utf-8') as f:
            model = validate_topology(yaml.load(f, Loader=UniqueLoader))
    except ImportError:
        parser.exit(1, 'Install PyYAML: python -m pip install PyYAML\n')
    except (OSError, ValueError, yaml.YAMLError) as exc:
        parser.exit(1, f'[ERROR] {exc}\n')
    print(f'[OK] topology: {len(model.species)} species, {len(model.reaction_ids)} reactions; elements conserved')
    print('Kinetics/thermodynamics/CFD are NOT evaluated by this foundation validator.')


if __name__ == '__main__':
    main()
