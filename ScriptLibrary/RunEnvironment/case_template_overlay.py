"""Compose leaf overrides while preserving all base template options/comments."""
import json
import re


def apply_template_overrides(text, overrides):
    lines = text.splitlines()
    positions = {}
    stack = []
    for index, line in enumerate(lines):
        match = re.match(r"^( *)([A-Za-z_][A-Za-z_0-9]*):(?:\s|$)", line)
        if not match:
            continue
        indent, key = len(match[1]), match[2]
        while stack and stack[-1][0] >= indent:
            stack.pop()
        path = '.'.join([item[1] for item in stack] + [key])
        if path in positions:
            raise ValueError(f'duplicate template key: {path}')
        positions[path] = (index, indent)
        stack.append((indent, key))
    for path, value in overrides.items():
        if not isinstance(path, str) or not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*(\.[A-Za-z_][A-Za-z_0-9]*)*', path):
            raise ValueError(f'invalid template override path: {path}')
        encoded = json.dumps(value, ensure_ascii=False, allow_nan=False)
        if path not in positions:
            if '.' in path:
                raise ValueError(f'template override key not found: {path}')
            lines.append(f'{path}: {encoded}')
            continue
        index, indent = positions[path]
        if any(other.startswith(path + '.') for other in positions):
            raise ValueError(f'override individual leaves, not section: {path}')
        comment = lines[index].partition(' #')[2]
        lines[index] = ' ' * indent + path.split('.')[-1] + ': ' + encoded
        if comment:
            lines[index] += ' #' + comment
    return '\n'.join(lines) + '\n'
