#!/usr/bin/env python3
"""Safe YAML loader with a dependency-free fallback for NSE build files."""

from __future__ import annotations

import ast
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any


class YamlFormatError(ValueError):
    """Raised when a YAML design file cannot be parsed safely."""


@dataclass(frozen=True)
class _Token:
    indent: int
    text: str
    line: int


def _without_comment(raw: str) -> str:
    quote = ""
    escaped = False
    for index, char in enumerate(raw):
        if escaped:
            escaped = False
            continue
        if char == "\\" and quote == '"':
            escaped = True
            continue
        if quote:
            if char == quote:
                quote = ""
            continue
        if char in {"'", '"'}:
            quote = char
        elif char == "#" and (index == 0 or raw[index - 1].isspace()):
            return raw[:index]
    return raw


def _tokenize(text: str) -> list[_Token]:
    tokens: list[_Token] = []
    for line_number, raw in enumerate(text.splitlines(), start=1):
        leading = raw[: len(raw) - len(raw.lstrip())]
        if "\t" in leading:
            raise YamlFormatError(f"line {line_number}: tabs are not allowed for indentation")
        clean = _without_comment(raw).rstrip()
        if not clean.strip():
            continue
        indent = len(clean) - len(clean.lstrip(" "))
        tokens.append(_Token(indent, clean.strip(), line_number))
    return tokens


def _split_mapping(text: str, line: int) -> tuple[str, str]:
    quote = ""
    depth = 0
    escaped = False
    for index, char in enumerate(text):
        if escaped:
            escaped = False
            continue
        if char == "\\" and quote == '"':
            escaped = True
            continue
        if quote:
            if char == quote:
                quote = ""
            continue
        if char in {"'", '"'}:
            quote = char
        elif char in "[{(":
            depth += 1
        elif char in "]})":
            depth -= 1
        elif char == ":" and depth == 0:
            key = text[:index].strip().strip("'\"")
            if not key:
                raise YamlFormatError(f"line {line}: empty mapping key")
            return key, text[index + 1 :].strip()
    raise YamlFormatError(f"line {line}: expected key: value")


def _split_inline(text: str, line: int) -> list[str]:
    result: list[str] = []
    quote = ""
    depth = 0
    start = 0
    for index, char in enumerate(text):
        if quote:
            if char == quote and (index == 0 or text[index - 1] != "\\"):
                quote = ""
            continue
        if char in {"'", '"'}:
            quote = char
        elif char in "[{(":
            depth += 1
        elif char in "]})":
            depth -= 1
        elif char == "," and depth == 0:
            result.append(text[start:index].strip())
            start = index + 1
    if quote or depth != 0:
        raise YamlFormatError(f"line {line}: unterminated inline collection")
    tail = text[start:].strip()
    if tail:
        result.append(tail)
    return result


_INTEGER = re.compile(r"^[+-]?\d+$")
_FLOAT = re.compile(r"^[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?$")


def _scalar(text: str, line: int) -> Any:
    lower = text.lower()
    if lower in {"null", "~"}:
        return None
    if lower in {"true", "yes", "on"}:
        return True
    if lower in {"false", "no", "off"}:
        return False
    if text.startswith(("'", '"')):
        try:
            return ast.literal_eval(text)
        except (SyntaxError, ValueError) as exc:
            raise YamlFormatError(f"line {line}: invalid quoted scalar") from exc
    if text.startswith("[") and text.endswith("]"):
        body = text[1:-1].strip()
        return [] if not body else [_scalar(item, line) for item in _split_inline(body, line)]
    if text.startswith("{") and text.endswith("}"):
        body = text[1:-1].strip()
        result: dict[str, Any] = {}
        for item in _split_inline(body, line) if body else []:
            key, value = _split_mapping(item, line)
            result[key] = _scalar(value, line)
        return result
    if _INTEGER.fullmatch(text):
        return int(text)
    if _FLOAT.fullmatch(text):
        return float(text)
    return text


class _SubsetParser:
    def __init__(self, tokens: list[_Token]):
        self.tokens = tokens
        self.index = 0

    def parse(self) -> Any:
        if not self.tokens:
            return {}
        if self.tokens[0].indent != 0:
            raise YamlFormatError("top-level YAML content must not be indented")
        value = self._block(0)
        if self.index != len(self.tokens):
            token = self.tokens[self.index]
            raise YamlFormatError(f"line {token.line}: unexpected indentation")
        return value

    def _block(self, indent: int) -> Any:
        if self.tokens[self.index].text.startswith("-"):
            return self._sequence(indent)
        return self._mapping(indent)

    def _sequence(self, indent: int) -> list[Any]:
        result: list[Any] = []
        while self.index < len(self.tokens):
            token = self.tokens[self.index]
            if token.indent != indent or not token.text.startswith("-"):
                break
            if token.text != "-" and not token.text.startswith("- "):
                raise YamlFormatError(f"line {token.line}: invalid sequence marker")
            item = token.text[1:].strip()
            self.index += 1
            if item:
                result.append(_scalar(item, token.line))
            elif self.index < len(self.tokens) and self.tokens[self.index].indent > indent:
                result.append(self._block(self.tokens[self.index].indent))
            else:
                result.append(None)
        return result

    def _mapping(self, indent: int) -> dict[str, Any]:
        result: dict[str, Any] = {}
        while self.index < len(self.tokens):
            token = self.tokens[self.index]
            if token.indent != indent or token.text.startswith("-"):
                break
            key, value_text = _split_mapping(token.text, token.line)
            if key in result:
                raise YamlFormatError(f"line {token.line}: duplicate key {key!r}")
            self.index += 1
            if value_text:
                result[key] = _scalar(value_text, token.line)
            elif self.index < len(self.tokens) and self.tokens[self.index].indent > indent:
                result[key] = self._block(self.tokens[self.index].indent)
            else:
                result[key] = None
        return result


def load_yaml(path: str | Path) -> Any:
    """Load UTF-8 YAML using PyYAML when available, then the safe subset parser."""
    source = Path(path)
    text = source.read_text(encoding="utf-8-sig")
    stripped = text.lstrip()
    if stripped.startswith(("{", "[")):
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            pass
    try:
        import yaml  # type: ignore
    except ImportError:
        return _SubsetParser(_tokenize(text)).parse()
    try:
        return yaml.safe_load(text)
    except yaml.YAMLError as exc:  # type: ignore[attr-defined]
        raise YamlFormatError(f"{source}: {exc}") from exc
