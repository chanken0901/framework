"""Load and resolve reusable environment-design choices."""

from __future__ import annotations

import copy
import json
import os
import re
from pathlib import Path
from typing import Any

from yaml_support import load_yaml


OPTION_CATEGORIES = (
    "source",
    "destination",
    "model",
    "target",
    "case",
    "execution",
    "scheduler",
    "archive",
)


# Compatibility aliases for schema-version-1 designs created before model and
# parallel selection were separated. These aliases are intentionally not
# exposed by the option catalog: new designs select ``nse`` or ``gpe`` and use
# the parallel section. The profile override preserves the exact old build.
LEGACY_MODEL_SELECTIONS = {
    "nse_cpu_mpi": ("nse", "cpu_mpi"),
    "nse_cpu_mpi_2decomp_fftw": ("nse", "cpu_mpi_2decomp_fftw"),
    "nse_cuda_single": ("nse", "cuda_single"),
    "gpe_cpu_serial_dft": ("gpe", "cpu_serial_dft"),
    "gpe_cpu_serial_fftw": ("gpe", "cpu_serial_fftw"),
    "gpe_cpu_mpi_dft": ("gpe", "cpu_mpi_dft"),
    "gpe_cpu_mpi_fftw": ("gpe", "cpu_mpi_fftw"),
    "gpe_cuda_single": ("gpe", "cuda_single"),
    "gpe_cuda_mpi_cufftmp": ("gpe", "cuda_mpi_cufftmp"),
}


class OptionCatalogError(RuntimeError):
    """Raised when an environment option catalog is invalid."""


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise OptionCatalogError(f"{label} must be a YAML mapping")
    return value


def _sequence(value: Any, label: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise OptionCatalogError(f"{label} must be a YAML list")
    return value


def _catalog_path(value: str, design_path: Path) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(value))
    unresolved = re.search(r"\$\{[^}]+\}|%[^%]+%", expanded)
    if unresolved:
        raise OptionCatalogError(
            "undefined environment variable in option catalog path: "
            f"{unresolved.group(0)}"
        )
    path = Path(expanded)
    if not path.is_absolute():
        path = design_path.parent / path
    return path.resolve()


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def load_option_catalogs(
    design: dict[str, Any],
    design_path: Path,
    builtin_path: Path,
) -> dict[str, Any]:
    """Merge the built-in catalog and optional design-local catalogs."""

    paths = [builtin_path.resolve()]
    for value in _sequence(design.get("option_catalogs"), "option_catalogs"):
        paths.append(_catalog_path(str(value), design_path))

    options: dict[str, dict[str, Any]] = {
        category: {} for category in OPTION_CATEGORIES
    }
    catalogs: list[dict[str, Any]] = []
    for path in paths:
        if not path.is_file():
            raise OptionCatalogError(f"option catalog not found: {path}")
        document = _mapping(load_yaml(path), f"option catalog {path}")
        if document.get("schema_version") != 1:
            raise OptionCatalogError(
                f"option catalog schema_version must be 1: {path}"
            )
        catalog_id = str(document.get("catalog_id") or path.stem)
        allow_override = bool(document.get("allow_override", False))
        catalog_options = _mapping(
            document.get("options"), f"option catalog {catalog_id}.options"
        )
        unknown = sorted(set(catalog_options) - set(OPTION_CATEGORIES))
        if unknown:
            raise OptionCatalogError(
                f"option catalog {catalog_id} has unsupported categories: {unknown}"
            )

        for category, category_value in catalog_options.items():
            category_options = _mapping(
                category_value, f"option catalog {catalog_id}.{category}"
            )
            for option_id, option_value in category_options.items():
                option = _mapping(
                    option_value,
                    f"option catalog {catalog_id}.{category}.{option_id}",
                )
                values = _mapping(
                    option.get("values"),
                    f"option catalog {catalog_id}.{category}.{option_id}.values",
                )
                option_key = str(option_id)
                if option_key in options[category] and not allow_override:
                    raise OptionCatalogError(
                        f"duplicate option {category}.{option_key}; "
                        "set allow_override: true in the later catalog to replace it"
                    )
                record = copy.deepcopy(option)
                record["values"] = copy.deepcopy(values)
                record["catalog_id"] = catalog_id
                record["catalog_path"] = str(path)
                options[category][option_key] = record

        catalogs.append(
            {
                "catalog_id": catalog_id,
                "path": str(path),
                "allow_override": allow_override,
            }
        )
    return {"catalogs": catalogs, "options": options}


def _lookup(document: dict[str, Any], dotted_key: str) -> Any:
    value: Any = document
    for key in dotted_key.split("."):
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def _validate_requirements(
    design: dict[str, Any],
    selections: dict[str, dict[str, Any]],
) -> None:
    for category, selected in selections.items():
        requirements = selected.get("requires", {})
        if requirements is None:
            continue
        for dotted_key, expected in _mapping(
            requirements, f"option {category}.{selected['id']}.requires"
        ).items():
            actual = _lookup(design, str(dotted_key))
            accepted = expected if isinstance(expected, list) else [expected]
            if actual not in accepted:
                raise OptionCatalogError(
                    f"option {category}.{selected['id']} requires "
                    f"{dotted_key}={expected!r}, but resolved value is {actual!r}"
                )


def resolve_design_options(
    design: dict[str, Any],
    design_path: Path,
    builtin_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Resolve selected catalog entries, then apply inline section overrides."""

    catalog = load_option_catalogs(design, design_path, builtin_path)
    select = _mapping(design.get("select", {}), "select")
    unknown = sorted(set(select) - set(OPTION_CATEGORIES))
    if unknown:
        raise OptionCatalogError(f"select has unsupported categories: {unknown}")

    resolved = copy.deepcopy(design)
    legacy_profile: str | None = None
    selected_model = select.get("model")
    if selected_model is not None:
        legacy = LEGACY_MODEL_SELECTIONS.get(str(selected_model))
        if legacy is not None:
            canonical_model, legacy_profile = legacy
            select = copy.deepcopy(select)
            select["model"] = canonical_model
            resolved["select"] = copy.deepcopy(select)
    selections: dict[str, dict[str, Any]] = {}
    for category, option_id_value in select.items():
        option_id = str(option_id_value)
        available = catalog["options"][category]
        if option_id not in available:
            raise OptionCatalogError(
                f"unknown option {category}.{option_id}; "
                f"available: {sorted(available)}"
            )
        option = available[option_id]
        inline = resolved.get(category, {})
        if inline is None:
            inline = {}
        inline_mapping = _mapping(inline, category)
        resolved[category] = _deep_merge(option["values"], inline_mapping)
        selections[category] = {
            "id": option_id,
            "catalog_id": option["catalog_id"],
            "catalog_path": option["catalog_path"],
            "requires": copy.deepcopy(option.get("requires", {})),
        }

    if legacy_profile is not None:
        solver = resolved.get("solver", {})
        if solver is None:
            solver = {}
        solver_mapping = _mapping(solver, "solver")
        configured_profile = solver_mapping.get("profile")
        if configured_profile not in {None, "", legacy_profile}:
            raise OptionCatalogError(
                "legacy select.model conflicts with solver.profile: "
                f"{configured_profile!r} != {legacy_profile!r}"
            )
        resolved["solver"] = _deep_merge(
            {"profile": legacy_profile}, solver_mapping
        )

    _validate_requirements(resolved, selections)
    state = {
        "catalogs": catalog["catalogs"],
        "selections": selections,
        "options": catalog["options"],
    }
    return resolved, state


def validate_design_selections(
    design: dict[str, Any],
    state: dict[str, Any],
) -> None:
    """Validate selection requirements after any command-line overrides."""

    _validate_requirements(design, state["selections"])


def catalog_as_json(
    catalog: dict[str, Any],
    category: str | None = None,
) -> str:
    categories = OPTION_CATEGORIES if category is None else (category,)
    payload = {
        "schema_version": 1,
        "catalogs": catalog["catalogs"],
        "options": {
            name: catalog["options"][name]
            for name in categories
        },
    }
    return json.dumps(payload, ensure_ascii=False, indent=2)


def catalog_as_text(
    catalog: dict[str, Any],
    category: str | None = None,
) -> str:
    categories = OPTION_CATEGORIES if category is None else (category,)
    lines = ["実行環境の選択肢カタログ"]
    lines.append(
        "読込カタログ: "
        + ", ".join(item["catalog_id"] for item in catalog["catalogs"])
    )
    for name in categories:
        lines.extend(["", f"[{name}]"])
        category_options = catalog["options"][name]
        if not category_options:
            lines.append("  (候補なし)")
            continue
        for option_id, option in category_options.items():
            label = str(option.get("label") or option_id)
            description = str(option.get("description") or "")
            lines.append(f"  {option_id}: {label}")
            if description:
                lines.append(f"    {description}")
            lines.append(
                "    設定値: "
                + json.dumps(option["values"], ensure_ascii=False, sort_keys=True)
            )
            requirements = option.get("requires")
            if requirements:
                lines.append(
                    "    選択条件: "
                    + json.dumps(requirements, ensure_ascii=False, sort_keys=True)
                )
            lines.append(f"    カタログ: {option['catalog_id']}")
    return "\n".join(lines)
