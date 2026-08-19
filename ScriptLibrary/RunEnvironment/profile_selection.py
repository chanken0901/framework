"""Select a staged solver profile from editable case requirements."""

from __future__ import annotations

from typing import Any


class ProfileSelectionError(RuntimeError):
    """Raised when no staged solver profile can satisfy a case."""


def _canonical_selector(value: Any) -> str:
    return str(value).strip().lower().replace("-", "_").replace(" ", "_")


def case_profile_requirements(
    case: dict[str, Any], model: str
) -> set[str]:
    """Return build capabilities required by the editable case design."""

    requirements: set[str] = set()
    if model != "nse":
        return requirements

    flow = case.get("flow", {})
    if not isinstance(flow, dict):
        raise ProfileSelectionError("case flow section must be a YAML mapping")
    flow_type = _canonical_selector(flow.get("type", ""))
    if flow_type in {
        "hit",
        "hit_spectral",
        "homogeneous_isotropic_turbulence",
    }:
        requirements.add("hit_spectral")

    forcing = case.get("forcing", {})
    if not isinstance(forcing, dict):
        raise ProfileSelectionError("case forcing section must be a YAML mapping")
    forcing_type = _canonical_selector(
        forcing.get("type", forcing.get("scheme", "none"))
    )
    if forcing_type not in {"", "none"}:
        requirements.add("forcing_fft")
    return requirements


def profile_capabilities(
    manifest: dict[str, Any], profile_name: str
) -> set[str]:
    profiles = manifest.get("profiles")
    if not isinstance(profiles, dict) or not isinstance(
        profiles.get(profile_name), dict
    ):
        raise ProfileSelectionError(f"unknown solver profile: {profile_name}")
    profile = profiles[profile_name]
    raw_capabilities = profile.get("capabilities", [])
    if not isinstance(raw_capabilities, list):
        raise ProfileSelectionError(
            f"solver profile {profile_name!r} capabilities must be a list"
        )
    capabilities = {str(value) for value in raw_capabilities}

    # Schema-version-1 manifests did not declare capabilities explicitly.
    # Infer the NSE FFT capabilities so older generated environments return a
    # precise one-time regeneration request instead of selecting an invalid build.
    if str(manifest.get("model", "")).lower() == "nse":
        cmake = profile.get("cmake", {})
        if not isinstance(cmake, dict):
            raise ProfileSelectionError(
                f"solver profile {profile_name!r} cmake must be a mapping"
            )
        if str(cmake.get("NSE_INIT_FFT_BACKEND", "none")).lower() != "none":
            capabilities.add("hit_spectral")
        if str(cmake.get("NSE_FORCING_FFT_BACKEND", "none")).lower() != "none":
            capabilities.add("forcing_fft")
    return capabilities


def select_case_profile(
    case: dict[str, Any], manifest: dict[str, Any], lock: dict[str, Any]
) -> tuple[str, set[str]]:
    """Select a staged profile that satisfies the current case requirements."""

    base_profile = str(lock.get("profile", ""))
    if not base_profile:
        raise ProfileSelectionError("environment lock has no solver profile")
    available_value = lock.get("available_profiles", [base_profile])
    if not isinstance(available_value, list) or not available_value:
        raise ProfileSelectionError(
            "environment lock available_profiles must be a list"
        )
    available = [str(value) for value in available_value]
    if base_profile not in available:
        available.insert(0, base_profile)

    model = str(lock.get("model", manifest.get("model", ""))).lower()
    requirements = case_profile_requirements(case, model)
    solver = case.get("solver", {})
    if not isinstance(solver, dict):
        raise ProfileSelectionError("case solver section must be a YAML mapping")
    case_profile = solver.get("profile")
    if case_profile not in {None, ""}:
        requested = str(case_profile)
        if requested not in available:
            raise ProfileSelectionError(
                f"case solver.profile={requested!r} was not staged; regenerate "
                "the execution environment with that explicit profile"
            )
        candidates = [requested]
    elif bool(lock.get("profile_explicit", False)):
        candidates = [base_profile]
    else:
        candidates = available

    matches: list[tuple[int, int, str]] = []
    for index, profile_name in enumerate(candidates):
        capabilities = profile_capabilities(manifest, profile_name)
        if requirements.issubset(capabilities):
            # Prefer the least specialized satisfying profile. Candidate order
            # remains the tie breaker, so the baseline profile wins equal cases.
            matches.append((len(capabilities), index, profile_name))
    if matches:
        return min(matches)[2], requirements

    missing = ", ".join(sorted(requirements)) or "basic solver"
    raise ProfileSelectionError(
        f"no staged solver profile provides case capabilities: {missing}. "
        "Regenerate this execution environment from the current FrameWork; "
        "compatible profiles are now staged automatically"
    )
