#!/usr/bin/env python3
"""Validate the test-only cuFFTMp forcing extension contract."""

from __future__ import annotations

import sys
from pathlib import Path

def main() -> int:
    contract_path = Path(sys.argv[1])
    contract = contract_path.read_text(encoding="utf-8")
    required = {
        "backend: cufftmp",
        "status: test_only",
        "host_gather_allowed: false",
        "state_location: distributed_device_resident",
        "forward_distributed_fft",
        "project_local_spectral_pencil",
        "allreduce_forcing_denominators",
        "inverse_distributed_fft",
        "add_local_momentum_rhs",
        "runtime_available: false",
    }
    missing = sorted(item for item in required if item not in contract)
    if missing:
        raise AssertionError(f"missing cuFFTMp contract entries: {missing}")
    print("cuFFTMp forcing contract test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
