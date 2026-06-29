# SolverLibrary

SolverLibrary is a shared solver development repository.

## Solvers

- NSE

## Recommended structure

```text
SolverLibrary/
├── NSE/
│   └── src/
├── GPE/
│   └── src/
└── Shared/
```

## Role

This repository stores reusable solver source code.

Research projects import solver code from this library into their own `src/` directory, then generate a project-specific Makefile.
