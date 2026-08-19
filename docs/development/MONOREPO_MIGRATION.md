# Monorepo migration record

- Migration date: 2026-08-19
- Base repository: https://github.com/chanken0901/framework.git
- Base branch: main
- Migration branch: chore/monorepo-migration
- FrameWork source commit: 8acaa476aab60e9f7c2f24a48c4aa43284ce30e0
- ScriptLibrary source commit: 2960c9efb5ac1de41d97c1a01ff02835a474b880
- SolverLibrary source commit: f525d15d8ecf1b34bfd470df589520bfbae90f74
- Safety tag: pre-monorepo-2026-08-19

ScriptLibrary and SolverLibrary were imported as merge parents with their
complete commit graphs, then placed below their prefixes with git read-tree.
The former standalone repositories must remain available until the monorepo
has passed review and a release tag has been created.
