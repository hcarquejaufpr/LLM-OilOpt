# LLM-OilOpt

**Reproducible research on LLM-assisted mathematical formulation for oilfield optimization.**

This repository accompanies an academic study comparing natural-language prompting with a structured SDD/TDD workflow for mathematical optimization models. The research context includes MILP, piecewise-linear (PWL) approximations, SOS2 constraints, and reservoir simulation with the Egg Model / OPM Flow.

## Status

The [P2 v3 benchmark](benchmark/p2-v3/README.md) contains frozen problem data, selected run evidence from 8 paired blocks (16 sessions), blinded semantic scorecards, reference results, and [post-campaign physical validation](validation/p2-v3/). The [P1 input](benchmark/p1/README.md) is a separate earlier exploratory pilot input. P1 results must not be combined with P2 v3 outcomes.

## Planned organization

- [benchmark](benchmark/README.md): problem specifications, authorized input data, and reference formulations.
- [experiments](experiments/README.md): exploratory and controlled evaluations across tools, including Open Math OS.
- [validation](validation/README.md): numerical checks, audit protocols, and reproducibility evidence.
- [docs](docs/README.md): study methodology, publication references, and release documentation.

## Research boundaries

The P2 v3 experiment is complete and frozen. Exploratory Open Math OS tests must be identified as such and must not be retroactively included in its paired comparison.

Only files checked for provenance, redistribution permissions, confidentiality, and completeness should be copied here. All experimental claims require reproducible evidence.

## Citation and licensing

This repository is distributed under the [MIT License](LICENSE). Cite the versioned GitHub release and its Zenodo DOI once the archive is published. The source reservoir deck and third-party papers are not distributed here.
