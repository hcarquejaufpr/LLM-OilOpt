# P1: controlled input for mathematical formulation

## What this package is

This is an **exploratory Open Math OS pilot** based on the same *problem input* used in the original confirmatory P1 experiment. It is not a new confirmatory run, a physical reservoir model, or a proof of simulator validity.

P1 has eight injectors, of which four are fixed and four (`INJECT3`, `INJECT6`, `INJECT7`, `INJECT8`) are optimization decisions. Its proxy is an additive PWL hinge model of accumulated oil production for a 3600-day horizon. All quantities and numerical coefficients are in `problem.json`.

## What an independent agent may see

- This document
- [problem.json](problem.json)
- [Audit instructions](../../experiments/openmath/P1_PILOT_INSTRUCTIONS.md)

Reference solutions, confirmatory outcomes, original generated formulations and private evaluators are **not** part of this input package.

## Audit boundaries

- Formulate the specified optimization problem independently; preserve bounds and capacity on **requested** injection rates.
- Inspect whether the chosen MILP/PWL representation is mathematically equivalent to the stated hinge definition throughout the declared bounds.
- Identify whether binary segment selectors or native SOS2 constraints are needed for the proposed formulation. Do not treat an unrestricted convex-hull relaxation as automatically equivalent.
- Distinguish a mathematical optimum of the PWL proxy from realized rates and production predicted by OPM Flow.
- WBHP 420 bar is a **separate physical validation criterion**, not a MILP constraint provided in this problem.
- Do not invent additional physical relationships or reservoir data.

## Source traceability

Source repository: `hcarquejaufpr/Pesquisa_Mestrado` (private), branch `main` at selection time, path:
`05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/runs/CF01-A/shared/problema.json`.

Original Git blob SHA: `e1acdb6b4f31034f5fef922cf22a24ecccd44d46`. Embedded `source_sha256` is preserved as supplied; its meaning should be checked against the original campaign manifest before a public release. The copy is only JSON pretty-printed, without changing field values.

**Public release status:** not approved yet. Review source provenance, model/data redistribution rights, and the associated paper before changing visibility.
