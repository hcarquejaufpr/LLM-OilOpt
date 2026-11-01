# P2 v3 frozen benchmark

Preparation snapshot for the Revista Mundi article; not yet a published reproduction release.

This directory is separate from the historical P1 benchmark. It contains the frozen mathematical problem, structured instructions, semantic rubric, reference solver and reference results, 16 selected run artifacts, blind semantic scorecards, and aggregate outcomes from eight paired blocks (16 sessions). M1: A 7/8, B 8/8, exact McNemar p=1.0. This does not establish statistical superiority or equivalence. M2–M4 are exploratory. C01's 0/5 execution marker follows a serialization incident and must not be read as five independent mathematical failures.

The `evaluator_private` name reflects its role during the completed campaign: these files were withheld from the LLMs. Their inclusion here is for subsequent reproduction, not a new blinded evaluation.

The reference solver depends on NumPy and SciPy (HiGHS through scipy.optimize.milp). It resolves its public JSON relative to this directory. Read its command-line options before running. No solver, OPM or LLM campaign was executed while preparing this package.

Physical validation aggregates are in `../../validation/p2-v3/`. These are post-campaign OPM results for frozen P1–P5 policies, not LLM feedback. Zero MILP gap does not establish a global physical reservoir optimum.

## Before the public release

- Decide redistribution licenses with the authors; no license is granted by this preparation note.
- Complete the reproduction environment and analysis entrypoints; this selection alone does not reproduce all paper claims end to end.
- Reconcile the manuscript's E3 RMSE provenance before publishing associated development claims.
- Review source-local paths retained in the C01 traceback; do not silently alter historical evidence.
- Add author ORCIDs and Zenodo metadata, reserve the DOI and freeze a release only after package verification.

The original reservoir deck, third-party papers, manuscript drafts, credentials and full machine logs are not included. File provenance and exact SHA-256 hashes are in `../../validation/p2-v3/SHA256_MANIFEST.json`.
