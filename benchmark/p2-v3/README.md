# P2 v3 frozen benchmark

Versioned benchmark package for the Revista Mundi article.

This directory is separate from the historical P1 benchmark. It contains the frozen mathematical problem, structured instructions, semantic rubric, reference solver and reference results, 16 selected run artifacts, blind semantic scorecards, and aggregate outcomes from eight paired blocks (16 sessions). The explicit M1 block matrix is `ANALISE_M1_M4_P2_V3/M1_BLOCK_AUDIT.csv`; regenerate it and recalculate the exact McNemar test with `python3 ANALISE_M1_M4_P2_V3/audit_m1_mcnemar.py`. M1: A 7/8, B 8/8, exact two-sided McNemar p=1.0. This does not establish statistical superiority or equivalence. M2–M4 are exploratory. C01's 0/5 execution marker follows a serialization incident and must not be read as five independent mathematical failures.

The `evaluator_private` name reflects its role during the completed campaign: these files were withheld from the LLMs. Their inclusion here is for subsequent reproduction, not a new blinded evaluation.

The reference solver depends on NumPy and SciPy (HiGHS through scipy.optimize.milp). It resolves its public JSON relative to this directory. Read its command-line options before running. No solver, OPM or LLM campaign was executed while preparing this package.

Physical validation aggregates are in `../../validation/p2-v3/`. These are post-campaign OPM results for frozen P1–P5 policies, not LLM feedback. Zero MILP gap does not establish a global physical reservoir optimum.

## Scope and limitations

- The repository is released under MIT; third-party reservoir decks and papers are excluded.
- Complete the reproduction environment and analysis entrypoints; this selection alone does not reproduce all paper claims end to end.
- Reconcile the manuscript's E3 RMSE provenance before publishing associated development claims.
- Review source-local paths retained in the C01 traceback; do not silently alter historical evidence.
- The published article should cite the DOI for the version used; the Zenodo concept DOI resolves to the latest version.

The original reservoir deck, third-party papers, manuscript drafts, credentials and full machine logs are not included. File provenance and exact SHA-256 hashes are in `../../validation/p2-v3/SHA256_MANIFEST.json`.
