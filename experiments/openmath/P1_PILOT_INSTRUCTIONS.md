# Open Math OS — P1 exploratory pilot

You are investigating a **controlled mathematical formulation benchmark**, not auditing a supplied LLM solution. Read `benchmark/p1/problem.json` and `benchmark/p1/README.md`.

1. Inventory every given parameter, decision variable, index, bound, unit and constraint. Note missing specifications without inventing them.
2. Express the objective and constraints precisely, respecting the PWL hinge formula given in the input.
3. Propose at least two equivalent, solver-compatible ways to formulate the PWL terms, **if equivalence can be demonstrated**. Carefully distinguish exact formulations from relaxations.
4. Perform algebraic and numerical cross-checks at selected points, including the hinge and both sides of it. Execute any available independent tools and record commands, versions and output.
5. Solve only if an appropriate installed solver is available. Identify objective, control decisions, status, gap and solver details. Never report an unexecuted result as computed.
6. Classify every conclusion as mechanically verified, reasoned but unverified, or unresolved. Do not claim physical validation without OPM Flow execution.
7. Save a standalone report under `experiments/openmath/reports/` without changing original inputs.

**Blinding:** Do not use the original confirmatory run outputs, generated LLM solutions, ground-truth optimum, or external/private evaluator in this first pilot. This is an exploratory exercise, not a preregistered confirmatory comparison.

**Environment:** Work in a separate directory/worktree if necessary. Do not delete or overwrite any existing Kakeya or other project files. If GitHub access is unavailable, request the input files rather than requesting an account token.
