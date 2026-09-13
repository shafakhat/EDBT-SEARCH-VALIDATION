# EDBT Search-Algorithm Benchmark — Code and Data

Experimental code, raw per-trial data, and a single-file analysis pipeline for
a benchmark study of twelve search algorithms under the EDBT (Extended
Dynamic Behavioural Tree) measurement framework. The repository is
self-contained and fully deterministic: rerunning the experiment reproduces
the data files byte-for-byte, and rerunning the pipeline regenerates every
table and figure from them.

## What is measured

The EDBT framework extends the classical tree-search formalism with three
online-computable components — a **complexity profile** γ, a
**pruning-accuracy** component Δ, and a **convergence operator** Ω. This
repository implements the measurement instrumentation and applies it to every
algorithm on every trial.

Reference for the quantities computed by the code (columns of the tables and
figures the pipeline emits):

| Symbol | Meaning (as computed by this code) |
|---|---|
| N̄, σ | Mean and standard deviation of expanded nodes over 30 trials. |
| U | Useful nodes — nodes that end up on the delivered solution. |
| SR | Solution rate: share of trials that returned a solution. |
| δ_s | Local fidelity ratio of each pruning step, computed from the retained vs. full candidate successor set. 1 = no accuracy discarded. |
| Δ_d | Depth-weighted pruning accuracy of the delivered search (successful runs). |
| q | Solution-quality factor: 1 for searches with an optimality guarantee, 1/w for WA*(w). |
| BEI | Behavioural Efficiency Index = U·q / (N̄·log_b(N̄+1)). |
| ABEI | Accuracy-Aware BEI = BEI × Δ_d. |
| β̂, α̂ | Log-linear node-growth fit N̄ = α̂·b^(β̂·d); R² = fit quality. |
| ACP | Accuracy–Complexity Product bound N̄·Δ_d ≥ γ⁻(β̂, d) = b^(β̂·d); the reported ratio is ≥ 1 when the bound holds. |
| PCCT | Pruning Correctness Certificate Test: δ̄_s ≥ 1/(1+ε_h) with ε_h the algorithm's bounded-suboptimality slack (∞ for GBFS). |
| C̄/C* | Delivered cost over Dijkstra's optimum. |
| MAPE | Held-out error of the β̂ model (train d ∈ {4,6,8}, test d ∈ {10,12}). |

## Experimental setup

- Weighted obstacle grids, size 4d × 4d, obstacle rate 0.15; edge cost is
  the weight of the cell entered (asymmetric).
- Depths d ∈ {4, 6, 8, 10, 12}; 30 trials per (algorithm, depth)
  configuration; deterministic seed 100·d + trial.
- Twelve algorithms:
  - graph searches — BFS, Dijkstra, A*, WA*(1.5), WA*(2.0), GBFS, BiDijkstra;
  - tree searches — RBFS (node cap 2×10⁶), IDA* (node cap 10⁷),
    Alpha–Beta (minimax tree instance);
  - stress algorithms — NoisyPrune(p), an A*-variant that drops a child with
    fixed probability p and never re-relaxes it, for p = 0.3 and p = 0.6.
- BiDijkstra is additionally measured by a dedicated 150-trial re-run driver
  on the same grids and seeds; the pipeline merges those rows automatically.

All result tables (per-depth results, β̂ fits, ACP ratios, PCCT
certificates, BEI/ABEI, MAPE, solution rates, realised suboptimality, and
paired Wilcoxon tests with Holm correction) and all figures are **generated
into `out/`** by the pipeline — they are not checked into the repository.

## Repository contents

| File | Description |
|---|---|
| `edbt_pipeline.py` | Single-file analysis pipeline: reads the two CSVs, cross-validates against the run log, writes all LaTeX tables, number macros, and figures to `out/` |
| `UnifiedTreeExperimentL.java` | The experiment: all twelve algorithms, seeded grid generation, online metric computation, per-trial CSV dump, end-of-run cross-check |
| `BiDRerun.java` | Dedicated 150-trial BiDijkstra re-run driver (same grids, same seeds) |
| `edbt_pertrial.csv` | Raw per-trial records, 1,800 rows |
| `edbt_bidij_rerun.csv` | Raw per-trial records for the BiDijkstra re-run, 150 rows |
| `results_full_run.txt` | Console log of the full experiment (used for cross-validation only) |
| `README.md` | This file |

## Requirements

- **Java** 8 or newer (only to regenerate the CSVs from scratch).
- **Python** 3.8 or newer with `numpy`, `pandas`, `scipy`, `matplotlib`.

## Usage

Regenerate all tables and figures from the shipped data (takes seconds):

```bash
python3 edbt_pipeline.py
```

Outputs land in `out/` (`tables*.tex`, `numbers.tex`, `figures/` as
PDF + PNG), and the run ends with `all cross-validation checks PASS`,
confirming the recomputed statistics match the experiment's console log.

Regenerate the raw data from scratch:

```bash
javac UnifiedTreeExperimentL.java BiDRerun.java
java UnifiedTreeExperimentL      # writes edbt_pertrial.csv
java BiDRerun                    # writes edbt_bidij_rerun.csv
python3 edbt_pipeline.py
```

`java UnifiedTreeExperimentL --fast` skips the two exponential tree searches
(RBFS, IDA*) for a quick smoke test. On machines with a display, a results
window opens after the run; the CSV and console output are identical either
way.

## Notes

- Figures render in grayscale only.
- The pipeline prints a note whenever the BiDijkstra re-run rows replace the
  main-run rows.

## Citation

If you use this code or data, please cite the accompanying paper. The
citation will be added here once the paper is published.
