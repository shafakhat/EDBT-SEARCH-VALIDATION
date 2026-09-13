# EDBT: An Accuracy-Aware Extension of Behavioural Tree Search — Code, Data, and Results

This repository contains the **complete experimental code, raw data, analysis
pipeline, and all results** for the research paper *"EDBT: An Accuracy-Aware
Extension of Behavioural Tree Search"*. A formal citation will be added here
upon publication.

The classical five-tuple *T* = (V, E, Φ, Ψ, Π) of tree-based heuristic search
captures graph structure and pruning mechanics, but provides no
online-computable measure of how much solution quality a pruning decision
discards. The **Extended Dynamic Behavioural Tree (EDBT)** augments it with a
complexity profile **γ**, an online-computable pruning-accuracy component
**Δ**, and a convergence operator **Ω**. This repository implements that
framework and validates it empirically on **twelve search algorithms** — ten
classical searches plus two purpose-built stress algorithms (NoisyPrune with
p = 0.3 and p = 0.6) that violate the evaluation-ordered-pruning condition
with fixed probability — across five depths on weighted obstacle grids
(**1,800 runs** in total, plus a dedicated 150-trial BiDijkstra re-run).

Because the local accuracy ratio δ_s is measured from the retained and
candidate successor sets each algorithm actually produces, the measurements
are falsifiable in practice: **δ_s < 1 occurs not only under NoisyPrune and
imperfectly ordered Alpha–Beta pruning, but under every historical-relaxation
and tree-search operator tested, while δ_s ≡ 1 holds structurally only for
breadth-first search (BFS).**

---

## Key findings

1. **The Accuracy–Complexity Product (ACP) bound N̄·Δ_d ≥ γ⁻(β̂, d) holds in
   53 of 54 solved configurations.** Its sole violation — GBFS at d = 4
   (ratio 0.45) — shows the check is only as strong as its β̂ calibration.
2. **The Pruning Correctness Certificate Test (PCCT) is one-sided:** four
   algorithms pass (GBFS only vacuously) and eight fail, yet all failing
   optimal searches realise C̄/C* = 1.000 — a fail is inconclusive for
   solution quality, exactly as designed.
3. **The Accuracy-Aware Behavioural Efficiency Index (ABEI) changes the
   ranking:** it discounts GBFS's raw advantage by **8.7×** at d = 10 and
   elevates BFS to second place at d = 12. Paired Wilcoxon tests
   (Holm-corrected) confirm these rank differences are systematic
   (p < 0.001 in 12 of 14 tests).

---

## The EDBT framework in one page

An EDBT is the eight-tuple **T⁺ = (V, E, Φ, Ψ, Π, γ, Δ, Ω)** — the classical
five-tuple plus:

| Component | Meaning |
|---|---|
| **γ** (complexity profile) | γ⁻(β, d) = b^(β·d): a lower bound on node growth at depth *d* with branching factor *b* = 4. β = 1 is the unpruned tree; β = 1/2 is the alpha–beta leaf-evaluation point. |
| **Δ** (pruning accuracy) | How much solution quality the pruning decisions discard, measured online. |
| **Ω** (convergence operator) | Relates the accuracy demand ε to the depth at which the search converges. |

### The measured quantities

| Symbol | Definition (as computed by this code) |
|---|---|
| N̄, σ | Mean and standard deviation of expanded nodes over 30 trials. |
| U | Useful nodes — nodes that end up on the delivered solution. |
| SR | Solution rate: share of trials that returned a solution. |
| δ_s | Local fidelity ratio of each pruning step, computed from the retained vs. full candidate successor set (mean over steps and trials). 1 = no accuracy discarded. |
| Δ_d | Depth-weighted pruning accuracy of the delivered search (mean over successful runs). |
| q | Solution-quality factor: 1 for searches with an optimality guarantee, 1/w for WA*(w) (bounded suboptimality C ≤ w·C*). |
| **BEI** | Behavioural Efficiency Index = U·q / (N̄·log_b(N̄+1)): the useful-node share per unit of search depth. |
| **ABEI** | Accuracy-Aware BEI = BEI × Δ_d: the same share discounted by measured pruning accuracy. |
| β̂, α̂ | Log-linear fit N̄ = α̂·b^(β̂·d) across all five depths; R² is the fit quality. |
| **ACP** | Accuracy–Complexity Product bound N̄·Δ_d ≥ γ⁻(β̂, d); the reported ratio must be ≥ 1. β̂ is refit incrementally as each depth arrives. |
| **PCCT** | Pruning Correctness Certificate Test: passes if δ̄_s ≥ 1/(1+ε_h), where ε_h is the algorithm's heuristic-suboptimality bound (ε_h = w−1 for WA*(w), 0 for admissible searches, ∞ for GBFS). A pass is a machine-checked sufficient condition; a fail is inconclusive for solution quality. |
| C̄/C* | Realised suboptimality: mean delivered cost over Dijkstra's optimum (successful runs only). |

---

## Experimental setup

- **Domain:** weighted obstacle grids of size 4d × 4d with obstacle rate 0.15;
  edge cost = weight of the cell entered (asymmetric). Depth
  d ∈ {4, 6, 8, 10, 12} → grids from 16×16 to 48×48.
- **Trials:** 30 per (algorithm, depth) configuration; grids are generated
  from the deterministic seed 100·d + trial, so every run is exactly
  reproducible.
- **Algorithms (12):**
  - *Graph searches:* BFS, Dijkstra, A*, WA*(1.5), WA*(2.0), GBFS, BiDijkstra;
  - *Tree searches:* RBFS (node cap 2×10⁶), IDA* (node cap 10⁷), Alpha–Beta
    (solves a minimax tree instance, not the grid path problem);
  - *Stress algorithms:* NoisyPrune(p) — an A*-variant that drops a child with
    fixed probability p at expansion and never re-relaxes it (incomplete by
    construction), for p = 0.3 and p = 0.6.
- **Node caps** make the exponential tree searches budget-terminated rather
  than run forever; a capped run counts as unsuccessful (SR contribution 0).
- **BiDijkstra** is additionally measured in a dedicated 150-trial re-run
  (same grids, same seeds) that isolates its asymmetric backward relaxation;
  the analysis pipeline merges those rows automatically.

**Everything below is regenerated by `python3 edbt_pipeline.py` from the raw
per-trial CSVs in this repository.**

---

## Results

### Headline: ranking by ABEI at d = 12 (48×48 grids)

| # | Algorithm | N̄ | BEI | ABEI | Δ_d | SR (%) | C̄/C* |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | GBFS | 114.7 | 0.2697 | 0.0301 | 0.112 | 100 | 1.651 |
| 2 | BFS | 1970.3 | 0.0088 | 0.0088 | 1.000 | 100 | 1.468 |
| 3 | BiDijkstra | 1481.4 | 0.0123 | 0.0084 | 0.680 | 100 | 1.000 |
| 4 | Dijkstra | 1968.3 | 0.0089 | 0.0060 | 0.678 | 100 | 1.000 |
| 5 | A* | 1965.4 | 0.0089 | 0.0060 | 0.675 | 100 | 1.000 |
| 6 | NP(0.3) | 1888.1 | 0.0083 | 0.0048 | 0.579 | 87 | 1.194 |
| 7 | WA*(1.5) | 1958.6 | 0.0060 | 0.0047 | 0.779 | 100 | 1.000 |
| 8 | WA*(2.0) | 1930.0 | 0.0046 | 0.0038 | 0.835 | 100 | 1.000 |
| 9 | NP(0.6) | 1563.7 | 0.0050 | 0.0023 | 0.452 | 40 | 1.497 |
| 10 | AlphaBeta | 467448.9 | 0.0000 | 0.0000 | 1.000 | 100 | — |
| 11 | RBFS | 2000000.0 | 0.0000 | 0.0000 | — | 0 | — |
| 12 | IDA* | 10000020.9 | 0.0000 | 0.0000 | — | 0 | — |

GBFS delivers a solution using ~115 expansions where BFS needs ~1,970 — but
its solutions cost 1.65× the optimum and its accuracy discount is the largest
of any algorithm. "—" marks quantities undefined when SR = 0 (or not
applicable: Alpha–Beta solves the minimax problem).

### Per-depth results

SR — solution rate; Δ_d and C̄/C* are averaged over successful runs only;
δ_s over all runs. ACP: N̄·Δ_d ≥ γ⁻(β̂, d) with β̂ refit incrementally.

<details><summary><b>Depth d = 4 (16×16 grids)</b></summary>

| Algorithm | N̄ ± σ | SR (%) | Δ_d | δ_s | BEI | ABEI | ACP | C̄/C* |
|---|---|---:|---:|---:|---:|---:|:---:|---:|
| GBFS | 37.0 ± 5.1 | 100 | 0.194 | 0.994 | 0.3498 | 0.0678 | **✗** | 1.548 |
| BiDijkstra | 158.8 ± 16.5 | 100 | 0.677 | 0.947 | 0.0539 | 0.0365 | ✓ | 1.000 |
| BFS | 222.3 ± 5.3 | 100 | 1.000 | 1.000 | 0.0358 | 0.0358 | ✓ | 1.432 |
| A* | 216.9 ± 9.8 | 100 | 0.679 | 0.983 | 0.0372 | 0.0253 | ✓ | 1.000 |
| Dijkstra | 220.2 ± 6.5 | 100 | 0.678 | 0.964 | 0.0365 | 0.0248 | ✓ | 1.000 |
| WA*(1.5) | 212.3 ± 13.9 | 100 | 0.775 | 0.986 | 0.0254 | 0.0197 | ✓ | 1.000 |
| NP(0.3) | 212.9 ± 11.2 | 80 | 0.604 | 0.992 | 0.0308 | 0.0186 | ✓ | 1.169 |
| WA*(2.0) | 203.1 ± 20.5 | 100 | 0.836 | 0.989 | 0.0201 | 0.0168 | ✓ | 1.000 |
| NP(0.6) | 173.8 ± 64.8 | 50 | 0.440 | 0.981 | 0.0257 | 0.0113 | ✓ | 1.564 |
| AlphaBeta | 192.5 ± 37.1 | 100 | 1.000 | 0.809 | 0.0068 | 0.0068 | ✓ | — |
| RBFS | 47617.4 ± 16643.5 | 100 | 0.679 | 0.855 | 0.0001 | 0.0001 | ✓ | 1.000 |
| IDA* | 8636986.9 ± 2468729.1 | 30 | 0.692 | 0.989 | 0.0000 | 0.0000 | ✓ | 1.000 |

</details>

<details><summary><b>Depth d = 6 (24×24 grids)</b></summary>

| Algorithm | N̄ ± σ | SR (%) | Δ_d | δ_s | BEI | ABEI | ACP | C̄/C* |
|---|---|---:|---:|---:|---:|---:|:---:|---:|
| GBFS | 55.8 ± 4.0 | 100 | 0.152 | 0.997 | 0.3190 | 0.0485 | ✓ | 1.621 |
| BFS | 497.6 ± 8.8 | 100 | 1.000 | 1.000 | 0.0211 | 0.0211 | ✓ | 1.460 |
| BiDijkstra | 366.7 ± 29.9 | 100 | 0.689 | 0.961 | 0.0303 | 0.0209 | ✓ | 1.000 |
| A* | 491.3 ± 14.1 | 100 | 0.692 | 0.988 | 0.0215 | 0.0149 | ✓ | 1.000 |
| Dijkstra | 495.6 ± 9.6 | 100 | 0.691 | 0.974 | 0.0213 | 0.0147 | ✓ | 1.000 |
| NP(0.3) | 470.6 ± 87.8 | 83 | 0.598 | 0.995 | 0.0193 | 0.0116 | ✓ | 1.183 |
| WA*(1.5) | 483.8 ± 23.3 | 100 | 0.784 | 0.991 | 0.0146 | 0.0115 | ✓ | 1.000 |
| WA*(2.0) | 466.7 ± 40.4 | 100 | 0.841 | 0.992 | 0.0115 | 0.0096 | ✓ | 1.000 |
| NP(0.6) | 336.0 ± 202.2 | 43 | 0.460 | 0.979 | 0.0157 | 0.0072 | ✓ | 1.536 |
| AlphaBeta | 1589.5 ± 300.6 | 100 | 1.000 | 0.828 | 0.0008 | 0.0008 | ✓ | — |
| RBFS | 310305.7 ± 90763.7 | 100 | 0.692 | 0.856 | 0.0000 | 0.0000 | ✓ | 1.000 |
| IDA* | 10000020.9 ± 6.9 | 0 | — | 0.990 | 0.0000 | 0.0000 | n/a | — |

</details>

<details><summary><b>Depth d = 8 (32×32 grids)</b></summary>

| Algorithm | N̄ ± σ | SR (%) | Δ_d | δ_s | BEI | ABEI | ACP | C̄/C* |
|---|---|---:|---:|---:|---:|---:|:---:|---:|
| GBFS | 79.0 ± 14.5 | 100 | 0.123 | 0.993 | 0.2871 | 0.0355 | ✓ | 1.623 |
| BFS | 879.6 ± 8.9 | 100 | 1.000 | 1.000 | 0.0146 | 0.0146 | ✓ | 1.416 |
| BiDijkstra | 645.7 ± 41.9 | 100 | 0.693 | 0.971 | 0.0210 | 0.0146 | ✓ | 1.000 |
| A* | 874.7 ± 9.8 | 100 | 0.692 | 0.991 | 0.0148 | 0.0103 | ✓ | 1.000 |
| Dijkstra | 878.0 ± 8.9 | 100 | 0.692 | 0.981 | 0.0148 | 0.0102 | ✓ | 1.000 |
| WA*(1.5) | 866.4 ± 20.7 | 100 | 0.796 | 0.993 | 0.0100 | 0.0080 | ✓ | 1.000 |
| NP(0.3) | 869.2 ± 10.9 | 90 | 0.575 | 0.996 | 0.0136 | 0.0078 | ✓ | 1.179 |
| WA*(2.0) | 846.1 ± 48.2 | 100 | 0.853 | 0.994 | 0.0077 | 0.0066 | ✓ | 1.000 |
| NP(0.6) | 668.1 ± 332.9 | 40 | 0.432 | 0.985 | 0.0090 | 0.0039 | ✓ | 1.545 |
| AlphaBeta | 12965.5 ± 2574.6 | 100 | 1.000 | 0.828 | 0.0001 | 0.0001 | ✓ | — |
| RBFS | 1263091.2 ± 349423.8 | 93 | 0.694 | 0.857 | 0.0000 | 0.0000 | ✓ | 1.000 |
| IDA* | 10000020.2 ± 5.9 | 0 | — | 0.991 | 0.0000 | 0.0000 | n/a | — |

</details>

<details><summary><b>Depth d = 10 (40×40 grids)</b></summary>

| Algorithm | N̄ ± σ | SR (%) | Δ_d | δ_s | BEI | ABEI | ACP | C̄/C* |
|---|---|---:|---:|---:|---:|---:|:---:|---:|
| GBFS | 95.1 ± 9.2 | 100 | 0.114 | 0.997 | 0.2796 | 0.0320 | ✓ | 1.678 |
| BFS | 1366.7 ± 11.9 | 100 | 1.000 | 1.000 | 0.0111 | 0.0111 | ✓ | 1.495 |
| BiDijkstra | 1015.5 ± 67.5 | 100 | 0.690 | 0.975 | 0.0158 | 0.0109 | ✓ | 1.000 |
| A* | 1360.8 ± 14.8 | 100 | 0.692 | 0.993 | 0.0113 | 0.0078 | ✓ | 1.000 |
| Dijkstra | 1364.8 ± 12.9 | 100 | 0.690 | 0.984 | 0.0113 | 0.0078 | ✓ | 1.000 |
| NP(0.3) | 1347.9 ± 19.1 | 90 | 0.600 | 0.997 | 0.0104 | 0.0062 | ✓ | 1.177 |
| WA*(1.5) | 1353.9 ± 18.6 | 100 | 0.790 | 0.994 | 0.0076 | 0.0060 | ✓ | 1.000 |
| WA*(2.0) | 1325.2 ± 41.8 | 100 | 0.850 | 0.995 | 0.0058 | 0.0050 | ✓ | 1.000 |
| NP(0.6) | 998.1 ± 542.8 | 40 | 0.467 | 0.987 | 0.0067 | 0.0031 | ✓ | 1.502 |
| AlphaBeta | 83760.9 ± 15426.0 | 100 | 1.000 | 0.829 | 0.0000 | 0.0000 | ✓ | — |
| RBFS | 2000000.0 ± 0.0 | 0 | — | 0.866 | 0.0000 | 0.0000 | n/a | — |
| IDA* | 10000020.5 ± 5.7 | 0 | — | 0.992 | 0.0000 | 0.0000 | n/a | — |

</details>

**Depth d = 12 (48×48 grids)** — see the headline table above.

### Pruning-efficiency fits: N̄ = α̂·b^(β̂·d), b = 4

| Algorithm | β̂ | α̂ | R² |
|---|---:|---:|---:|
| IDA* | 0.0106 | 8636986.84 | 0.5000 |
| GBFS | 0.1008 | 23.15 | 0.9658 |
| BFS | 0.1938 | 89.12 | 0.9746 |
| Dijkstra | 0.1945 | 88.12 | 0.9741 |
| NP(0.3) | 0.1954 | 84.73 | 0.9718 |
| A* | 0.1957 | 86.40 | 0.9738 |
| WA*(1.5) | 0.1974 | 83.92 | 0.9738 |
| NP(0.6) | 0.1978 | 63.74 | 0.9863 |
| BiDijkstra | 0.1978 | 62.76 | 0.9741 |
| WA*(2.0) | 0.2001 | 79.21 | 0.9738 |
| RBFS | 0.3368 | 14200.45 | 0.8461 |
| AlphaBeta | 0.7053 | 4.38 | 0.9980 |

All grid algorithms cluster at β̂ ≈ 0.19–0.20 (R² ≥ 0.97) — far below the
unpruned β = 1 and well below the alpha–beta point β = 1/2; Alpha–Beta sits
at β̂ = 0.71 with the best fit of all (R² = 0.998). RBFS and IDA* saturate
their node caps, so their fits reflect the caps, not the algorithms.

### ACP bound ratio N̄·Δ_d / γ⁻(β̂, d) — values ≥ 1 satisfy the bound

| Algorithm | d=4 | d=6 | d=8 | d=10 | d=12 |
|---|---:|---:|---:|---:|---:|
| BFS | 13.89 | 44.35 | 56.17 | 67.42 | 78.35 |
| Dijkstra | 9.33 | 30.02 | 38.22 | 45.86 | 52.42 |
| A* | 9.20 | 29.23 | 37.22 | 44.87 | 51.13 |
| WA*(1.5) | 10.28 | 32.08 | 41.41 | 49.62 | 57.17 |
| WA*(2.0) | 10.61 | 32.31 | 41.60 | 50.18 | 57.78 |
| GBFS | **0.45** | 2.48 | 2.14 | 2.22 | 2.39 |
| BiDijkstra | 6.72 | 20.52 | 27.06 | 32.66 | 37.48 |
| RBFS | 2021.26 | 775.48 | 1246.57 | n/a | n/a |
| IDA* | 373373.52 | n/a | n/a | n/a | n/a |
| AlphaBeta | 12.03 | 2.82 | 2.86 | 3.23 | 3.75 |
| NP(0.3) | 8.03 | 26.06 | 29.99 | 37.35 | 42.38 |
| NP(0.6) | 4.77 | 21.38 | 19.53 | 24.01 | 26.34 |

53 of 54 solved configurations satisfy the bound. The single violation
(bold) is GBFS at d = 4, where β̂ is still at its default calibration value.
"n/a": SR = 0, so Δ_d and the ACP product are undefined.

### PCCT certificate test

| Algorithm | ε_h | δ_s^cert | δ̄_s^obs | PCCT |
|---|---:|---:|---:|---|
| BFS | 0.00 | 1.000 | 1.000 | pass |
| GBFS | ∞ | 0 | 0.995 | vacuous |
| NP(0.3) | 0.00 | 1.000 | 0.995 | fail |
| WA*(2.0) | 1.00 | 0.500 | 0.993 | pass |
| WA*(1.5) | 0.50 | 0.667 | 0.992 | pass |
| IDA* | 0.00 | 1.000 | 0.991 | fail |
| A* | 0.00 | 1.000 | 0.990 | fail |
| NP(0.6) | 0.00 | 1.000 | 0.985 | fail |
| Dijkstra | 0.00 | 1.000 | 0.978 | fail |
| BiDijkstra | 0.00 | 1.000 | 0.967 | fail |
| RBFS | 0.00 | 1.000 | 0.862 | fail |
| AlphaBeta | 0.00 | 1.000 | 0.825 | fail |

Four algorithms pass and eight fail — yet every failing optimal search still
delivers C̄/C* = 1.000, demonstrating that the certificate is a one-sided
sufficient condition whose failure is inconclusive for solution quality.

### ABEI vs BEI at d = 10 — the accuracy discount

| Algorithm | N̄ | Δ_d | BEI | ABEI | BEI/ABEI |
|---|---:|---:|---:|---:|---:|
| GBFS | 95.1 | 0.114 | 0.2796 | 0.0320 | **8.7×** |
| BFS | 1366.7 | 1.000 | 0.0111 | 0.0111 | 1.0× |
| BiDijkstra | 1015.5 | 0.690 | 0.0158 | 0.0109 | 1.4× |
| A* | 1360.8 | 0.692 | 0.0113 | 0.0078 | 1.4× |
| Dijkstra | 1364.8 | 0.690 | 0.0113 | 0.0078 | 1.4× |
| NP(0.3) | 1347.9 | 0.600 | 0.0104 | 0.0062 | 1.7× |
| WA*(1.5) | 1353.9 | 0.790 | 0.0076 | 0.0060 | 1.3× |
| WA*(2.0) | 1325.2 | 0.850 | 0.0058 | 0.0050 | 1.2× |
| NP(0.6) | 998.1 | 0.467 | 0.0067 | 0.0031 | 2.1× |
| AlphaBeta | 83760.9 | 1.000 | 0.0000 | 0.0000 | 1.0× |
| RBFS | 2000000.0 | — | 0.0000 | 0.0000 | — |
| IDA* | 10000020.5 | — | 0.0000 | 0.0000 | — |

### BEI and ABEI by depth

| Algorithm | d=4 BEI | d=4 ABEI | d=6 BEI | d=6 ABEI | d=8 BEI | d=8 ABEI | d=10 BEI | d=10 ABEI | d=12 BEI | d=12 ABEI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BFS | 0.0358 | 0.0358 | 0.0211 | 0.0211 | 0.0146 | 0.0146 | 0.0111 | 0.0111 | 0.0088 | 0.0088 |
| Dijkstra | 0.0365 | 0.0248 | 0.0213 | 0.0147 | 0.0148 | 0.0102 | 0.0113 | 0.0078 | 0.0089 | 0.0060 |
| A* | 0.0372 | 0.0253 | 0.0215 | 0.0149 | 0.0148 | 0.0103 | 0.0113 | 0.0078 | 0.0089 | 0.0060 |
| WA*(1.5) | 0.0254 | 0.0197 | 0.0146 | 0.0115 | 0.0100 | 0.0080 | 0.0076 | 0.0060 | 0.0060 | 0.0047 |
| WA*(2.0) | 0.0201 | 0.0168 | 0.0115 | 0.0096 | 0.0077 | 0.0066 | 0.0058 | 0.0050 | 0.0046 | 0.0038 |
| GBFS | 0.3498 | 0.0678 | 0.3190 | 0.0485 | 0.2871 | 0.0355 | 0.2796 | 0.0320 | 0.2697 | 0.0301 |
| BiDijkstra | 0.0539 | 0.0365 | 0.0303 | 0.0209 | 0.0210 | 0.0146 | 0.0158 | 0.0109 | 0.0123 | 0.0084 |
| RBFS | 0.0001 | 0.0001 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| IDA* | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| AlphaBeta | 0.0068 | 0.0068 | 0.0008 | 0.0008 | 0.0001 | 0.0001 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| NP(0.3) | 0.0308 | 0.0186 | 0.0193 | 0.0116 | 0.0136 | 0.0078 | 0.0104 | 0.0062 | 0.0083 | 0.0048 |
| NP(0.6) | 0.0257 | 0.0113 | 0.0157 | 0.0072 | 0.0090 | 0.0039 | 0.0067 | 0.0031 | 0.0050 | 0.0023 |

### Held-out prediction error of the β̂ model (train d ∈ {4,6,8}, test d ∈ {10,12})

| Algorithm | MAPE (%) |
|---|---:|
| IDA* | 14.45 |
| GBFS | 35.39 |
| NP(0.6) | 47.07 |
| AlphaBeta | 57.25 |
| BFS | 58.46 |
| Dijkstra | 59.41 |
| BiDijkstra | 59.80 |
| A* | 60.63 |
| WA*(1.5) | 61.22 |
| WA*(2.0) | 62.80 |
| NP(0.3) | 63.97 |
| RBFS | 981.86 |

RBFS and IDA* saturate their caps over the test depths, so their MAPE
reflects cap saturation.

### Solution rate (%) by depth

| Algorithm | d=4 | d=6 | d=8 | d=10 | d=12 |
|---|---:|---:|---:|---:|---:|
| BFS | 100 | 100 | 100 | 100 | 100 |
| Dijkstra | 100 | 100 | 100 | 100 | 100 |
| A* | 100 | 100 | 100 | 100 | 100 |
| WA*(1.5) | 100 | 100 | 100 | 100 | 100 |
| WA*(2.0) | 100 | 100 | 100 | 100 | 100 |
| GBFS | 100 | 100 | 100 | 100 | 100 |
| BiDijkstra | 100 | 100 | 100 | 100 | 100 |
| RBFS | 100 | 100 | 93 | 0 | 0 |
| IDA* | 30 | 0 | 0 | 0 | 0 |
| AlphaBeta | 100 | 100 | 100 | 100 | 100 |
| NP(0.3) | 80 | 83 | 90 | 90 | 87 |
| NP(0.6) | 50 | 43 | 40 | 40 | 40 |

### Realised suboptimality C̄/C* (C* from Dijkstra)

| Algorithm | d=4 | d=6 | d=8 | d=10 | d=12 |
|---|---:|---:|---:|---:|---:|
| BFS | 1.432 | 1.460 | 1.416 | 1.495 | 1.468 |
| Dijkstra | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| A* | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| WA*(1.5) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| WA*(2.0) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| GBFS | 1.548 | 1.621 | 1.623 | 1.678 | 1.651 |
| BiDijkstra | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| RBFS | 1.000 | 1.000 | 1.000 | — | — |
| IDA* | 1.000 | — | — | — | — |
| NP(0.3) | 1.169 | 1.183 | 1.179 | 1.177 | 1.194 |
| NP(0.6) | 1.564 | 1.536 | 1.545 | 1.502 | 1.497 |

BFS returns hop-minimal paths, which need not be cost-minimal. Alpha–Beta
solves the minimax problem on a separate tree instance and is not comparable.

### Statistical tests — paired Wilcoxon signed-rank (30 paired trials, two-sided, Holm-corrected across the 14 tests)

| Measure | Comparison (A vs B) | med(A) | med(B) | p_Holm | n |
|---|---|---:|---:|---:|---:|
| N̄, d=10 | A* vs Dijkstra | 1360.5 | 1367.5 | <0.001 | 30 |
| N̄, d=10 | BiDijkstra vs Dijkstra | 1007.5 | 1367.5 | <0.001 | 30 |
| N̄, d=10 | GBFS vs BFS | 92.5 | 1368.5 | <0.001 | 30 |
| ABEI, d=10 | BFS vs Dijkstra | 0.0111 | 0.0078 | <0.001 | 30 |
| ABEI, d=10 | BFS vs BiDijkstra | 0.0111 | 0.0106 | 0.452 | 30 |
| ABEI, d=10 | GBFS vs BFS | 0.0317 | 0.0111 | <0.001 | 30 |
| ABEI, d=10 | A* vs WA*(1.5) | 0.0078 | 0.0060 | <0.001 | 30 |
| N̄, d=12 | A* vs Dijkstra | 1967.0 | 1968.0 | <0.001 | 30 |
| N̄, d=12 | BiDijkstra vs Dijkstra | 1502.5 | 1968.0 | <0.001 | 30 |
| N̄, d=12 | GBFS vs BFS | 114.5 | 1970.0 | <0.001 | 30 |
| ABEI, d=12 | BFS vs Dijkstra | 0.0088 | 0.0060 | <0.001 | 30 |
| ABEI, d=12 | BFS vs BiDijkstra | 0.0088 | 0.0083 | 0.013 | 30 |
| ABEI, d=12 | GBFS vs BFS | 0.0324 | 0.0088 | <0.001 | 30 |
| ABEI, d=12 | A* vs WA*(1.5) | 0.0061 | 0.0046 | <0.001 | 30 |

---

## Figures

All figures are grayscale, overlap-checked, and regenerated by the pipeline
into `out/figures/` (PDF for print, PNG for preview; the PNGs shipped in
`figures/` are the ones shown below).

**Node growth vs depth** — graph searches (log scale) and the capped tree
searches:

![Node growth](figures/fig_nodes.png)

**Fidelity by depth** — the local accuracy ratio δ_s and depth-weighted
accuracy Δ_d for every algorithm:

![Fidelity](figures/fig_fidelity.png)

**ACP bound ratios** (log scale; ≥ 1 satisfies the bound; GBFS at d = 4 is
the sole violation):

![ACP](figures/fig_acp.png)

**BEI (open circles) vs ABEI (filled squares) at d = 12** — the accuracy
discount reorders the field:

![BEI vs ABEI](figures/fig_bei_abei.png)

**Solution rate by depth and realised suboptimality at d = 12:**

![Solution rate and cost](figures/fig_sr_cost.png)

---

## Repository contents

| File | Description |
|---|---|
| `edbt_pipeline.py` | Single-file analysis pipeline: reads the two CSVs, cross-validates against the run log, and regenerates all LaTeX tables, number macros, and figures into `out/` |
| `UnifiedTreeExperimentL.java` | The experiment itself: all twelve algorithms, grid generation (seeded, reproducible), online metric computation, per-trial CSV dump, and a cross-check against the paper's printed tables |
| `BiDRerun.java` | Dedicated 150-trial BiDijkstra re-run driver (same grids, same seeds) |
| `edbt_pertrial.csv` | Raw per-trial records, 1,800 rows (12 algorithms × 5 depths × 30 trials) |
| `edbt_bidij_rerun.csv` | Raw per-trial records for the dedicated BiDijkstra re-run, 150 rows |
| `results_full_run.txt` | Console log of the full Java experiment (used by the pipeline for cross-validation only) |
| `figures/` | The five result figures as PNG (regenerated by the pipeline as PDF+PNG into `out/figures/`) |
| `README.md` | This file |

## Requirements

- **Java** 8 or newer (only needed to regenerate the CSVs from scratch).
- **Python** 3.8 or newer with `numpy`, `pandas`, `scipy`, `matplotlib`
  (e.g. `pip install numpy pandas scipy matplotlib`).

## Reproducing the results

### From the shipped data (seconds)

```bash
python3 edbt_pipeline.py
```

This writes to `out/`: `tables.tex`, `tables_depth.tex`, `tables_theory.tex`,
`tables_abei.tex`, `tables_final.tex`, `numbers.tex`, and
`figures/fig_nodes`, `fig_fidelity`, `fig_acp`, `fig_bei_abei`, `fig_sr_cost`
(PDF + PNG). The pipeline re-derives every statistic from the CSVs and
cross-validates them against the Java experiment's printed console log
(`results_full_run.txt`); it reports `all cross-validation checks PASS`.

### From scratch

```bash
javac UnifiedTreeExperimentL.java BiDRerun.java
java UnifiedTreeExperimentL      # writes edbt_pertrial.csv (full run, 12 algorithms)
java BiDRerun                    # writes edbt_bidij_rerun.csv (BiDijkstra re-run)
python3 edbt_pipeline.py         # rebuild all tables/figures from the new CSVs
```

Grids are generated from the deterministic seed `100·depth + trial`, so both
programs reproduce byte-identical CSVs. `java UnifiedTreeExperimentL --fast`
skips the two exponential-tree algorithms (RBFS, IDA*) for a quick smoke
test — run without the flag to reproduce the paper's data. On machines with
a display, a results window opens after the run; the CSV and console log are
identical either way.

## Notes

- The BiDijkstra rows shipped in `edbt_pertrial.csv` are replaced by the
  dedicated re-run records (`edbt_bidij_rerun.csv`) during analysis; the
  pipeline prints a note whenever this happens.
- Figures are rendered in grayscale only, as in the paper.

## Citation

If you use this code or data, please cite the paper. A formal citation will
be added here upon publication.
