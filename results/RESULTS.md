# Executed synthetic benchmark

Seed: 20260925. Created: 2026-09-25T13:59:21.303369+00:00.

All numbers below concern generated toy data. They are not evidence of biological validity.

| Scenario | Model | Brier (lower better) | Log loss | Joint 90% set coverage | Mean set mass |
|---|---|---:|---:|---:|---:|
| matched | full | 0.0754 | 0.2437 | 0.943 | 0.942 |
| matched | fixed_mapping | 0.0760 | 0.2468 | 0.940 | 0.944 |
| matched | fixed_network | 0.0769 | 0.2488 | 0.940 | 0.947 |
| matched | fixed_both | 0.0775 | 0.2531 | 0.923 | 0.952 |
| matched | no_network | 0.1155 | 0.4348 | 0.783 | 0.954 |
| noise_misspecification | full | 0.1538 | 0.5430 | 0.707 | 0.946 |
| noise_misspecification | fixed_mapping | 0.1533 | 0.5597 | 0.700 | 0.949 |
| noise_misspecification | fixed_network | 0.1572 | 0.5967 | 0.667 | 0.956 |
| noise_misspecification | fixed_both | 0.1553 | 0.6166 | 0.653 | 0.958 |
| noise_misspecification | no_network | 0.1778 | 0.7576 | 0.600 | 0.958 |

Coverage uses a joint credible set over the 16 possible core-membership states, not a confidence interval for a gene.
Discrete sets often contain more than 90% posterior mass; compare coverage with mean set mass as well as 0.90.
The stress scenario increases observation noise by 1.8 without informing inference.
Calibration bins are descriptive: gene outcomes within each replicate are dependent.
Monte Carlo standard errors and Wilson coverage intervals are recorded in summary.json.

Largest marginal change in the ten-dataset quadrature check: 6.37e-07.
This is a local convergence diagnostic, not a proof for every possible input.

There is no train/test split: no model parameter was fitted or tuned. Independent datasets are prior-predictive draws.
