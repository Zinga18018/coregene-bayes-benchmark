# Model and evaluation

## Scope

This is a small Bayesian model for synthetic evidence integration. The graph is stipulated; all hyperparameters are fixed before simulation. Core membership is a simulated binary label. There is no real disease endpoint, tissue, gene symbol, variant coordinate, participant or clinical actionability outcome.

The implementation is bounded to four core indicators, four peripheral activities, eight variant evidence values and three uncertain network links. A score is called a posterior only because it follows from the probability model below. It is not a posterior probability of a real gene being causal or therapeutically actionable.

## Generative model

Let `C_k ∈ {0,1}` denote latent synthetic core membership for `k=0,…,3`.

`C_k ~ Bernoulli(π)`, independently, with `π=0.25`.

For network links `j=0,1,2`, independently:

`W_j ~ Beta(a_j, b_j)` with shape pairs `(2,3), (3,2), (2,2)`.

Define peripheral activity:

`A_j = (1-W_j) C_j + W_j C_(j+1)` for `j=0,1,2`, and `A_3=C_3`.

Each variant index `v` has two candidate peripheral genes `(g_v,h_v)` and an independent continuous mixture weight `Q_v ~ Beta(α_v,β_v)`. The eight fixed gene pairs and Beta shapes are stored in `Config` and copied into every run's metadata. These are synthetic mixture weights, not estimated fine-mapping probabilities or data-derived annotations.

Observed variant evidence:

`Z_v | C,W,Q_v ~ Normal(δ [Q_v A_(g_v) + (1-Q_v) A_(h_v)], σ_z²)`.

Observed direct core evidence:

`Y_k | C_k ~ Normal(γ C_k, σ_y²)`.

Defaults are `δ=2.8`, `σ_z=0.8`, `γ=1.2`, `σ_y=1.0`. The direct evidence is an illustrative second likelihood, not a real rare-variant burden or drug-target assay. Its conditional independence from variant evidence is an explicit simulation assumption.

After generation, each observed value is independently removed with probability `0.15` (missing completely at random). Inference receives only `Z` and `Y`; saved truth and sampled weights are used for evaluation and audit, never as inference inputs. Missing values have likelihood contribution one. Missingness that depends on unobserved values is not modelled.

## Posterior computation

The 16 binary vectors `c` are enumerated. For each vector:

`p(c | Z,Y) ∝ p(c) × ∫ p(w) [∏_v ∫ p(q_v) p(Z_v | c,w,q_v) dq_v] dw × ∏_k p(Y_k | c_k)`.

Products include only observed values. **The shared network weights are integrated after multiplying variant likelihoods.** Integrating each variant's network uncertainty separately would incorrectly remove the dependence induced by the shared `W`.

Gauss–Legendre nodes on `[0,1]` are reweighted by the Beta density and normalised. The default uses eight nodes per network dimension (`8³=512` joint settings) and sixteen nodes per mapping dimension. The three network dimensions use a product rule; independent variant mappings are integrated one at a time inside the conditional likelihood. All likelihood accumulation and normalisation use log space.

Binary-state enumeration is exact; continuous integration is approximate. Beta shapes below one are rejected because endpoint singularities require a different quadrature treatment. For the first ten matched datasets, the run compares orders `(8,16)` with `(12,24)`. This diagnostic estimates local numerical stability; it is not a global error bound. Runtime grows exponentially with core indicators and network dimensions, so the implementation is not genome scale.

Marginal core probability is `P(C_k=1 | Z,Y) = Σ_c c_k p(c | Z,Y)`.

## Comparisons

| Model | Mapping uncertainty | Network uncertainty | Network topology |
|---|---|---|---|
| `full` | Integrated | Integrated | Stipulated chain |
| `fixed_mapping` | Replaced by Beta mean | Integrated | Stipulated chain |
| `fixed_network` | Integrated | Replaced by Beta mean | Stipulated chain |
| `fixed_both` | Replaced by Beta mean | Replaced by Beta mean | Stipulated chain |
| `no_network` | Integrated | Removed (`W=0`) | Identity |

`fixed_both` is the plug-in baseline. `no_network` is a deliberately misspecified structural baseline, not an equivalent model without an uncertainty component. All models use the same prior, noise assumptions, direct evidence, and identical observed datasets. The majority non-core accuracy is included as a reference because only 25% of simulated memberships are core on average.

## Evaluation

**Matched scenario:** 300 independent draws from the specified generative model.

**Noise stress scenario:** 150 new independent draws with both observation-noise standard deviations multiplied by `1.8`, while inference retains the original standard deviations. No attempt is made to correct or hide the resulting coverage loss.

There is no training, parameter fitting, threshold selection or hyperparameter search. Accordingly a train/test split would not serve the usual model-selection purpose here. Evaluation is prior predictive, not fixed-parameter frequentist calibration and not validation on unseen biological cohorts.

Metrics:

- Brier score: mean `(p_k-C_k)²`. Reported Monte Carlo standard error is calculated from independent replicate-level mean scores, respecting within-replicate gene dependence.
- Binary log loss: mean `-C_k log(p_k)-(1-C_k) log(1-p_k)`. Numerical probabilities are clipped to `[10^-12,1-10^-12]` only for this reporting calculation.
- Accuracy at the fixed threshold `0.5`, accompanied by observed prevalence and majority non-core accuracy.
- Ten descriptive probability bins with counts, mean predicted probability and observed core fraction. Empty bins stay empty. No independent-gene binomial confidence intervals are asserted.
- Joint 90% credible set: sort the 16 states by posterior probability and include states until accumulated mass reaches 0.9. Ties use stable enumeration order. Empirical coverage is the fraction of independent datasets whose true state lies in that set. A Wilson 95% interval reflects finite replicate count. Mean actual set mass and size are reported because discreteness makes the set mass exceed 0.9.

Small differences between methods are not presented as significance claims. The stress scenario demonstrates sensitivity to a particular misspecification; it does not exhaust possible failure modes.

## Independent checks

The tests compare Beta quadrature against analytic first and second moments, recover the prior when all observations are missing, recover the analytic independent Bernoulli posterior when only direct Gaussian evidence is available, and compare fixed-weight inference against a separately written scalar enumeration. Additional checks cover zero signal, fractional propagation from integer truth, reproducibility, invalid observations, quadrature convergence and credible-set mass.

## Boundaries before real genetics

This model omits LD, allelic direction/harmonisation, ancestry structure, sample overlap, horizontal pleiotropy, tissue-specific effects, batch effects, shared annotation errors, multimodal evidence dependence, feedback/cyclic regulation, and causal identifiability. It uses an acyclic positive convex mixture with fixed effect scales. There is no modelling of inhibition or protective versus risk effects. The Beta distributions are stipulated; their concentration does not come from real evidence quality.

It does not estimate `P(therapeutically actionable | disease)`. A high synthetic core probability cannot establish causality, intervention efficacy, druggability, direction of treatment effect or safety. Moving toward a real target-prioritisation method requires additional evidence models, external data validation, and expert review.

The conceptual alignment is with uncertainty propagation in the [official UNIL project description](https://wp.unil.ch/sgg/2027-fbm-phd-fellowship-applications/). Numerical quadrature uses the [NumPy Gauss–Legendre implementation](https://numpy.org/doc/2.3/reference/generated/numpy.polynomial.legendre.leggauss.html). No external dataset was downloaded or analysed.
