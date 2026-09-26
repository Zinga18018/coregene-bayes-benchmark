"""Exact binary-state enumeration, with quadrature for continuous latent weights.

This is a toy statistical model, not a genomic association or causal inference tool.
See METHOD.md for its generative assumptions and limitations.
"""

from dataclasses import asdict, dataclass
from itertools import product
import math

import numpy as np


def logsumexp(values, axis=None):
    """Numerically stable log(sum(exp(values))) for finite log densities."""
    values = np.asarray(values, dtype=float)
    maximum = np.max(values, axis=axis, keepdims=True)
    result = maximum + np.log(np.sum(np.exp(values - maximum), axis=axis, keepdims=True))
    return np.squeeze(result, axis=axis) if axis is not None else float(result.item())


@dataclass(frozen=True)
class Config:
    """Fixed synthetic design. No parameters are fitted to benchmark labels."""

    core_probability: float = 0.25
    signal_scale: float = 2.8
    variant_noise: float = 0.8
    direct_scale: float = 1.2
    direct_noise: float = 1.0
    mapping_pairs: tuple = ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 3), (2, 0), (3, 1))
    mapping_beta: tuple = ((2., 2.), (4., 2.), (2., 4.), (5., 2.), (2., 2.), (2., 3.), (3., 2.), (4., 3.))
    network_beta: tuple = ((2., 3.), (3., 2.), (2., 2.))

    def __post_init__(self):
        if not np.isfinite([self.core_probability, self.signal_scale, self.variant_noise,
                            self.direct_scale, self.direct_noise]).all():
            raise ValueError("All scalar parameters must be finite")
        if not 0 < self.core_probability < 1:
            raise ValueError("core_probability must be strictly between zero and one")
        if self.variant_noise <= 0 or self.direct_noise <= 0:
            raise ValueError("Noise standard deviations must be positive")
        if len(self.mapping_pairs) != len(self.mapping_beta) or len(self.network_beta) != 3:
            raise ValueError("This bounded prototype uses four genes and three network edges")
        for pair in self.mapping_pairs:
            if len(pair) != 2 or any(g not in range(4) for g in pair):
                raise ValueError("Mapping endpoints must be gene indices 0..3")
        for a, b in self.mapping_beta + self.network_beta:
            if not np.isfinite([a, b]).all() or a < 1 or b < 1:
                raise ValueError("This quadrature implementation supports Beta shapes >= 1")

    def to_dict(self):
        return asdict(self)


def beta_quadrature(a, b, order):
    """Gauss-Legendre nodes reweighted by a Beta density on [0, 1].

    Supported Beta shapes are >= 1. Weights are normalised; general smooth
    expectations are approximate, so benchmark output includes convergence checks.
    """
    if a < 1 or b < 1 or order < 2:
        raise ValueError("Require Beta shapes >= 1 and quadrature order >= 2")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    log_pdf = ((a - 1) * np.log(nodes) + (b - 1) * np.log1p(-nodes)
               + math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b))
    weights = weights * np.exp(log_pdf)
    return nodes, weights / weights.sum()


def states_and_prior(config):
    states = np.asarray(list(product((0., 1.), repeat=4)))
    count = states.sum(axis=1)
    log_prior = (count * np.log(config.core_probability)
                 + (4 - count) * np.log1p(-config.core_probability))
    return states, log_prior


def propagated_activity(states, network_weights):
    """Network rows 0..2 mix own and next core; row 3 retains its own core."""
    states = np.asarray(states, dtype=float)
    weights = np.asarray(network_weights, dtype=float)
    activity = np.broadcast_to(states, (len(weights),) + states.shape).copy()
    activity[:, :, :3] = ((1 - weights[:, None, :]) * states[None, :, :3]
                          + weights[:, None, :] * states[None, :, 1:])
    return activity


def simulate(config, rng, noise_multiplier=1.0, missing_rate=0.15):
    """Draw independent latent weights and truth, then noisy observed evidence.

    Missingness is independent of every value and label (MCAR). A noise multiplier
    other than one intentionally violates the inference model for stress testing.
    """
    if not 0 <= missing_rate <= 1 or noise_multiplier <= 0:
        raise ValueError("Invalid missingness or noise multiplier")
    truth = rng.binomial(1, config.core_probability, size=4)
    network = np.array([rng.beta(a, b) for a, b in config.network_beta])
    mapping = np.array([rng.beta(a, b) for a, b in config.mapping_beta])
    activity = propagated_activity(truth[None, :], network[None, :])[0, 0]
    means = np.array([config.signal_scale * (q * activity[g] + (1 - q) * activity[h])
                      for q, (g, h) in zip(mapping, config.mapping_pairs)])
    variant = means + rng.normal(0, config.variant_noise * noise_multiplier, len(means))
    direct = config.direct_scale * truth + rng.normal(0, config.direct_noise * noise_multiplier, 4)
    variant[rng.random(len(variant)) < missing_rate] = np.nan
    direct[rng.random(4) < missing_rate] = np.nan
    return {"truth": truth, "variant": variant, "direct": direct,
            "network_weights": network, "mapping_weights": mapping}


class PreparedModel:
    """Precompute latent means; infer a posterior over all 16 binary core states."""

    MODES = ("full", "fixed_mapping", "fixed_network", "fixed_both", "no_network")

    def __init__(self, config=None, mode="full", network_order=8, mapping_order=16):
        self.config = config or Config()
        if mode not in self.MODES:
            raise ValueError(f"Unknown mode: {mode}")
        if network_order < 2 or mapping_order < 2:
            raise ValueError("Quadrature orders must be at least two")
        self.mode = mode
        self.network_order, self.mapping_order = network_order, mapping_order
        self.states, self.log_prior = states_and_prior(self.config)
        if mode in ("fixed_network", "fixed_both"):
            network = np.array([[a / (a + b) for a, b in self.config.network_beta]])
            network_probability = np.ones(1)
        elif mode == "no_network":
            network, network_probability = np.zeros((1, 3)), np.ones(1)
        else:
            grids = [beta_quadrature(a, b, network_order) for a, b in self.config.network_beta]
            indices = np.array(list(product(range(network_order), repeat=3)))
            network = np.column_stack([grids[j][0][indices[:, j]] for j in range(3)])
            network_probability = np.prod(np.column_stack([grids[j][1][indices[:, j]]
                                                           for j in range(3)]), axis=1)
        self.log_network_probability = np.log(network_probability)
        activity = propagated_activity(self.states, network)
        self.variant_means, self.log_mapping_probability = [], []
        for (g, h), (a, b) in zip(self.config.mapping_pairs, self.config.mapping_beta):
            if mode in ("fixed_mapping", "fixed_both"):
                q, probability = np.array([a / (a + b)]), np.ones(1)
            else:
                q, probability = beta_quadrature(a, b, mapping_order)
            means = self.config.signal_scale * (q * activity[:, :, g, None]
                                                + (1 - q) * activity[:, :, h, None])
            self.variant_means.append(means)
            self.log_mapping_probability.append(np.log(probability))

    def infer(self, variant, direct):
        variant, direct = np.asarray(variant, float), np.asarray(direct, float)
        if variant.shape != (len(self.config.mapping_pairs),) or direct.shape != (4,):
            raise ValueError("Expected eight variant observations and four direct observations")
        if np.isinf(variant).any() or np.isinf(direct).any():
            raise ValueError("Use NaN for missing evidence; infinity is invalid")
        c = self.config
        log_likelihood = np.zeros((len(self.log_network_probability), len(self.states)))
        for value, means, log_weights in zip(variant, self.variant_means, self.log_mapping_probability):
            if np.isnan(value):
                continue  # Marginalise an unobserved datum: its likelihood contribution is one.
            log_density = (-0.5 * ((value - means) / c.variant_noise) ** 2
                           - np.log(c.variant_noise * np.sqrt(2 * np.pi)))
            log_likelihood += logsumexp(log_density + log_weights, axis=2)
        integrated = logsumexp(log_likelihood + self.log_network_probability[:, None], axis=0)
        observed = ~np.isnan(direct)
        if observed.any():
            integrated += np.sum(-0.5 * ((direct[observed] - c.direct_scale * self.states[:, observed])
                                         / c.direct_noise) ** 2
                                 - np.log(c.direct_noise * np.sqrt(2 * np.pi)), axis=1)
        log_joint = self.log_prior + integrated
        probability = np.exp(log_joint - logsumexp(log_joint))
        marginal = probability @ self.states
        return {"state_probability": probability, "core_probability": marginal}

    def credible_set(self, state_probability, mass=0.9):
        if not 0 < mass <= 1:
            raise ValueError("Credible-set mass must be in (0, 1]")
        order = np.argsort(-np.asarray(state_probability), kind="stable")
        cumulative = np.cumsum(np.asarray(state_probability)[order])
        size = min(int(np.searchsorted(cumulative, mass)) + 1, len(order))
        chosen = order[:size]
        return chosen, float(np.asarray(state_probability)[chosen].sum())
