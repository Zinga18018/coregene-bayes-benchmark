"""Independent limits and failure cases, rather than duplicated implementation."""

import unittest

import numpy as np

from coregene_bayes.model import Config, PreparedModel, beta_quadrature, propagated_activity, simulate


class BayesianModelTests(unittest.TestCase):
    def test_integer_truth_preserves_fractional_network_activity(self):
        activity = propagated_activity(np.array([[0, 1, 0, 1]]), np.array([[.25, .5, .75]]))
        np.testing.assert_allclose(activity[0, 0], [.25, .5, .75, 1.])

    def test_beta_quadrature_matches_analytic_moments(self):
        a, b = 2., 3.
        q, w = beta_quadrature(a, b, 12)
        self.assertAlmostEqual(float(w @ q), a/(a+b), places=12)
        self.assertAlmostEqual(float(w @ (q*q)), a*(a+1)/((a+b)*(a+b+1)), places=12)

    def test_no_evidence_recovers_prior_in_every_mode(self):
        for mode in PreparedModel.MODES:
            model = PreparedModel(mode=mode)
            posterior = model.infer(np.full(8, np.nan), np.full(4, np.nan))
            np.testing.assert_allclose(posterior["core_probability"], .25, atol=1e-12)
            self.assertAlmostEqual(float(posterior["state_probability"].sum()), 1.)

    def test_direct_only_has_analytic_independent_bernoulli_posterior(self):
        model = PreparedModel()
        y = np.array([-1., .3, 1.2, 2.8])
        c = model.config
        log_odds = (np.log(c.core_probability/(1-c.core_probability))
                    + (2*y*c.direct_scale-c.direct_scale**2)/(2*c.direct_noise**2))
        analytic = 1/(1+np.exp(-log_odds))
        posterior = model.infer(np.full(8, np.nan), y)
        np.testing.assert_allclose(posterior["core_probability"], analytic, atol=1e-12)

    def test_zero_signal_is_uninformative(self):
        model = PreparedModel(Config(signal_scale=0., direct_scale=0.))
        posterior = model.infer(np.arange(8)-4., np.arange(4))
        np.testing.assert_allclose(posterior["core_probability"], .25, atol=1e-12)

    def test_fixed_weights_agree_with_independent_enumeration(self):
        model = PreparedModel(mode="fixed_both")
        c = model.config
        z, y = np.linspace(-1, 2, 8), np.array([.1, 1.4, -.4, 1.8])
        masses = []
        for state in model.states:
            activity = state.copy()
            for i, (a, b) in enumerate(c.network_beta):
                activity[i] = b/(a+b)*state[i] + a/(a+b)*state[i+1]
            likelihood = 1.
            for value, (g, h), (a, b) in zip(z, c.mapping_pairs, c.mapping_beta):
                mean = c.signal_scale*(a/(a+b)*activity[g] + b/(a+b)*activity[h])
                likelihood *= np.exp(-.5*((value-mean)/c.variant_noise)**2)
            likelihood *= np.exp(-.5*np.sum(((y-c.direct_scale*state)/c.direct_noise)**2))
            masses.append(likelihood*c.core_probability**sum(state)*(1-c.core_probability)**(4-sum(state)))
        expected = np.array(masses)/sum(masses)
        np.testing.assert_allclose(model.infer(z, y)["state_probability"], expected, atol=1e-12)

    def test_seed_reproduces_synthetic_data(self):
        a = simulate(Config(), np.random.default_rng(2))
        b = simulate(Config(), np.random.default_rng(2))
        for key in a:
            np.testing.assert_equal(a[key], b[key])

    def test_uncertainty_changes_posterior(self):
        sample = simulate(Config(), np.random.default_rng(17), missing_rate=0)
        a = PreparedModel(mode="full").infer(sample["variant"], sample["direct"])
        b = PreparedModel(mode="fixed_both").infer(sample["variant"], sample["direct"])
        self.assertGreater(np.max(np.abs(a["core_probability"]-b["core_probability"])), .001)

    def test_more_quadrature_nodes_agree_for_nontrivial_example(self):
        sample = simulate(Config(), np.random.default_rng(71), missing_rate=0)
        a = PreparedModel(network_order=8, mapping_order=16).infer(sample["variant"], sample["direct"])
        b = PreparedModel(network_order=12, mapping_order=24).infer(sample["variant"], sample["direct"])
        np.testing.assert_allclose(a["core_probability"], b["core_probability"], atol=1e-4)

    def test_rejects_invalid_shapes_and_infinity(self):
        model = PreparedModel()
        with self.assertRaises(ValueError):
            model.infer([1], [1])
        with self.assertRaises(ValueError):
            model.infer(np.full(8, np.inf), np.zeros(4))
        with self.assertRaises(ValueError):
            Config(variant_noise=0)
        with self.assertRaises(ValueError):
            beta_quadrature(.5, 1, 8)

    def test_joint_credible_set_attains_requested_mass(self):
        model = PreparedModel()
        p = model.infer(np.zeros(8), np.zeros(4))["state_probability"]
        selected, mass = model.credible_set(p, .9)
        self.assertGreaterEqual(mass, .9)
        if len(selected) > 1:
            self.assertLess(float(p[selected[:-1]].sum()), .9)


if __name__ == "__main__":
    unittest.main()
