"""Every probe's self-check, run without an economicspace checkout.

Each probe checks itself before it measures: it reproduces a shipped or
upstream figure, or recovers a known answer. Those checks need no catalog, so
CI can run them on every push. A check that needs a local download skips
when the download is absent.
"""

import os
import unittest

import numpy as np

from tests._util import builder_rows, fixture, quiet, synthetic_neos, synthetic_toliou

import probe_covariance
import probe_lambert
import probe_nhats
import probe_perihelion
import probe_plane_change
import probe_prospector
import probe_spin
import probe_taxonomy
import probe_taxonomy_prior
import probe_recovery

TOLIOU_CACHE = probe_perihelion.DEFAULT_CACHE
HAVE_TOLIOU = os.path.exists(os.path.join(TOLIOU_CACHE, "Toliou+_2021_MNRAS.dat.gz"))


class ShippedModel(unittest.TestCase):
    """The probes that price delta-v must reproduce the docstring's five figures."""

    def test_nhats_reference(self):
        self.assertTrue(quiet(probe_nhats.check_reference))

    def test_plane_change_reference(self):
        self.assertTrue(quiet(probe_plane_change.check_reference))

    def test_lambert_reference(self):
        self.assertTrue(quiet(probe_lambert.check_reference))

    def test_lambert_solver(self):
        self.assertTrue(quiet(probe_lambert.check_solver))

    def test_vectorised_and_scalar_shipped_agree(self):
        r = probe_nhats.REFERENCE
        vec, _ = probe_nhats.dv_out_shipped(r.semi_major_axis_au.to_numpy(float),
                                            r.eccentricity.to_numpy(float),
                                            r.inclination_deg.to_numpy(float))
        scalar = [probe_lambert.shipped_dv_out(a, e, i) for a, e, i in
                  zip(r.semi_major_axis_au, r.eccentricity, r.inclination_deg)]
        np.testing.assert_allclose(vec, scalar, atol=1e-12)


class Prospector(unittest.TestCase):
    """Adapted prospector code must reproduce prospector's own output."""

    def test_intercept_and_coplanar(self):
        self.assertTrue(quiet(probe_prospector.check_upstream))

    def test_cohesion(self):
        self.assertTrue(quiet(probe_spin.check_upstream))

    def test_shape_integrals_sum_to_two(self):
        # Holsapple's A_x + A_y + A_z = 2; the sphere's are 2/3 each.
        self.assertAlmostEqual(3 * probe_spin.A_SPHERE, 2.0, places=15)

    @unittest.skipUnless(HAVE_TOLIOU, "Toliou table not downloaded")
    def test_perihelion_upstream(self):
        grid = probe_perihelion.Grid(probe_perihelion.load_table(TOLIOU_CACHE))
        self.assertTrue(quiet(probe_perihelion.check_upstream, grid))


class Perihelion(unittest.TestCase):
    """The fast grid lookup against the full scan, on a synthetic grid."""

    def test_lookup_matches_full_scan(self):
        grid = probe_perihelion.Grid(synthetic_toliou())
        self.assertTrue(quiet(probe_perihelion.check_lookup, grid, synthetic_neos(), 300))

    def test_lookup_uses_fallback_for_empty_cells(self):
        grid = probe_perihelion.Grid(synthetic_toliou(fill=0.2))
        n = synthetic_neos()
        _, n_fallback = grid.lookup(n.semi_major_axis_au, n.eccentricity,
                                    n.inclination_deg, n.absolute_magnitude_h)
        self.assertGreater(n_fallback, 0)

    def test_expected_q_min_of_a_step(self):
        cum = np.zeros(26)
        cum[7:] = 1.0                      # all mass in the 0.35-0.40 AU bin
        self.assertAlmostEqual(float(probe_perihelion.expected_q_min(cum)), 0.375, places=12)


class Taxonomy(unittest.TestCase):

    def test_mnlogit_recovers_known_coefficients(self):
        self.assertTrue(quiet(probe_taxonomy_prior.check_recovery))

    def test_builder_rules_match_the_release(self):
        self.assertTrue(quiet(probe_taxonomy_prior.check_builder, builder_rows()))

    def test_neo_assumed_albedo_flips_at_2_5_au(self):
        a = np.array([1.0, 2.49, 2.51, 3.5])
        cls = probe_taxonomy_prior.rule_class(
            probe_taxonomy_prior.assumed_albedo(a, np.ones(4, bool)))
        self.assertEqual(list(cls), [1, 1, 0, 0])          # S, S, C, C

    def test_entry_for_falls_back_by_root(self):
        tax = {"S": {"group": "S-complex"}, "Ch": {"group": "C-complex"},
               "C": {"group": "C-complex"}, "Unknown": {"group": None}}
        self.assertEqual(probe_taxonomy.entry_for(tax, "Sq")["group"], "S-complex")
        self.assertEqual(probe_taxonomy.entry_for(tax, "Chx")["group"], "C-complex")
        self.assertIsNone(probe_taxonomy.entry_for(tax, "Zz")["group"])


class Recovery(unittest.TestCase):

    def setUp(self):
        with open(fixture("neodys_2000SG344.ke1"), encoding="utf-8") as fh:
            self.orbit = probe_covariance.parse(fh.read())

    def test_earth_sits_where_the_almanac_puts_it(self):
        self.assertTrue(quiet(probe_recovery.check_earth))

    def test_linear_mean_anomaly_growth_matches_draws(self):
        rng = np.random.default_rng(0)
        draws = probe_recovery.draws_from(self.orbit, 20000, rng)
        for years in (1.0, 10.0):
            lin = probe_recovery.sigma_m_linear(self.orbit, years)
            mc = probe_recovery.sigma_m_mc(draws, years)
            self.assertLess(abs(mc - lin) / lin, 0.03, years)

    def test_phase_function_is_one_at_opposition(self):
        self.assertAlmostEqual(float(probe_recovery.phase_function(0.0)), 1.0, places=12)

    def test_assess_finds_a_window_for_a_bright_earth_crosser(self):
        # 2000 SG344 at H 24.8 with a near-Earth orbit; it comes within V 24
        # at its 2028-2030 approaches in the two-body model.
        r = probe_recovery.assess(self.orbit, 20.0, probe_recovery.mjd(
            probe_recovery.datetime.date(2026, 11, 1)), probe_recovery.mjd(probe_recovery.END),
            200, np.random.default_rng(1))
        self.assertTrue(r["window"])
        self.assertLess(r["search_deg"], 1.0)


class Covariance(unittest.TestCase):

    def setUp(self):
        with open(fixture("neodys_2000SG344.ke1"), encoding="utf-8") as fh:
            self.orbit = probe_covariance.parse(fh.read())

    def test_parse_reads_elements_and_covariance(self):
        self.assertIsNotNone(self.orbit)
        self.assertAlmostEqual(self.orbit["kep"][0], 0.977483, places=6)
        self.assertEqual(self.orbit["cov"].shape, (6, 6))
        np.testing.assert_allclose(self.orbit["cov"], self.orbit["cov"].T)

    def test_covariance_matches_rms_line(self):
        self.assertLess(probe_covariance.rms_gap(self.orbit), 1e-3)

    def test_packed_designation(self):
        self.assertEqual(probe_covariance.packed("2000 SG344"), "2000SG344")
        self.assertEqual(probe_covariance.packed(" 433 "), "433")


if __name__ == "__main__":
    unittest.main()
