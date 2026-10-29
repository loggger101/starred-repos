"""Every self-check must FAIL on a planted defect.

A check that cannot fail proves nothing. Each test here breaks one thing the
check is meant to guard (a recorded figure, a constant, an exponent, a parse
order) and asserts that the check notices. These are the defects planted by
hand when each check was written, kept so they are planted on every push.
"""

import os
import unittest
from unittest import mock

import numpy as np

from tests._util import builder_rows, fixture, quiet, synthetic_neos, synthetic_toliou

import orbital
import probe_covariance
import probe_lambert
import probe_nhats
import probe_perihelion
import probe_plane_change
import probe_prospector
import probe_spin
import probe_taxonomy_prior
import probe_recovery

HAVE_TOLIOU = os.path.exists(os.path.join(probe_perihelion.DEFAULT_CACHE,
                                          "Toliou+_2021_MNRAS.dat.gz"))


def shifted(frame, column, row, by):
    out = frame.copy()
    out.loc[row, column] += by
    return out


class ShippedModel(unittest.TestCase):

    def test_nhats_reference_catches_a_moved_figure(self):
        bad = shifted(probe_nhats.REFERENCE, "docstring_dv_out", 3, 0.02)
        with mock.patch.object(probe_nhats, "REFERENCE", bad):
            self.assertFalse(quiet(probe_nhats.check_reference))

    def test_plane_change_reference_catches_a_moved_figure(self):
        bad = shifted(probe_plane_change.REFERENCE, "docstring_dv_out", 0, 0.02)
        with mock.patch.object(probe_plane_change, "REFERENCE", bad):
            self.assertFalse(quiet(probe_plane_change.check_reference))

    def test_nhats_reference_catches_a_moved_constant(self):
        with mock.patch.object(probe_nhats, "V_EARTH_KM_S", 29.9):
            self.assertFalse(quiet(probe_nhats.check_reference))

    def test_lambert_reference_catches_a_moved_figure(self):
        bad = shifted(probe_nhats.REFERENCE, "docstring_dv_out", 2, 0.02)
        with mock.patch.object(probe_nhats, "REFERENCE", bad):
            self.assertFalse(quiet(probe_lambert.check_reference))

    def test_lambert_solver_catches_a_wrong_velocity(self):
        real = orbital.lambert_izzo

        def off(r1, r2, tof, **kw):
            v1, v2 = real(r1, r2, tof, **kw)
            return v1 * (1 + 1e-6), v2
        with mock.patch.object(probe_lambert.ob, "lambert_izzo", off):
            self.assertFalse(quiet(probe_lambert.check_solver))


class Prospector(unittest.TestCase):

    def test_intercept_catches_a_moved_upstream_figure(self):
        bad = probe_prospector.UPSTREAM_INTERCEPT + np.array([0, 0, 1e-4, 0, 0])
        with mock.patch.object(probe_prospector, "UPSTREAM_INTERCEPT", bad):
            self.assertFalse(quiet(probe_prospector.check_upstream))

    def test_coplanar_check_catches_a_second_earth_speed(self):
        with mock.patch.object(probe_prospector, "V_EARTH_KM_S", 29.785):
            self.assertFalse(quiet(probe_prospector.check_upstream))

    def test_cohesion_catches_upstream_without_its_bug(self):
        with mock.patch.object(probe_spin, "A_X_UPSTREAM", 2.0 / 3.0):
            self.assertFalse(quiet(probe_spin.check_upstream))

    def test_cohesion_catches_a_wrong_closed_form(self):
        real = probe_spin.cohesionless_period_h

        def bad(rho, a_x=probe_spin.A_SPHERE):
            if a_x == probe_spin.A_SPHERE:
                rho = np.asarray(rho, float)
                s = probe_spin.friction_coefficient()
                omega2 = 4 * np.pi * s * rho * probe_spin.G / (1 / np.sqrt(3) + s)
                return 2 * np.pi / np.sqrt(omega2) / probe_spin.SECONDS_PER_HOUR
            return real(rho, a_x)
        with mock.patch.object(probe_spin, "cohesionless_period_h", bad):
            self.assertFalse(quiet(probe_spin.check_upstream))

    @unittest.skipUnless(HAVE_TOLIOU, "Toliou table not downloaded")
    def test_perihelion_catches_a_moved_upstream_figure(self):
        grid = probe_perihelion.Grid(probe_perihelion.load_table(probe_perihelion.DEFAULT_CACHE))
        bad = shifted(probe_perihelion.UPSTREAM, "q_min_upstream", 2, 1e-5)
        with mock.patch.object(probe_perihelion, "UPSTREAM", bad):
            self.assertFalse(quiet(probe_perihelion.check_upstream, grid))


class Perihelion(unittest.TestCase):

    def test_lookup_catches_a_shifted_axis(self):
        grid = probe_perihelion.Grid(synthetic_toliou())
        grid.axes[1] = grid.axes[1] + 0.03
        self.assertFalse(quiet(probe_perihelion.check_lookup, grid, synthetic_neos(), 300))


class Taxonomy(unittest.TestCase):

    def test_recovery_catches_rotated_labels(self):
        real = probe_taxonomy_prior.fit_mnlogit
        with mock.patch.object(probe_taxonomy_prior, "fit_mnlogit",
                               lambda X, y, k, **kw: real(X, (y + 1) % k, k, **kw)):
            self.assertFalse(quiet(probe_taxonomy_prior.check_recovery))

    def test_builder_check_catches_a_moved_split(self):
        with mock.patch.object(probe_taxonomy_prior, "C_BELOW_ALBEDO", 0.11):
            self.assertFalse(quiet(probe_taxonomy_prior.check_builder, builder_rows()))

    def test_builder_check_catches_a_moved_neo_bin(self):
        bad = [(0.0, 1.3, 0.187), (1.3, 2.0, 0.17), (2.0, 2.5, 0.09),
               (2.5, 2.82, 0.062), (2.82, 1e9, 0.037)]
        with mock.patch.object(probe_taxonomy_prior, "ALBEDO_BY_A_NEO", bad):
            self.assertFalse(quiet(probe_taxonomy_prior.check_builder, builder_rows()))


class Recovery(unittest.TestCase):

    def test_earth_check_catches_a_shifted_mean_anomaly(self):
        bad = dict(probe_recovery.EARTH, M0=probe_recovery.EARTH["M0"] + np.radians(0.1))
        with mock.patch.object(probe_recovery, "EARTH", bad):
            self.assertFalse(quiet(probe_recovery.check_earth))

    def test_linear_growth_without_the_a_m_correlation_disagrees(self):
        with open(fixture("neodys_2000SG344.ke1"), encoding="utf-8") as fh:
            orbit = probe_covariance.parse(fh.read())
        draws = probe_recovery.draws_from(orbit, 20000, np.random.default_rng(0))
        mc = probe_recovery.sigma_m_mc(draws, 10.0)
        no_cross = dict(orbit, cov=orbit["cov"].copy())
        no_cross["cov"][0, 5] = no_cross["cov"][5, 0] = 0.0
        lin = probe_recovery.sigma_m_linear(no_cross, 10.0)
        self.assertGreater(abs(mc - lin) / lin, 0.05)


class Covariance(unittest.TestCase):

    def test_rms_check_catches_column_major_fill(self):
        with open(fixture("neodys_2000SG344.ke1"), encoding="utf-8") as fh:
            orbit = probe_covariance.parse(fh.read())
        vals = orbit["cov"][np.triu_indices(6)]
        wrong = np.zeros((6, 6))
        wrong[np.tril_indices(6)] = vals
        orbit["cov"] = wrong + wrong.T - np.diag(np.diag(wrong))
        with np.errstate(invalid="ignore"):
            gap = probe_covariance.rms_gap(orbit)
        self.assertFalse(gap < 1e-3)


if __name__ == "__main__":
    unittest.main()
