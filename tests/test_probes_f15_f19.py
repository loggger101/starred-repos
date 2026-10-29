"""Checks behind F19 (the oracle's two-body error, via Horizons), and their planted defects.

F15's recovery checks live in test_checks.py and test_planted.py.
"""

import unittest
from unittest import mock

import numpy as np

from tests._util import fixture

import orbital as ob
import probe_horizons as ph


def kepler_table(a_au=1.2, e=0.6, days=400, step=2.0):
    """A Horizons-shaped vector table (JD, r km, v km/s) for a two-body orbit."""
    t = np.arange(0.0, days, step)
    r, v = ob.propagate_elements(a_au * ob.AU_KM, e, 0.2, 0.3, 0.4, 0.0, t * ob.DAY_S)
    return np.column_stack([ph.JD_START + t, r.T, v.T])


class Horizons(unittest.TestCase):

    def setUp(self):
        with open(fixture("horizons_vectors.txt"), encoding="utf-8") as fh:
            self.table = ph.parse_horizons(fh.read())

    def test_parse_reads_the_vector_block(self):
        self.assertEqual(self.table.shape, (30, 7))
        self.assertAlmostEqual(self.table[0, 0], ph.JD_START, places=6)
        # Horizons' own first row for 2000 SG344 on 2027-01-01, as served 2026-10-29.
        self.assertEqual(self.table[0, 1], 4.056028413546421e7)
        self.assertEqual(self.table[0, 2], 1.318753402696495e8)

    def test_parse_refuses_an_answer_without_data(self):
        with self.assertRaises(LookupError):
            ph.parse_horizons("No matches found.")

    def test_hermite_track_returns_its_own_nodes(self):
        track = ph.Track(self.table)
        r, v = track.state(track.t[5:8])
        np.testing.assert_allclose(r, self.table[5:8, 1:4], rtol=0, atol=1e-6)
        np.testing.assert_allclose(v, self.table[5:8, 4:7], rtol=0, atol=1e-9)

    def test_interpolation_error_falls_as_the_fourth_power(self):
        err, ratio = ph.interpolation_error_km(kepler_table())
        self.assertTrue(10 < ratio < 20, ratio)
        self.assertLess(err / (0.25 * 365.25 * ob.DAY_S), 1e-3)     # the probe's own bar

    def test_planted_a_linear_interpolant_is_caught(self):
        def linear_state(self, t_s):
            t_s = np.atleast_1d(np.asarray(t_s, float))
            r = np.column_stack([np.interp(t_s, self.t, self.r[:, k]) for k in range(3)])
            return r, np.zeros_like(r)
        with mock.patch.object(ph.Track, "state", linear_state):
            err, ratio = ph.interpolation_error_km(kepler_table())
        # A linear interpolant's error falls as the square, not the fourth power,
        # and is far larger; either breaks the check's premise.
        self.assertTrue(ratio < 6 or err > 1000, (ratio, err))

    def test_elements_row_maps_to_radians(self):
        row = [ph.JD_START, 0.0669, 0.9122, 0.1132, 191.69, 275.61, 0.0, 1.0197, 329.75,
               325.6, 0.97756, 1.0429, 353.03]
        a, e, i, om, w, ma = ph.elements_from(row)
        self.assertAlmostEqual(a / ob.AU_KM, 0.97756)
        self.assertAlmostEqual(ma, np.radians(329.75))
        self.assertAlmostEqual(e, 0.0669)

    def test_command_for_numbers_and_designations(self):
        self.assertEqual(ph.command_for("433"), "433;")
        self.assertEqual(ph.command_for("2000 SG344"), "DES=2000 SG344;")


if __name__ == "__main__":
    unittest.main()
