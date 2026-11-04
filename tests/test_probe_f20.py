"""Checks behind F20 (JPL's outbound legs), and their planted defects."""

import json
import unittest
import urllib.error
from unittest import mock

import numpy as np
import pandas as pd

from tests._util import fixture, quiet

import probe_nhats_legs as pl
import snapshots


def recorded():
    """Eight U <= 2 trajectories: JPL's answers, and the catalog row of each."""
    with open(fixture("nhats_legs.json"), encoding="utf-8") as fh:
        entries = json.load(fh)
    rows = [dict(pl.parse_object(json.dumps(e["answer"])), **e["elements"]) for e in entries]
    return pd.DataFrame(rows)


class Parse(unittest.TestCase):

    def test_reads_every_leg_of_the_cheapest_trajectory(self):
        with open(fixture("nhats_legs.json"), encoding="utf-8") as fh:
            answer = json.load(fh)[0]["answer"]
        legs = pl.parse_object(json.dumps(answer))
        self.assertEqual(legs["dv_total"], float(answer["min_dv_traj"]["dv_total"]))
        self.assertEqual(legs["dur_out"], answer["min_dv_traj"]["dur_out"])
        self.assertEqual(legs["occ"], int(answer["occ"]))

    def test_an_answer_with_no_trajectory_is_none(self):
        self.assertIsNone(pl.parse_object('{"des": "x"}'))

    def test_julian_date(self):
        self.assertEqual(pl.julian_date("2000-01-01"), 2451544.5)
        self.assertEqual(pl.julian_date("2020-01-01"), pl.LAUNCH_JD[0])
        self.assertEqual(pl.julian_date("2046-01-01"), pl.LAUNCH_JD[1])


class Legs(unittest.TestCase):

    def setUp(self):
        self.m = recorded()

    def test_legs_sum_and_fit_the_parking_orbit(self):
        self.assertTrue(quiet(pl.check_legs, self.m))

    def test_planted_a_dropped_leg_is_caught(self):
        self.assertFalse(quiet(pl.check_legs, self.m.assign(vrel_dep_neo=0.0)))

    def test_planted_the_wrong_parking_orbit_is_caught(self):
        with mock.patch.object(pl, "JPL_PARK_ALT_KM", 200.0):
            self.assertFalse(quiet(pl.check_legs, self.m))

    def test_200_km_departure_costs_more_than_400_km(self):
        self.assertGreater(pl.parking_dv(4.0, 200.0), pl.parking_dv(4.0, 400.0))


class Window(unittest.TestCase):

    def setUp(self):
        self.m = recorded()

    def test_jpl_legs_lie_inside_the_search(self):
        self.assertTrue(quiet(pl.check_window, self.m))

    def test_the_latest_launch_recorded_is_inside(self):
        self.assertTrue(quiet(pl.check_window, self.m.assign(launch_jd=2468341.5)))

    def test_planted_the_first_window_is_caught(self):
        # It stopped at 2045-01-01, and 11 of 415 launches fell after that.
        with mock.patch.object(pl, "LAUNCH_JD", (2458849.5, 2467981.5)):
            self.assertFalse(quiet(pl.check_window, self.m.assign(launch_jd=2468341.5)))

    def test_planted_a_leg_longer_than_the_search_flies_is_caught(self):
        self.assertFalse(quiet(pl.check_window, self.m.assign(dur_out=pl.MAX_TOF_D + 1)))


class Replay(unittest.TestCase):

    def setUp(self):
        self.m = recorded()

    def test_the_solver_flies_jpl_trajectories(self):
        self.assertTrue(quiet(pl.check_replay, self.m))

    def test_planted_a_shifted_epoch_is_caught(self):
        self.assertFalse(quiet(pl.check_replay, self.m.assign(
            element_epoch_jd=self.m.element_epoch_jd + 30.0)))

    def test_planted_a_moved_earth_is_caught(self):
        with mock.patch.dict(pl.EARTH, M0=pl.EARTH["M0"] + np.radians(2.0)):
            self.assertFalse(quiet(pl.check_replay, self.m))

    def test_planted_swapped_legs_are_caught(self):
        swapped = self.m.assign(v_dep_earth=self.m.vrel_arr_neo, vrel_arr_neo=self.m.v_dep_earth)
        self.assertFalse(quiet(pl.check_replay, swapped))

    def test_no_well_known_orbit_is_a_failure_not_a_pass(self):
        self.assertFalse(quiet(pl.check_replay, self.m.assign(occ=5)))

    def fly_one_badly(self, answer):
        """The last orbit, a poor one the median gate never sees, flown as `answer`."""
        m = self.m.assign(occ=[2] * (len(self.m) - 1) + [6])
        real = pl.fly_jpl
        last = m.iloc[-1].designation

        def fly(row):
            return answer if row.designation == last else real(row)

        with mock.patch.object(pl, "fly_jpl", fly):
            return quiet(pl.check_replay, m)

    def test_a_poor_orbit_may_be_wrong(self):
        self.assertTrue(self.fly_one_badly((30.0, 20.0)))

    def test_planted_a_solver_answer_no_orbit_gives_is_caught(self):
        # What the unguarded solver gave on a 25 day flight to 2021 GQ10.
        self.assertFalse(self.fly_one_badly((793.0, 790.0)))

    def test_planted_an_unsolved_flight_is_caught(self):
        self.assertFalse(self.fly_one_badly((np.nan, np.nan)))


class Coverage(unittest.TestCase):

    def test_a_full_answer_passes(self):
        self.assertTrue(quiet(pl.check_coverage, 400, []))

    def test_planted_a_throttled_run_is_refused(self):
        # The first live run: 243 of 415 lost to 403s and quietly dropped.
        self.assertFalse(quiet(pl.check_coverage, 415, ["HTTP 403"] * 243))


class Throttle(unittest.TestCase):

    def throttled(self, *codes):
        """A `snapshots.fetch` that raises each code in turn, then answers."""
        steps = iter(codes)

        def fetch(url, timeout):
            code = next(steps, None)
            if code is None:
                return b"{}"
            raise urllib.error.HTTPError(url, code, "x", None, None)
        return fetch

    def test_a_403_is_waited_out_and_retried(self):
        with mock.patch.object(snapshots, "fetch", self.throttled(403, 403)), \
                mock.patch.object(pl.time, "sleep") as sleep:
            self.assertEqual(quiet(pl.fetch_object, "u", 1.5, 120.0), "{}")
        self.assertEqual([c.args[0] for c in sleep.call_args_list], [1.5, 120.0, 240.0])

    def test_a_404_is_not_retried(self):
        with mock.patch.object(snapshots, "fetch", self.throttled(404)), \
                mock.patch.object(pl.time, "sleep"), self.assertRaises(urllib.error.HTTPError):
            quiet(pl.fetch_object, "u", 1.5, 120.0)

    def test_a_body_that_stays_throttled_is_counted_with_its_cause(self):
        with mock.patch.object(snapshots, "fetch", self.throttled(*[403] * 9)), \
                mock.patch.object(pl.time, "sleep"):
            legs, failed = quiet(pl.fetch_legs, ["2000 SG344"])
        self.assertTrue(legs.empty)
        self.assertEqual(failed, ["HTTP 403"])


class Oracles(unittest.TestCase):

    def setUp(self):
        self.m = recorded()

    def test_transfer_cost_is_jpl_leg_on_jpl_dates(self):
        r = next(self.m.itertuples())
        cost = pl.transfer_cost(r, r.launch_jd, r.dur_out)
        jpl = pl.leo_departure_dv_km_s(r.v_dep_earth) + r.vrel_arr_neo
        self.assertLess(abs(cost - jpl), 0.3)

    def test_dated_oracle_finds_jpl_leg_or_better(self):
        # It searches a superset of JPL's outbound legs, so it must not land
        # far above one. Coarse grid and two bodies, to keep this quick.
        for r in list(self.m.itertuples())[:2]:
            jpl = pl.leo_departure_dv_km_s(r.v_dep_earth) + r.vrel_arr_neo
            with mock.patch.object(pl, "LAUNCH_JD", (r.launch_jd - 60.0, r.launch_jd + 60.0)):
                self.assertLess(pl.dated_oracle(r, step_d=20.0), jpl + 0.3)

    def test_oracle_legs_price_what_the_oracle_returned(self):
        r = next(self.m.itertuples())
        best, at = pl.oracle_dv_out(r, n_dep=8, n_tof=8)
        v_inf, match, _ = pl.oracle_legs(r, at)
        self.assertAlmostEqual(pl.leo_departure_dv_km_s(v_inf) + match, best, places=6)

    def test_arrival_anomaly_at_the_apsides(self):
        row = next(self.m.itertuples())._replace(mean_anomaly_deg=0.0)
        deg, r, which = pl.arrival_anomaly(row, row.element_epoch_jd)
        self.assertAlmostEqual(deg, 0.0)
        self.assertAlmostEqual(r, row.semi_major_axis_au * (1 - row.eccentricity))
        self.assertEqual(which, "perihelion")

    def test_opik_u_vanishes_on_earths_own_orbit(self):
        self.assertAlmostEqual(float(pl.opik_u(1.0, 0.0, 0.0)), 0.0)


if __name__ == "__main__":
    unittest.main()
