# -*- coding: utf-8 -*-
"""JPL's outbound legs against the shipped estimator and F4's oracle.

F6 and F10 used NHATS's summary table, which gives one number per body, the
round-trip total. So those findings could rank an estimator and not size it,
and F10 could only suspect why the shipped model ranks at 0.60. NHATS's
object mode answers with the trajectory itself, leg by leg:

    launch, dur_out     when it leaves and how long the outbound leg flies
    v_dep_earth         departure v-infinity
    vrel_arr_neo        the rendezvous burn at the asteroid
    vrel_dep_neo        the burn leaving it
    vrel_arr_earth      the burn at Earth on return, usually zero

That makes three things measurable for the first time, on a sample of NHATS
bodies. The outbound leg can be repriced from the shipped LEO, so JPL can size
the estimator directly. The two legs can be compared one at a time, which is
what F10's match-burn weighting could only infer. And the arrival date says
where on its orbit JPL actually meets each body, which tests whether the
closed form's apsis rule is the problem.

The leg belongs to the cheapest ROUND trip of at most 450 days launched in
2020 to 2045. A one-way optimum can only be cheaper, so JPL's outbound leg is
an upper bound on the one-way cost and F4's oracle is the lower one. A body
the oracle prices above JPL's leg is a place the oracle is not a minimum.

F4's oracle has no calendar: it scans 40 departures across one synodic period,
capped at 12 years, from an arbitrary start. So a second one prices the same
objective on real dates, inside JPL's launch window and with Earth where it
was: if it ranks far closer to JPL than F4's does, phasing sets a ceiling that
no date-free estimator can reach.

Three checks run before anything is measured: the five docstring figures, the
Lambert solver's round trips, and three on JPL's answers. Their legs must sum to
their total, their departure burn must fit one parking orbit, which the API's
documentation does not name, and each must lie inside the window the dated
oracle searches. Then F4's Lambert solver has to fly JPL's
own trajectories: given only the launch date and flight time, it must arrive
at JPL's two legs from the catalog's elements. That pins the solver, the frame
and `element_epoch_jd` to something this repo did not build.

JPL answers an impatient client with HTTP 403, not 429, and the first run here
lost 243 of 415 bodies that way without noticing. So requests are paced, a
throttled one is waited out, and a run that still loses more than 2% stops.

This script READS ONLY. It fetches one public read-only API (recorded, or
replayed, through snapshots.py), does not import master, does not build a
stage, and does not write to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_nhats_legs.py
    py ../starred-repos/probe_nhats_legs.py --n 400 --seed 1
"""

import argparse
import datetime
import json
import time
import urllib.error
import urllib.parse
from collections import Counter

import numpy as np
import pandas as pd

import orbital as ob
import snapshots
from probe_lambert import EARTH, check_solver, oracle_dv_out
from probe_nhats import (
    DEFAULT_CATALOG,
    MU_EARTH_KM3_S2,
    NHATS_URL,
    V_EARTH_KM_S,
    check_reference,
    dv_out_shipped,
    fetch_nhats,
    leo_departure_dv_km_s,
)
from probe_plane_change import best_split
from probe_prospector import intercept_dv, shipped_legs

PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"
J2000_JD = 2451545.0
R_EARTH_KM = 6_378.14

# The parking orbit JPL's departure burn fits; see check_legs. The API's
# documentation does not name one.
JPL_PARK_ALT_KM = 400.0

THROTTLED = (403, 429, 503)

# NHATS's default search: launches from 2020 to 2045, at most 450 days in all,
# at least 8 of them at the body. That leaves 442 days of flight to split, and
# the outbound leg can take nearly all of it. The window runs through the end
# of 2045, not to its start: 2020-01-01 to 2046-01-01 in Julian dates. See
# check_window.
LAUNCH_JD = (2458849.5, 2468346.5)
MAX_TOF_D = 440.0

# Solar escape speed at 1 AU, 42.1 km/s, against Earth's own 29.8 head on: the
# fastest a bound transfer leaves. See check_replay.
MAX_VINF_KM_S = 72.0

LEG_COLUMNS = ("dv_total", "dv_dep_park", "vrel_arr_neo", "vrel_dep_neo", "vrel_arr_earth",
               "v_dep_earth", "c3", "dur_out", "dur_ret", "dur_total")
ELEMENTS = ["designation", "semi_major_axis_au", "eccentricity", "inclination_deg",
            "longitude_asc_node_deg", "arg_perihelion_deg", "mean_anomaly_deg",
            "element_epoch_jd"]


def julian_date(iso):
    """Julian date at 00:00 of a YYYY-MM-DD string."""
    return datetime.date.fromisoformat(iso).toordinal() + 1721424.5


def parse_object(text):
    """The legs of one object-mode answer's cheapest trajectory, or None if it has none."""
    d = json.loads(text)
    t = d.get("min_dv_traj")
    if not isinstance(t, dict):
        return None
    out = {k: float(t[k]) for k in LEG_COLUMNS}
    out["launch_jd"] = julian_date(t["launch"])
    out["occ"] = pd.to_numeric(d.get("occ"), errors="coerce")
    out["computed"] = d.get("computed")
    return out


def fetch_object(url, pause_s, backoff_s, tries=4):
    """One object-mode answer, waiting out JPL's rate limit.

    At two requests a second the API starts answering 403 after a couple of
    hundred, for bodies its own summary table lists. So a 403 here means
    throttled, not unknown, and a body must not be counted missing for it.
    """
    for attempt in range(tries):
        if not snapshots.replay_folder():
            time.sleep(pause_s if attempt == 0 else backoff_s * attempt)
        try:
            return snapshots.fetch(url, timeout=60).decode("utf-8", "replace")
        except urllib.error.HTTPError as err:
            if err.code not in THROTTLED or attempt == tries - 1:
                raise
            print(f"  throttled ({err.code}); waiting {backoff_s * (attempt + 1):.0f} s", flush=True)
    raise AssertionError("unreachable")


def fetch_legs(designations, pause_s=1.5, backoff_s=120.0):
    """Object-mode legs for each designation; (frame, what went wrong per failure)."""
    rows, failed = [], []
    for k, des in enumerate(designations, 1):
        url = f"{NHATS_URL}?des={urllib.parse.quote(des)}"
        try:
            legs = parse_object(fetch_object(url, pause_s, backoff_s))
        except urllib.error.HTTPError as err:
            legs, why = None, f"HTTP {err.code}"
        except (urllib.error.URLError, OSError, ValueError, KeyError) as err:
            legs, why = None, type(err).__name__
        else:
            why = "no trajectory"
        if legs is None:
            failed.append(why)
        else:
            rows.append(dict(legs, designation=des))
        if k % 100 == 0:
            print(f"  {k}/{len(designations)} ...", flush=True)
    return pd.DataFrame(rows), failed


def check_coverage(wanted, failed, tol=0.02):
    """Every body asked for is in JPL's own list, so nearly all must answer.

    A run that lost bodies to throttling is not the random sample it claims
    to be, and it would replay as if those bodies had no trajectory.
    """
    share = len(failed) / max(wanted, 1)
    ok = share <= tol
    kinds = ", ".join(f"{k} {v}" for k, v in sorted(Counter(failed).items())) or "none"
    print(f"COVERAGE  {wanted - len(failed):,} of {wanted:,} bodies answered; failures: {kinds}"
          f" ({'OK' if ok else 'MISMATCH -- this is not the sample it claims to be'})")
    return ok


def parking_dv(v_inf, alt_km):
    """Burn from a circular parking orbit at `alt_km` onto this v-infinity."""
    v_c = np.sqrt(MU_EARTH_KM3_S2 / (R_EARTH_KM + alt_km))
    return np.sqrt(2.0 * v_c ** 2 + np.asarray(v_inf, float) ** 2) - v_c


def check_legs(m, sum_tol=0.0021, park_tol=0.005):
    """JPL's legs must sum to its total, and its departure must fit one parking orbit.

    Each published figure is rounded to 0.001, so four of them can sum to at
    most 0.002 from a rounded total.
    """
    parts = m["dv_dep_park"] + m["vrel_arr_neo"] + m["vrel_dep_neo"] + m["vrel_arr_earth"]
    sum_gap = float((parts - m["dv_total"]).abs().max())
    park_gap = float(np.abs(parking_dv(m["v_dep_earth"], JPL_PARK_ALT_KM) - m["dv_dep_park"]).max())
    ok = sum_gap <= sum_tol and park_gap <= park_tol
    print(f"JPL LEGS  {len(m):,} trajectories: legs against dv_total, worst {sum_gap:.4f} km/s;"
          f" departure against a {JPL_PARK_ALT_KM:.0f} km parking orbit, worst {park_gap:.4f}"
          f" ({'OK' if ok else 'MISMATCH -- these are not the legs of the total'})")
    return ok


def check_window(m):
    """Every JPL leg must lie inside what the dated oracle searches.

    The dated oracle claims to search a superset of JPL's outbound legs. A
    first version stopped at 2045-01-01 and left 11 of 415 launches outside.
    """
    outside = (m["launch_jd"] < LAUNCH_JD[0]) | (m["launch_jd"] > LAUNCH_JD[1])
    long_leg = m["dur_out"] > MAX_TOF_D
    ok = not outside.any() and not long_leg.any()
    print(f"JPL WINDOW  {len(m):,} trajectories: {int(outside.sum())} launched outside the dated"
          f" oracle's window, {int(long_leg.sum())} with a longer outbound leg than it flies"
          f" ({'OK' if ok else 'MISMATCH -- the oracle cannot reach these legs'})")
    return ok


def earth_at(jd):
    """Keplerian Earth at a Julian date, on probe_lambert's J2000 elements."""
    r, v = ob.propagate_elements(EARTH["a"], EARTH["e"], EARTH["i"], EARTH["raan"], EARTH["argp"],
                                 EARTH["M0"], (jd - J2000_JD) * ob.DAY_S)
    return r.ravel(), v.ravel()


def body_at(row, jd):
    """Two-body state of a catalog row at a Julian date, from its own epoch."""
    r, v = ob.propagate_elements(
        float(row.semi_major_axis_au) * ob.AU_KM, float(row.eccentricity),
        np.radians(float(row.inclination_deg)), np.radians(float(row.longitude_asc_node_deg)),
        np.radians(float(row.arg_perihelion_deg)), np.radians(float(row.mean_anomaly_deg)),
        (jd - float(row.element_epoch_jd)) * ob.DAY_S)
    return r.ravel(), v.ravel()


def fly_jpl(row):
    """F4's solver on JPL's own launch date and flight time: (v_inf, arrival speed)."""
    t1 = float(row.launch_jd)
    t2 = t1 + float(row.dur_out)
    r_e, v_e = earth_at(t1)
    r_a, v_a = body_at(row, t2)
    v1, v2 = ob.lambert_izzo(r_e, r_a, float(row.dur_out) * ob.DAY_S)
    if v1 is None:
        return np.nan, np.nan
    return float(np.linalg.norm(v1 - v_e)), float(np.linalg.norm(v2 - v_a))


def check_replay(m, tol_km_s=0.05):
    """The solver, from the catalog's elements, must land on JPL's two legs.

    Only orbits with U <= 2 gate it: a worse orbit's catalog elements have
    moved since JPL computed the trajectory, which is the orbit's error and
    not the solver's. All of them are reported. Every orbit must still fly: an
    orbit's error moves an answer, but no error of position makes a transfer
    leave faster than a bound one can, or leaves it unsolved.
    """
    flown = [fly_jpl(r) for r in m.itertuples()]
    vinf = np.array([f[0] for f in flown])
    dep = np.abs(vinf - m["v_dep_earth"].to_numpy(float))
    arr = np.abs(np.array([f[1] for f in flown]) - m["vrel_arr_neo"].to_numpy(float))
    good = (m["occ"] <= 2).to_numpy()
    wild = ~(vinf < MAX_VINF_KM_S)
    print(f"REPLAY  F4's solver flown on JPL's dates, {len(m):,} trajectories: median |gap|"
          f" departure {np.nanmedian(dep):.3f}, arrival {np.nanmedian(arr):.3f} km/s")
    if not good.any():
        print("  no U <= 2 trajectory to gate on (MISMATCH -- the check has no data)")
        return False
    worst = max(float(np.nanmedian(dep[good])), float(np.nanmedian(arr[good])))
    ok = worst < tol_km_s
    print(f"  U <= 2, n={int(good.sum()):,}: departure {np.nanmedian(dep[good]):.3f},"
          f" arrival {np.nanmedian(arr[good]):.3f}, p90 {np.nanpercentile(dep[good], 90):.3f} and"
          f" {np.nanpercentile(arr[good], 90):.3f} km/s"
          f" ({'OK' if ok else 'MISMATCH -- the solver does not fly JPL trajectories'})")
    print(f"  every orbit: {int(wild.sum())} of {len(m):,} unsolved or leaving above the"
          f" {MAX_VINF_KM_S:.0f} km/s a bound transfer can reach"
          f" ({'OK' if not wild.any() else 'MISMATCH -- the solver is answering nonsense'})")
    return ok and not wild.any()


def arrival_anomaly(row, jd):
    """Where a body is at `jd`: (degrees from the nearest apsis, r in AU, which apsis)."""
    a = float(row.semi_major_axis_au) * ob.AU_KM
    e = float(row.eccentricity)
    n = np.sqrt(ob.MU_SUN_KM3_S2 / a ** 3)
    dt = (jd - float(row.element_epoch_jd)) * ob.DAY_S
    E = ob.kepler_E(e, np.radians(float(row.mean_anomaly_deg)) + n * dt)
    nu = np.degrees(np.abs(float(ob.true_anomaly(e, E))))
    return (min(nu, 180.0 - nu), float(row.semi_major_axis_au) * (1.0 - e * np.cos(float(E))),
            "perihelion" if nu < 90.0 else "aphelion")


def oracle_legs(row, at):
    """The legs of the oracle's best transfer, and its arrival date on the body's clock.

    `at` is what `oracle_dv_out` returned, (departure, flight time) in years.
    Its clock starts at J2000 for Earth and at the catalog's epoch for the
    body, because the oracle is date-free.
    """
    td, tf = at[0] * 365.25, at[1] * 365.25
    arrive = float(row.element_epoch_jd) + td + tf
    r_e, v_e = earth_at(J2000_JD + td)
    r_a, v_a = body_at(row, arrive)
    v1, v2 = ob.lambert_izzo(r_e, r_a, tf * ob.DAY_S)
    return float(np.linalg.norm(v1 - v_e)), float(np.linalg.norm(v2 - v_a)), arrive


def transfer_cost(row, jd, tof_d):
    """Outbound cost from the shipped LEO, launching at `jd` and flying `tof_d` days."""
    r_e, v_e = earth_at(jd)
    r_a, v_a = body_at(row, jd + tof_d)
    try:
        v1, v2 = ob.lambert_izzo(r_e, r_a, tof_d * ob.DAY_S)
    except (ValueError, FloatingPointError):
        return np.inf
    if v1 is None or not np.all(np.isfinite(v1)) or not np.all(np.isfinite(v2)):
        return np.inf
    return float(leo_departure_dv_km_s(np.linalg.norm(v1 - v_e)) + np.linalg.norm(v2 - v_a))


def dated_oracle(row, step_d=20.0):
    """F4's objective on real dates: JPL's launch window and flight-time limit.

    `oracle_dv_out` is date-free, so it can meet a body at a phase JPL's 2020
    to 2045 window never offers, and its 12 year cap can miss one that it does.
    This flies the same one-way leg from the catalog's real phase, on a coarse
    grid and then a fine one around the best cell, so what is left between it
    and JPL is not the calendar.
    """
    best, best_at = np.inf, None
    for t1 in np.arange(LAUNCH_JD[0], LAUNCH_JD[1], step_d):
        for tof in np.arange(step_d, MAX_TOF_D + 1, step_d):
            dv = transfer_cost(row, t1, tof)
            if dv < best:
                best, best_at = dv, (t1, tof)
    if best_at is None:
        return np.inf
    t0, f0 = best_at
    for t1 in np.linspace(t0 - step_d, t0 + step_d, 11):
        for tof in np.linspace(max(f0 - step_d, 2.0), f0 + step_d, 11):
            best = min(best, transfer_cost(row, t1, tof))
    return best


def shipped_winner_legs(m):
    """The shipped model's (v_inf, match) at whichever apsis it picks, per body."""
    a = m.semi_major_axis_au.to_numpy(float)
    e = m.eccentricity.to_numpy(float)
    i = m.inclination_deg.to_numpy(float)
    legs = shipped_legs(a, e, i)
    per = np.asarray(m["apsis"]) == "perihelion"
    return (np.where(per, legs["perihelion"][0], legs["aphelion"][0]),
            np.where(per, legs["perihelion"][1], legs["aphelion"][1]))


def price(m):
    """Every figure the report reads, one column per estimator."""
    a = m.semi_major_axis_au.to_numpy(float)
    e = m.eccentricity.to_numpy(float)
    i = m.inclination_deg.to_numpy(float)
    m = m.copy()
    m["shipped"], m["apsis"] = dv_out_shipped(a, e, i)
    m["ship_vinf"], m["ship_match"] = shipped_winner_legs(m)
    m["ship_depart"] = leo_departure_dv_km_s(m["ship_vinf"])
    # JPL's outbound leg repriced from the shipped 200 km LEO, so it is on the
    # shipped scale and the two compare as figures, not only as ranks.
    m["jpl_depart"] = leo_departure_dv_km_s(m["v_dep_earth"])
    m["jpl_out"] = m["jpl_depart"] + m["vrel_arr_neo"]
    m["split_c"] = np.minimum(best_split(a, i, a * (1.0 + e))[0],
                              best_split(a, i, a * (1.0 - e))[0])
    m["prospector"] = intercept_dv(a, e, i)
    m["opik_u"] = opik_u(a, e, i)
    return m


def opik_u(a, e, i_deg):
    """Encounter speed with Earth from the Tisserand parameter, km/s (Opik 1951)."""
    t = 1.0 / a + 2.0 * np.sqrt(a * (1.0 - e ** 2)) * np.cos(np.radians(i_deg))
    return np.sqrt(np.clip(3.0 - t, 0.0, None)) * V_EARTH_KM_S


def spearman(m, x, y):
    return m[[x, y]].corr(method="spearman").iloc[0, 1]


def report(m, label):
    """Sizes, legs, ranks and meeting points for one group of bodies."""
    gap = m["jpl_out"] - m["shipped"]
    orc = m["oracle"] - m["shipped"]
    print(f"\n{label}  (n = {len(m):,})")
    print("  outbound leg, km/s, all from the shipped 200 km LEO:")
    print(f"    shipped closed form   median {m.shipped.median():7.3f}")
    print(f"    F4's Lambert oracle   median {m.oracle.median():7.3f}"
          f"   oracle minus shipped median {orc.median():+.3f}")
    print(f"    the oracle on dates   median {m.dated.median():7.3f}"
          f"   dated minus shipped  median {(m.dated - m.shipped).median():+.3f}")
    print(f"    JPL's outbound leg    median {m.jpl_out.median():7.3f}"
          f"   JPL minus shipped    median {gap.median():+.3f}"
          f"  p10 {gap.quantile(0.1):+.3f}  p90 {gap.quantile(0.9):+.3f}")
    print(f"    shipped below JPL's leg on {(gap > 0).mean() * 100:.1f}%,"
          f" below the oracle on {(orc > 0).mean() * 100:.1f}%")
    over = m["oracle"] - m["jpl_out"]
    print(f"    oracle above JPL's leg (not a minimum there): {int((over > 0.01).sum()):,}"
          f" ({(over > 0.01).mean() * 100:.1f}%), median excess"
          f" {over[over > 0.01].median() if (over > 0.01).any() else 0.0:.3f} km/s")
    above = m["dated"] - m["jpl_out"]
    print(f"    dated oracle above JPL's leg: {int((above > 0.01).sum()):,}"
          f" ({(above > 0.01).mean() * 100:.1f}%);"
          f" a one-way leg costs JPL's minus dated, median {-above.median():+.3f} km/s less")

    print("  leg by leg, km/s:")
    print(f"    departure  shipped median {m.ship_depart.median():6.3f}"
          f"   JPL median {m.jpl_depart.median():6.3f}"
          f"   rho {spearman(m, 'ship_depart', 'jpl_depart'):.4f}")
    print(f"    rendezvous shipped median {m.ship_match.median():6.3f}"
          f"   JPL median {m.vrel_arr_neo.median():6.3f}"
          f"   rho {spearman(m, 'ship_match', 'vrel_arr_neo'):.4f}")
    print(f"    the oracle's own legs: departure median"
          f" {leo_departure_dv_km_s(m.oracle_vinf).median():.3f},"
          f" rendezvous median {m.oracle_match.median():.3f}"
          f"; its rendezvous against JPL's rho {spearman(m, 'oracle_match', 'vrel_arr_neo'):.4f}")

    print("  rank agreement (Spearman):")
    print(f"    {'estimator':16s} {'JPL outbound':>13s} {'JPL round trip':>15s}")
    for name, col in (("shipped", "shipped"), ("shipped depart", "ship_depart"),
                      ("shipped match", "ship_match"), ("split C", "split_c"),
                      ("prospector", "prospector"), ("Opik U", "opik_u"),
                      ("F4 oracle", "oracle"), ("dated oracle", "dated")):
        print(f"    {name:16s} {spearman(m, col, 'jpl_out'):13.4f}"
              f" {spearman(m, col, 'dv_total'):15.4f}")
    print(f"    {'JPL outbound':16s} {'':13s} {spearman(m, 'jpl_out', 'dv_total'):15.4f}")

    saving = m["shipped"] - m["split_c"]
    print(f"  split C saves median {saving.median():.3f} km/s, p90 {saving.quantile(0.9):.3f};"
          f" against JPL's leg it moves the median gap from {gap.median():+.3f}"
          f" to {(m['jpl_out'] - m['split_c']).median():+.3f}")

    near = m["apsis_deg"]
    print(f"  where JPL meets the body: degrees of true anomaly from the nearest apsis,"
          f" median {near.median():.1f} (uniform would be 45.0)")
    print(f"    within 30 deg of an apsis {(near < 30).mean() * 100:.1f}%"
          f" (uniform 33.3%); heliocentric distance median {m.r_arr.median():.3f} AU")
    same = (m["jpl_apsis"] == m["apsis"]).mean() * 100
    print(f"    JPL's nearest apsis is perihelion on {(m.jpl_apsis == 'perihelion').mean() * 100:.1f}%;"
          f" the shipped model picks perihelion on {(m.apsis == 'perihelion').mean() * 100:.1f}%;"
          f" they agree on {same:.1f}%")
    print(f"    F4's oracle, under no round-trip limit: nearest apsis is perihelion on"
          f" {(m.oracle_apsis == 'perihelion').mean() * 100:.1f}%,"
          f" {m.oracle_apsis_deg.median():.1f} deg from it at the median,"
          f" at {m.oracle_r_arr.median():.3f} AU; best flight"
          f" {m.oracle_tof_yr.median():.2f} yr against JPL's {m.dur_out.median() / 365.25:.2f}")


def grid_check(m, n=30, seed=1):
    """How far F4's 40 x 34 grid sits above a 80 x 68 one on these orbits."""
    s = m.sample(n=min(n, len(m)), random_state=seed)
    fine = np.array([oracle_dv_out(r, n_dep=80, n_tof=68)[0] for r in s.itertuples()])
    drop = s["oracle"].to_numpy(float) - fine
    print(f"\nGRID  F4's oracle against twice its resolution, n={len(s)}: median drop"
          f" {np.median(drop):.3f} km/s, p90 {np.percentile(drop, 90):.3f}, max {drop.max():.3f}")


def main():
    ap = argparse.ArgumentParser(description="JPL's outbound legs against the shipped estimator")
    ap.add_argument("--catalog", default=DEFAULT_CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--n", type=int, default=400, help="NHATS bodies to sample")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--top", type=int, default=100, help="ranking rows to add NHATS bodies from")
    ap.add_argument("--out", default=None, help="write the joined frame here")
    args = ap.parse_args()

    if not check_reference() or not check_solver():
        return
    print("\n" + snapshots.describe())

    cat = pd.read_csv(args.catalog, usecols=ELEMENTS + ["is_neo"], low_memory=False)
    cat["designation"] = cat["designation"].astype(str).str.strip()
    cat = cat[(cat["semi_major_axis_au"] > 0) & (cat["eccentricity"] >= 0)
              & (cat["eccentricity"] < 1)].dropna(subset=ELEMENTS)

    listed = fetch_nhats()
    joined = listed[["designation", "nhats_dv"]].merge(cat[ELEMENTS], on="designation", how="inner")
    print(f"NHATS  : {len(listed):,} NEAs from JPL; {len(joined):,} joined the catalog")

    sample = joined.sample(n=min(args.n, len(joined)), random_state=args.seed)
    top = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"],
                      nrows=args.top)
    top["designation"] = top["designation"].astype(str).str.strip()
    top_in = joined[joined["designation"].isin(top["designation"])]
    print(f"sample : {len(sample):,} (seed {args.seed}), plus {len(top_in)} from the top"
          f" {args.top} rows of the ranking")

    wanted = list(dict.fromkeys(list(sample["designation"]) + list(top_in["designation"])))
    legs, failed = fetch_legs(wanted)
    print()
    if not check_coverage(len(wanted), failed):
        return
    m = legs.merge(cat[ELEMENTS], on="designation", how="inner")

    if not check_legs(m) or not check_window(m) or not check_replay(m):
        return

    print("\nrunning both oracles; this is the slow part ...", flush=True)
    solved = [oracle_dv_out(r) for r in m.itertuples()]
    m["oracle"] = [s[0] for s in solved]
    m["oracle_at"] = [s[1] for s in solved]
    m = m[np.isfinite(m["oracle"])].reset_index(drop=True)
    legs_at = [oracle_legs(r, r.oracle_at) for r in m.itertuples()]
    m["oracle_tof_yr"] = [a[1] for a in m["oracle_at"]]
    m["oracle_vinf"] = [x[0] for x in legs_at]
    m["oracle_match"] = [x[1] for x in legs_at]
    there = [arrival_anomaly(r, x[2]) for r, x in zip(m.itertuples(), legs_at)]
    m["oracle_apsis_deg"] = [w[0] for w in there]
    m["oracle_r_arr"] = [w[1] for w in there]
    m["oracle_apsis"] = [w[2] for w in there]
    where = [arrival_anomaly(r, float(r.launch_jd) + float(r.dur_out)) for r in m.itertuples()]
    m["apsis_deg"] = [w[0] for w in where]
    m["r_arr"] = [w[1] for w in where]
    m["jpl_apsis"] = [w[2] for w in where]
    dated = []
    for k, r in enumerate(m.itertuples(), 1):
        dated.append(dated_oracle(r))
        if k % 100 == 0:
            print(f"  dated oracle {k}/{len(m)} ...", flush=True)
    m["dated"] = dated
    m = price(m[np.isfinite(m["dated"])])

    report(m[m["designation"].isin(sample["designation"])], "NHATS SAMPLE")
    in_top = m[m["designation"].isin(top_in["designation"])]
    if len(in_top) >= 5:
        report(in_top, f"TOP {args.top} ROWS OF THE RANKING, IN NHATS")
    grid_check(m)

    if args.out:
        m.to_csv(args.out, index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
