# -*- coding: utf-8 -*-
"""How much of F4's gap is the oracle's own two-body error?

F4 measured the shipped closed form about 1.3 km/s below a Lambert porkchop
oracle. SPACE-MAP.md section 4 pointed out that the oracle propagates every
asteroid two-body over flight times of up to four years, so it carries an
error nobody had measured, in both directions. SECOND-PASS.md proposed
bounding it with brahe and the sb441 kernels, which reach only the 373 most
massive asteroids. THIRD-PASS.md proposed ASSIST, which needs a Linux build.

JPL Horizons does the integration itself. For any body it serves heliocentric
state vectors from its full n-body solution (planets, Moon, large asteroids,
relativity), so the oracle can be run twice on the same grid: once with the
target on Horizons' trajectory, and once two-body from Horizons' osculating
elements at the start of the window. The difference is the oracle's two-body
error, body by body, with nothing else changed. Earth is Horizons' in both runs.

The sample is two groups:

    F4 sample   the first 40 bodies of probe_lambert's seed-1 sample of 400
                on today's catalog, which is mostly main belt like F4's
    NEOs        30 near-Earth objects, seed 1, where planetary encounters
                bend orbits most

Each body's window opens on 2027-01-01 TDB and runs for one synodic period
of departures (capped at 12 years, as probe_lambert does) plus four years of
flight. Horizons vectors are fetched every two days and interpolated with
cubic Hermite splines, which use the velocities as well as the positions.

Checks before measuring: probe_lambert's docstring and solver checks; the
interpolation, measured on nodes a sparser track never saw at four- and
eight-day spacing and carried down to the two days used by its fourth-power
scaling, which must be too small to move a 0.25-year Lambert solve by 0.001
km/s; and two-body propagation from Horizons' osculating elements, which must land on
Horizons' own position at the start of the window to 1 km.

Network reads go through snapshots.py (about 28 MB for the default sample).
This script READS ONLY and writes nothing to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_horizons.py
    py ../starred-repos/probe_horizons.py --main 80 --neo 60
"""

import argparse
import io
import time
import urllib.parse

import numpy as np
import pandas as pd

import orbital as ob
import snapshots
from probe_lambert import check_reference, check_solver, leo_departure_dv_km_s

API = "https://ssd.jpl.nasa.gov/api/horizons.api"
START = "2027-01-01"
JD_START = 2461406.5                     # 2027-01-01 00:00 TDB
CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PAUSE_S = 0.5
STEP_D = 2


def horizons(command, ephem, stop, step, units):
    q = {"format": "text", "COMMAND": f"'{command}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": f"'{ephem}'", "CENTER": "'500@10'", "START_TIME": f"'{START}'",
         "STOP_TIME": f"'{stop}'", "STEP_SIZE": f"'{step}'", "REF_PLANE": "'ECLIPTIC'",
         "REF_SYSTEM": "'J2000'", "OUT_UNITS": f"'{units}'", "CSV_FORMAT": "'YES'"}
    if ephem == "VECTORS":
        q.update({"VEC_TABLE": "'2'", "VEC_LABELS": "'NO'"})
    text = snapshots.fetch(API + "?" + urllib.parse.urlencode(q), timeout=300).decode("utf-8", "replace")
    return parse_horizons(text)


def parse_horizons(text):
    """The CSV block between $$SOE and $$EOE: JD then the numeric columns (the date dropped)."""
    if "$$SOE" not in text:
        raise LookupError(text[:300])
    body = text.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    rows = [r.rstrip(",").split(",") for r in body.strip().splitlines()]
    return np.array([[float(r[0])] + [float(x) for x in r[2:]] for r in rows])


def command_for(designation):
    d = str(designation).strip()
    return f"{d};" if d.isdigit() else f"DES={d};"


def stop_date(days):
    return (pd.Timestamp(START) + pd.Timedelta(days=int(np.ceil(days)) + STEP_D)).strftime("%Y-%m-%d")


class Track:
    """A Horizons vector table, interpolated with cubic Hermite splines."""

    def __init__(self, table):
        self.t = (table[:, 0] - JD_START) * ob.DAY_S
        self.r = table[:, 1:4]
        self.v = table[:, 4:7]

    def state(self, t_s):
        t_s = np.atleast_1d(np.asarray(t_s, float))
        k = np.clip(np.searchsorted(self.t, t_s) - 1, 0, len(self.t) - 2)
        h = (self.t[k + 1] - self.t[k])[:, None]
        s = ((t_s - self.t[k]) / h[:, 0])[:, None]
        r0, r1, v0, v1 = self.r[k], self.r[k + 1], self.v[k], self.v[k + 1]
        h00, h10, h01, h11 = 2 * s**3 - 3 * s**2 + 1, s**3 - 2 * s**2 + s, -2 * s**3 + 3 * s**2, s**3 - s**2
        r = h00 * r0 + h10 * h * v0 + h01 * r1 + h11 * h * v1
        d00, d10, d01, d11 = 6 * s**2 - 6 * s, 3 * s**2 - 4 * s + 1, -6 * s**2 + 6 * s, 3 * s**2 - 2 * s
        v = (d00 * r0 + d10 * h * v0 + d01 * r1 + d11 * h * v1) / h
        return r, v


def interpolation_gap(table, stride=2):
    """Worst position error (km) at nodes a track built from every `stride`-th node never saw."""
    sparse = Track(table[::stride])
    keep = np.ones(len(table), bool)
    keep[::stride] = False
    keep[-stride:] = False
    held = table[keep]
    r, _ = sparse.state((held[:, 0] - JD_START) * ob.DAY_S)
    return float(np.max(np.linalg.norm(r - held[:, 1:4], axis=1)))


def interpolation_error_km(table):
    """The interpolation error at the spacing actually used, extrapolated down.

    A cubic Hermite spline's error goes as the fourth power of the spacing, so
    the error is measured at twice and four times the spacing used (on nodes the
    sparser track never saw) and carried down one more halving by the ratio
    between them. The ratio is returned too: near 16 says the fourth-power
    regime holds, which is what licenses the extrapolation.
    """
    g2, g4 = interpolation_gap(table, 2), interpolation_gap(table, 4)
    ratio = g4 / g2 if g2 > 0 else np.inf
    return g2 / max(ratio, 1.0), ratio


def elements_from(row):
    """Horizons ELEMENTS row (AU-D): (a km, e, i, raan, argp, M0) in radians."""
    _, ec, _, inc, om, w, _, _, ma, _, a, _, _ = row
    return (a * ob.AU_KM, ec, np.radians(inc), np.radians(om), np.radians(w), np.radians(ma))


def porkchop(earth, target_state, syn_yr, n_dep=40, n_tof=34, tof_lo_yr=0.25, tof_hi_yr=4.0):
    """probe_lambert.oracle_dv_out's grid and objective, with the states supplied."""
    t_dep = np.linspace(0.0, syn_yr * 365.25 * ob.DAY_S, n_dep, endpoint=False)
    tofs = np.linspace(tof_lo_yr, tof_hi_yr, n_tof) * 365.25 * ob.DAY_S
    best = np.inf
    rE_all, vE_all = earth.state(t_dep)
    for td, rE, vE in zip(t_dep, rE_all, vE_all):
        rA_all, vA_all = target_state(td + tofs)
        for tf, rA, vA in zip(tofs, rA_all, vA_all):
            try:
                v1, v2 = ob.lambert_izzo(rE, rA, tf)
            except (ValueError, FloatingPointError):
                continue
            if v1 is None or not (np.all(np.isfinite(v1)) and np.all(np.isfinite(v2))):
                continue
            dv = leo_departure_dv_km_s(np.linalg.norm(v1 - vE)) + np.linalg.norm(v2 - vA)
            best = min(best, dv)
    return best


def synodic_yr(a_au):
    p = a_au ** 1.5
    denom = abs(1.0 / p - 1.0)
    return min(1.0 / denom, 12.0) if denom > 1e-9 else 12.0


def sample(path, n_main, n_neo, seed):
    cols = ["designation", "semi_major_axis_au", "eccentricity", "inclination_deg",
            "longitude_asc_node_deg", "arg_perihelion_deg", "mean_anomaly_deg", "is_neo"]
    df = pd.read_csv(path, usecols=cols, low_memory=False)
    df = df[(df["semi_major_axis_au"] > 0.3) & (df["semi_major_axis_au"] < 6.0)
            & (df["eccentricity"] >= 0) & (df["eccentricity"] < 0.9)].dropna(subset=cols[1:7])
    df["designation"] = df["designation"].astype(str).str.strip()
    f4 = df.sample(n=min(400, len(df)), random_state=seed).head(n_main).assign(group="F4 sample")
    neo = df[df["is_neo"].astype(str).str.lower().isin(["true", "1"])]
    neo = neo[~neo.designation.isin(f4.designation)].sample(n=n_neo, random_state=seed)
    return pd.concat([f4, neo.assign(group="NEOs")], ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description="the two-body error in F4's Lambert oracle")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--main", type=int, default=40)
    ap.add_argument("--neo", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    if not check_reference() or not check_solver():
        return
    print("\n" + snapshots.describe())
    bodies = sample(args.catalog, args.main, args.neo, args.seed)
    span = max(synodic_yr(a) for a in bodies.semi_major_axis_au) * 365.25 + 4.0 * 365.25
    earth_table = horizons("399", "VECTORS", stop_date(span), f"{STEP_D} d", "KM-S")
    earth = Track(earth_table)
    gap_e, ratio_e = interpolation_error_km(earth_table)
    print(f"Earth: {len(earth_table):,} Horizons vectors; interpolation error at {STEP_D}-day spacing"
          f" {gap_e:.3f} km (ratio {ratio_e:.1f})")

    rows, worst_interp, worst_start, failed = [], gap_e, 0.0, []
    ratios = [ratio_e]
    for b in bodies.itertuples():
        syn = synodic_yr(b.semi_major_axis_au)
        stop = stop_date(syn * 365.25 + 4.0 * 365.25)
        try:
            vec = horizons(command_for(b.designation), "VECTORS", stop, f"{STEP_D} d", "KM-S")
            el = horizons(command_for(b.designation), "ELEMENTS", "2027-01-02", "1 d", "AU-D")[0]
        except (LookupError, OSError, ValueError) as err:
            failed.append((b.designation, str(err)[:80]))
            continue
        finally:
            if not snapshots.replay_folder():
                time.sleep(PAUSE_S)
        track = Track(vec)
        elements = elements_from(el)
        err, ratio = interpolation_error_km(vec)
        worst_interp = max(worst_interp, err)
        ratios.append(ratio)

        def two_body(t_s, elements=elements):
            r, v = ob.propagate_elements(*elements, np.asarray(t_s, float))
            return r.T, v.T

        r0, _ = two_body([0.0])
        worst_start = max(worst_start, float(np.linalg.norm(r0[0] - vec[0, 1:4])))
        t_all = track.t
        drift = np.linalg.norm(two_body(t_all)[0] - track.r, axis=1) / ob.AU_KM
        true_dv = porkchop(earth, track.state, syn)
        kep_dv = porkchop(earth, two_body, syn)
        rows.append({"designation": b.designation, "group": b.group,
                     "a": b.semi_major_axis_au, "e": b.eccentricity, "i": b.inclination_deg,
                     "years": t_all[-1] / ob.DAY_S / 365.25, "drift_max_au": drift.max(),
                     "oracle_two_body": kep_dv, "oracle_horizons": true_dv,
                     "delta": true_dv - kep_dv})

    r = pd.DataFrame(rows)
    # A position error spread over the shortest flight time the grid flies is
    # the velocity error it can put into a Lambert solve.
    dv_err = worst_interp / (0.25 * 365.25 * ob.DAY_S)
    ok_i = dv_err < 1e-3
    ok_s = worst_start < 1.0
    print(f"\nCHECKS  interpolation: worst error at {STEP_D}-day spacing {worst_interp:.1f} km,"
          f" so at most {dv_err:.1e} km/s in a 0.25-year Lambert solve; fourth-power ratio"
          f" median {np.median(ratios):.1f} ({'OK' if ok_i else 'MISMATCH -- interpolation too coarse'})")
    print(f"        two-body from Horizons' elements at the window start: worst {worst_start:.4f} km"
          f" ({'OK' if ok_s else 'MISMATCH -- frame, units or epoch disagree'})")
    print(f"  bodies solved {len(r)}; Horizons refused {len(failed)}"
          + (f": {failed[:3]}" if failed else ""))
    if not (ok_i and ok_s) or r.empty:
        return

    print("\nTHE ORACLE'S TWO-BODY ERROR  (oracle with Horizons' trajectory minus two-body, km/s)")
    for g, s in list(r.groupby("group")) + [("all", r)]:
        d = s.delta[np.isfinite(s.delta)]
        print(f"  {g:10s} n={len(d):3}  median {d.median():+.4f}  mean {d.mean():+.4f}"
              f"  median |delta| {d.abs().median():.4f}  p90 |delta| {d.abs().quantile(0.9):.4f}"
              f"  max |delta| {d.abs().max():.3f}")
        print(f"  {'':10s}        two-body position drift over the window: median"
              f" {s.drift_max_au.median():.4f} AU, max {s.drift_max_au.max():.3f} AU")
    print("\n  per body:")
    print(r.sort_values("delta", key=np.abs, ascending=False).to_string(
        index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
