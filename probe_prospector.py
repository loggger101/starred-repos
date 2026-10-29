# -*- coding: utf-8 -*-
"""Prospector's closed form against the shipped one, with NHATS as referee.

THIRD-PASS.md found that Karmanplus/prospector (Apache-2.0) carries an
impulsive estimator, `intercept_dv`, built the way economicspace's is: a
Hohmann half-ellipse from 1 AU to an apsis of the target, and one burn to match
it there. It differs in two places, and FINDINGS.md has argued about both:

    PLANE CHANGE   prospector charges the inclination on the match burn, at the
                   target's local speed. economicspace charges it on departure.
                   In probe_plane_change.py's terms that is placement B against
                   the shipped placement A, and F4 measured C, the best split
                   between them, making the estimator worse.
    MEETING POINT  prospector also tries the 1 AU crossing, for orbits with
                   q <= 1 <= Q. economicspace prices the two apsides only.

Four estimators, all one-way:

    shipped         placement A, both apsides, departure from LEO. A copy of
                    `_transfer_legs_for_apsis` (from probe_nhats.py).
    prospector      `intercept_dv` as upstream wrote it: heliocentric, from
                    Earth's orbit with no escape from LEO, three meeting points.
                    Its figures are not on the shipped scale, so only its ranks
                    compare.
    B_leo           prospector's construction at the two apsides with the
                    departure burn priced from LEO and the shipped constants.
                    Against `shipped`, the plane change is the only difference.
    B_leo_cross     B_leo plus the 1 AU crossing. Against `B_leo`, the crossing
                    is the only difference.

NHATS is a round trip, so it ranks the estimators and cannot size them. Two
checks run before anything is measured: the five figures in
`asteroid_transfer_dv_km_s`'s docstring, which pin `shipped` to the shipped
model, and prospector's own `intercept_dv` on the same five bodies, which pin
the adapted copy to upstream.

`intercept_dv` is adapted from prospector's `prospector/solvers/edelbaum.py`;
CITATIONS.md section 3 carries its notice and says what changed.

This script READS ONLY. It fetches one public read-only API, does not import
master, does not build a stage, and does not write to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_prospector.py
    py ../starred-repos/probe_prospector.py --catalog path/to/asteroid_catalog.csv
"""

import argparse

import numpy as np
import pandas as pd

from probe_nhats import (
    DEFAULT_CATALOG,
    REFERENCE,
    V_EARTH_KM_S,
    check_reference,
    dv_out_shipped,
    fetch_nhats,
    leo_departure_dv_km_s,
    load_catalog,
)

# Mirrored from prospector/constants.py (0.4.3), for the `prospector` column
# only; B_leo and B_leo_cross use the shipped V_EARTH_KM_S.
MU_SUN_KM3_S2 = 1.32712440018e11
AU_KM         = 1.495978707e8

# prospector 0.4.3's own `intercept_dv` on REFERENCE's five bodies, run from
# its source on 2026-10-29. The adapted copy below must reproduce them.
UPSTREAM_INTERCEPT = np.array([11.460480, 6.171367, 4.709884, 7.586167, 4.115497])

BANDS = (("i < 5 deg", 0.0, 5.0), ("5-15 deg", 5.0, 15.0), ("i > 15 deg", 15.0, 180.0))


def meeting_legs(a, e, i_deg, v1):
    """Prospector's three meetings, as (depart, match) pairs in units of km/s.

    Adapted from `intercept_dv`. Changed from upstream: it returns the legs
    rather than the cheapest total, so the departure can be priced either
    heliocentrically or from LEO; the circular speed `v1` is a parameter; and
    the departure orbit is fixed at a0 = 1 AU, i0 = 0, upstream's defaults.
    The crossing has no match burn: it departs straight onto the target's
    velocity, and its leg is inf where the orbit does not reach 1 AU.
    """
    a = np.asarray(a, float)
    e = np.asarray(e, float)
    cdi = np.cos(np.radians(np.abs(np.asarray(i_deg, float))))
    q = a * (1 - e)
    Q = a * (1 + e)

    def speed_at(r, sma):
        return v1 * np.sqrt(np.clip(2.0 / r - 1.0 / sma, 0.0, None))

    def at_apse(R):
        a_t = 0.5 * (1.0 + R)
        depart = np.abs(speed_at(1.0, a_t) - v1)
        v_transfer = speed_at(R, a_t)
        v_target = speed_at(R, a)
        match = np.sqrt(np.clip(
            v_transfer ** 2 + v_target ** 2 - 2 * v_transfer * v_target * cdi, 0.0, None))
        return depart, match

    crosses = (q <= 1.0) & (1.0 <= Q)
    v_target_at_1 = speed_at(1.0, a)
    v_target_tangential = v1 * np.sqrt(np.clip(a * (1 - e ** 2), 0.0, None))
    cross = np.sqrt(np.clip(
        v_target_at_1 ** 2 + v1 ** 2 - 2 * v1 * v_target_tangential * cdi, 0.0, None))
    cross = np.where(crosses, cross, np.inf)

    return {"perihelion": at_apse(q), "aphelion": at_apse(Q), "crossing": (cross, 0.0)}


def intercept_dv(a, e, i_deg):
    """Prospector's estimator as upstream computes it: heliocentric, cheapest of three."""
    legs = meeting_legs(a, e, i_deg, np.sqrt(MU_SUN_KM3_S2 / AU_KM))
    return np.minimum.reduce([d + m for d, m in legs.values()])


def b_leo(a, e, i_deg, crossing):
    """Prospector's construction with the departure priced from LEO.

    Returns (dv, which_meeting). On the shipped constants, so against
    `dv_out_shipped` the plane change, and with `crossing` the extra meeting
    point, are the only differences.
    """
    legs = meeting_legs(a, e, i_deg, V_EARTH_KM_S)
    names = ["aphelion", "perihelion"] + (["crossing"] if crossing else [])
    totals = np.vstack([leo_departure_dv_km_s(legs[n][0]) + legs[n][1] for n in names])
    k = np.argmin(totals, axis=0)
    return totals[k, np.arange(totals.shape[1])], np.asarray(names)[k]


def shipped_legs(a, e, i_deg):
    """The shipped model's (v_inf, match) per apsis, in km/s.

    Placement A: the whole plane change in the departure v_inf by the law of
    cosines, and a coplanar match, which is prospector's match at i = 0.
    """
    a = np.asarray(a, float)
    e = np.asarray(e, float)
    flat = meeting_legs(a, e, 0.0, V_EARTH_KM_S)
    out = {}
    for name, R in (("aphelion", a * (1 + e)), ("perihelion", a * (1 - e))):
        v_t = np.sqrt(np.clip(2.0 - 2.0 / (1.0 + R), 1e-12, None))
        v_inf = np.sqrt(np.clip(
            v_t ** 2 + 1.0 - 2.0 * v_t * np.cos(np.radians(i_deg)), 0.0, None)) * V_EARTH_KM_S
        out[name] = (v_inf, flat[name][1])
    return out


def check_upstream():
    """The adapted `intercept_dv` must reproduce prospector's own figures."""
    r = REFERENCE[["body", "semi_major_axis_au", "eccentricity", "inclination_deg"]].copy()
    r["upstream"] = UPSTREAM_INTERCEPT
    r["adapted"] = intercept_dv(r.semi_major_axis_au, r.eccentricity, r.inclination_deg)
    r["delta"] = r["adapted"] - r["upstream"]

    print("\nUPSTREAM  (the adapted intercept_dv must reproduce prospector 0.4.3's)")
    print(r.to_string(index=False, float_format=lambda x: f"{x:10.6f}"))
    worst = r["delta"].abs().max()
    ok = worst < 1e-5
    print(f"\n  worst |adapted - upstream| = {worst:.2e} km/s "
          f"({'OK' if ok else "MISMATCH -- this is not prospector's estimator"})")

    # With no inclination there is no plane change to place, so B_leo must be
    # the shipped model exactly. If it is not, something besides the placement
    # differs and the comparison below would not isolate F1.
    flat = np.zeros(len(r))
    gap = np.abs(b_leo(r.semi_major_axis_au, r.eccentricity, flat, crossing=False)[0]
                 - dv_out_shipped(r.semi_major_axis_au.to_numpy(float),
                                  r.eccentricity.to_numpy(float), flat)[0]).max()
    flat_ok = gap < 1e-9
    print(f"  at i = 0, worst |B_leo - shipped| = {gap:.2e} km/s "
          f"({'OK' if flat_ok else 'MISMATCH -- B_leo differs by more than the placement'})")
    return ok and flat_ok


def spearman(m, x, y="nhats_dv"):
    return m[[x, y]].corr(method="spearman").iloc[0, 1]


def report(m):
    """Rank agreement with JPL, overall and in F4's inclination bands."""
    cols = ["shipped", "prospector", "B_leo", "B_leo_cross"]
    rt = m["nhats_dv"]

    print(f"\nNHATS  (n = {len(m):,} joined)")
    print(f"  {'estimator':12s} {'median':>8s} {'rho':>8s}"
          + "".join(f" {lab:>12s}" for lab, _, _ in BANDS)
          + f" {'> round trip':>14s}")
    print("  " + " " * 12 + " " * 18
          + "".join(f" {'n=' + format(int(((m.inclination_deg >= lo) & (m.inclination_deg < hi)).sum()), ','):>12s}"
                    for _, lo, hi in BANDS))
    for c in cols:
        row = f"  {c:12s} {m[c].median():8.3f} {spearman(m, c):8.4f}"
        for _, lo, hi in BANDS:
            s = m[(m.inclination_deg >= lo) & (m.inclination_deg < hi)]
            row += f" {spearman(s, c):12.4f}"
        over = "-" if c == "prospector" else f"{int((m[c] > rt).sum()):,}"
        print(row + f" {over:>14s}")

    print("\n  meeting point chosen by B_leo_cross:")
    for name, n in m["meet_cross"].value_counts().items():
        print(f"    {name:11s} {n:6,}  ({n / len(m) * 100:.1f}%)")

    gain = m["B_leo"] - m["B_leo_cross"]
    won = gain > 0
    if won.any():
        print(f"  where the crossing wins it saves median {gain[won].median():.3f}"
              f" km/s, max {gain.max():.3f}")

    print("\n  agreement between the estimators (Spearman):")
    for x, y in (("shipped", "B_leo"), ("shipped", "prospector"), ("B_leo", "B_leo_cross")):
        print(f"    {x:10s} vs {y:12s} {spearman(m, x, y):.4f}")

    d = m["B_leo"] - m["shipped"]
    print(f"\n  B_leo minus shipped : median {d.median():+.3f}  mean {d.mean():+.3f}"
          f"  p10 {np.percentile(d, 10):+.3f}  p90 {np.percentile(d, 90):+.3f} km/s")
    print(f"  B_leo cheaper than shipped on {int((d < 0).sum()):,} of {len(m):,} bodies")


def diagnose(m):
    """Which of prospector's differences carries its rank agreement.

    Its estimator differs from the shipped one in three ways: the plane change,
    the crossing, and a heliocentric departure where the shipped prices it from
    LEO. B_leo and B_leo_cross isolate the first two; this isolates the third,
    then sweeps the weight on the match burn to show what the departure pricing
    is doing to the balance between the two legs. A weight of 2 is roughly a
    round trip, which NHATS is.
    """
    a = m.semi_major_axis_au.to_numpy(float)
    e = m.eccentricity.to_numpy(float)
    i = m.inclination_deg.to_numpy(float)
    rt = m["nhats_dv"]
    A = shipped_legs(a, e, i)
    B = meeting_legs(a, e, i, np.sqrt(MU_SUN_KM3_S2 / AU_KM))
    apsides = ("aphelion", "perihelion")

    def rho(x):
        return pd.Series(x, index=m.index).corr(rt, method="spearman")

    rebuilt = np.minimum.reduce([leo_departure_dv_km_s(A[n][0]) + A[n][1] for n in apsides])
    gap = np.abs(rebuilt - m["shipped"].to_numpy(float)).max()
    print(f"\nDIAGNOSIS  (shipped rebuilt from its legs: worst gap {gap:.1e} km/s)")

    print("  departure priced heliocentrically, apsides only:")
    print(f"    placement A  {rho(np.minimum.reduce([A[n][0] + A[n][1] for n in apsides])):.4f}")
    print(f"    placement B  {rho(np.minimum.reduce([B[n][0] + B[n][1] for n in apsides])):.4f}"
          f"   (prospector without the crossing)")

    print("  shipped from LEO, match burn weighted by w:")
    for w in (0.25, 0.5, 0.75, 1.0, 2.0):
        x = np.minimum.reduce([leo_departure_dv_km_s(A[n][0]) + w * A[n][1] for n in apsides])
        print(f"    w = {w:4.2f}   {rho(x):.4f}")


def main():
    """Both reference checks, then the four estimators against NHATS."""
    ap = argparse.ArgumentParser(description="prospector vs shipped, refereed by NHATS")
    ap.add_argument("--catalog", default=DEFAULT_CATALOG)
    ap.add_argument("--out", default=None, help="write the joined frame here")
    args = ap.parse_args()

    if not check_reference() or not check_upstream():
        return

    cat = load_catalog(args.catalog)
    print(f"\ncatalog: {len(cat):,} usable rows")
    n = fetch_nhats()
    print(f"NHATS  : {len(n):,} NEAs from JPL")
    m = n.merge(cat, on="designation", how="inner")
    print(f"joined : {len(m):,}  ({len(m) / len(n) * 100:.1f}% of NHATS)")

    a = m.semi_major_axis_au.to_numpy(float)
    e = m.eccentricity.to_numpy(float)
    i = m.inclination_deg.to_numpy(float)
    m["shipped"], m["apsis"] = dv_out_shipped(a, e, i)
    m["prospector"] = intercept_dv(a, e, i)
    m["B_leo"], m["meet"] = b_leo(a, e, i, crossing=False)
    m["B_leo_cross"], m["meet_cross"] = b_leo(a, e, i, crossing=True)
    m = m[np.isfinite(m[["shipped", "prospector", "B_leo", "B_leo_cross"]]).all(axis=1)]
    report(m)
    diagnose(m)

    if args.out:
        m.to_csv(args.out, index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
