# -*- coding: utf-8 -*-
"""Will anyone be able to find the bodies at the top of the ranking again?

F14 found that orbit uncertainty barely moves delta-v for the bodies at the
top of the profitability ranking, and that what it moves is where they are.
A mission needs the body to be found again before launch, to refine its orbit,
and a metre-scale body is only visible during a close approach. This probe
asks, for each NEO in the top rows: in the years before 2036, is there a
window when the body is bright enough and far enough from the Sun to observe,
and if so, how large a patch of sky would a telescope have to search?

For each body it takes NEODyS's Keplerian elements and full 6x6 covariance
(through probe_covariance.fetch, so recorded or replayed by snapshots.py), and:

    window      propagates the nominal orbit two-body, day by day, against a
                Keplerian Earth, and finds days with V <= 24 (a Rubin-class
                limit, optimistic for today's surveys) at solar elongation
                >= 60 degrees. V uses the H, G magnitude system with G = 0.15
    p_observable  the share of draws from the covariance that have any such day
    search radius the 95th percentile of the angular distance, as seen from
                Earth, between the draws and the nominal position on the
                nominal orbit's brightest window day

    recoverable   a window exists and the search radius is under 1 degree
    search        a window exists and the radius is 1 to 10 degrees
    lost          a window exists and the radius is over 10 degrees
    no window     the nominal orbit is never observable before 2036

Two-body propagation ignores the planets, and every one of these bodies makes
close Earth approaches, which stretch an orbit's uncertainty further. So the
search radii here are lower bounds: a body this probe calls lost is lost, and
one it calls recoverable may not be.

Checks before measuring: the Keplerian Earth must sit at heliocentric
ecliptic longitude 100.376 degrees at J2000 (the Sun's true longitude there,
less 180), and on each body the linear growth of the mean anomaly's
uncertainty over ten years must agree with the spread of 4,000 draws to 5%. The 10-year drift
of the mean anomaly is also compared with the catalog's condition code U,
which the MPC defines from that drift.

This script READS ONLY. It reads the catalog and the profitability catalog,
fetches NEODyS (recorded or replayed), and writes nothing to
`asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_recovery.py
    py ../starred-repos/probe_recovery.py --top 1000 --draws 300
"""

import argparse
import datetime
import time

import numpy as np
import pandas as pd

import orbital as ob
import snapshots
from probe_covariance import PAUSE_S, fetch
from probe_lambert import EARTH

CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"

MJD_J2000 = 51544.5
END = datetime.date(2036, 1, 1)
V_LIMIT = 24.0
MIN_ELONGATION_DEG = 60.0
G_SLOPE = 0.15
GAUSS_N_DEG_DAY = 0.9856076686           # mean motion at 1 AU, degrees per day
# The Sun's true geocentric longitude at J2000 is its mean longitude 280.460
# plus the equation of centre, 1.915 sin(357.529) + 0.020 sin(715.058) =
# -0.084, so 280.376; the Earth's heliocentric longitude is 180 less.
EARTH_LON_J2000_DEG = 100.376


def mjd(d):
    return (d - datetime.date(1858, 11, 17)).days


def earth_state(mjd_days):
    dt = (np.asarray(mjd_days, float) - MJD_J2000) * ob.DAY_S
    r, _ = ob.propagate_elements(EARTH["a"], EARTH["e"], EARTH["i"], EARTH["raan"],
                                 EARTH["argp"], EARTH["M0"], dt)
    return r / ob.AU_KM


def body_positions(kep, epoch_mjd, mjd_days):
    """Heliocentric positions (AU) for element rows `kep` (n, 6; deg) at `mjd_days`.

    Returns shape (3, n, t).
    """
    kep = np.atleast_2d(kep)
    a = kep[:, 0:1] * ob.AU_KM
    rad = np.radians(kep[:, 2:6])
    dt = (np.asarray(mjd_days, float)[None, :] - epoch_mjd) * ob.DAY_S
    r, _ = ob.propagate_elements(a, kep[:, 1:2], rad[:, 0:1], rad[:, 1:2], rad[:, 2:3],
                                 rad[:, 3:4], dt)
    return r / ob.AU_KM


def phase_function(alpha):
    """H, G system (Bowell et al. 1989), G = G_SLOPE."""
    t = np.tan(alpha / 2)
    phi1 = np.exp(-3.33 * t ** 0.63)
    phi2 = np.exp(-1.87 * t ** 1.22)
    return (1 - G_SLOPE) * phi1 + G_SLOPE * phi2


def observability(r_body, r_earth, h):
    """(V magnitude, solar elongation in degrees, geocentric unit vector)."""
    geo = r_body - r_earth
    delta = np.linalg.norm(geo, axis=0)
    r = np.linalg.norm(r_body, axis=0)
    cos_alpha = np.sum(r_body * geo, axis=0) / (r * delta)
    alpha = np.arccos(np.clip(cos_alpha, -1, 1))
    v = h + 5 * np.log10(r * delta) - 2.5 * np.log10(np.maximum(phase_function(alpha), 1e-30))
    rho_e = np.linalg.norm(r_earth, axis=0)
    cos_elong = np.sum(-r_earth * geo, axis=0) / (rho_e * delta)
    elong = np.degrees(np.arccos(np.clip(cos_elong, -1, 1)))
    return v, elong, geo / delta


def draws_from(orbit, n, rng):
    """Element draws (n, 6) from NEODyS's Gaussian; unbound or unphysical draws dropped."""
    # NEODyS prints nine significant digits, and the rounding leaves some
    # matrices with eigenvalues around -1e-15 of their largest. They are
    # clipped to zero; on the 69 bodies of 2026-10-29 the worst was -1.5e-11.
    cov = 0.5 * (orbit["cov"] + orbit["cov"].T)
    w, v = np.linalg.eigh(cov)
    cov = (v * np.clip(w, 0.0, None)) @ v.T
    x = rng.multivariate_normal(orbit["kep"], cov, size=n, method="eigh", check_valid="ignore")
    ok = (x[:, 0] > 0) & (x[:, 1] >= 0) & (x[:, 1] < 1)
    return x[ok]


def sigma_m_linear(orbit, years):
    """Linear 1-sigma of the mean anomaly (deg) after `years`: M = M0 + n(a) dt."""
    a = orbit["kep"][0]
    n = GAUSS_N_DEG_DAY * a ** -1.5
    dn_da = -1.5 * n / a
    dt = years * 365.25
    c = orbit["cov"]
    var = c[5, 5] + (dn_da * dt) ** 2 * c[0, 0] + 2 * dn_da * dt * c[0, 5]
    return float(np.sqrt(max(var, 0.0)))


def sigma_m_mc(draws, years):
    n = GAUSS_N_DEG_DAY * draws[:, 0] ** -1.5
    m = draws[:, 5] + n * years * 365.25
    centred = (m - np.median(m) + 180.0) % 360.0 - 180.0
    return float(np.std(centred))


def check_earth():
    lon = np.degrees(np.arctan2(*earth_state([MJD_J2000])[1::-1, 0])) % 360
    r = np.linalg.norm(earth_state(np.arange(MJD_J2000, MJD_J2000 + 366)), axis=0)
    ok = abs(lon - EARTH_LON_J2000_DEG) < 0.02 and 0.982 < r.min() and r.max() < 1.018
    print(f"EARTH  longitude at J2000 {lon:.3f} deg (expect {EARTH_LON_J2000_DEG}); distance"
          f" {r.min():.4f}-{r.max():.4f} AU over a year"
          f" ({'OK' if ok else 'MISMATCH -- the Earth model is wrong'})")
    return ok


def assess(orbit, h, start_mjd, end_mjd, n_draws, rng, step_days=2.0):
    """Window, observability share and search radius for one body."""
    days = np.arange(start_mjd, end_mjd, step_days)
    r_e = earth_state(days)
    nominal = body_positions(orbit["kep"], orbit["mjd"], days)[:, 0, :]
    v, elong, _ = observability(nominal, r_e, h)
    visible = (v <= V_LIMIT) & (elong >= MIN_ELONGATION_DEG)

    draws = draws_from(orbit, n_draws, rng)
    pos = body_positions(draws, orbit["mjd"], days)                       # (3, n, t)
    vd, ed, _ = observability(pos, r_e[:, None, :], h)
    p_obs = float(np.mean(np.any((vd <= V_LIMIT) & (ed >= MIN_ELONGATION_DEG), axis=1)))

    out = {"window": bool(visible.any()), "p_observable": p_obs, "search_deg": np.nan,
           "window_date": None, "v_best": float(np.min(v)),
           "sigma_m_10yr_mc": sigma_m_mc(draws_from(orbit, 4000, rng), 10.0),
           "sigma_m_10yr_lin": sigma_m_linear(orbit, 10.0)}
    if visible.any():
        k = int(np.argmin(np.where(visible, v, np.inf)))
        _, _, u_nom = observability(nominal[:, k:k + 1], r_e[:, k:k + 1], h)
        _, _, u_draw = observability(pos[:, :, k], r_e[:, k:k + 1], h)
        sep = np.degrees(np.arccos(np.clip(np.sum(u_draw * u_nom, axis=0), -1, 1)))
        out["search_deg"] = float(np.percentile(sep, 95))
        out["window_date"] = (datetime.date(1858, 11, 17)
                              + datetime.timedelta(days=float(days[k]))).isoformat()
    return out


def verdict(r):
    if not r["window"]:
        return "no window"
    if r["search_deg"] < 1:
        return "recoverable"
    if r["search_deg"] <= 10:
        return "search"
    return "lost"


def main():
    ap = argparse.ArgumentParser(description="recovery risk for the top of the ranking")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--start", default=None,
                    help="first day searched, YYYY-MM-DD (default: today)")
    args = ap.parse_args()

    if not check_earth():
        return
    start = datetime.date.fromisoformat(args.start) if args.start else datetime.date.today()

    cat = pd.read_csv(args.catalog, low_memory=False,
                      usecols=["designation", "is_neo", "orbit_condition_code",
                               "absolute_magnitude_h", "observation_arc_days"])
    cat["designation"] = cat["designation"].astype(str).str.strip()
    cat = cat[cat["is_neo"].astype(str).str.lower().isin(["true", "1"])]
    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"],
                       nrows=args.top)
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    top = prof.merge(cat, on="designation", how="inner")
    print(f"\ntop {args.top} rows: {len(top)} NEOs; searching {start} to {END}")
    print(f"  {snapshots.describe()}")

    rng = np.random.default_rng(args.seed)
    rows, misses, lin_bad = [], 0, []
    for _, b in top.iterrows():
        orbit, _why = fetch(b.designation)
        if not snapshots.replay_folder():
            time.sleep(PAUSE_S)
        if orbit is None:
            misses += 1
            continue
        h = b.absolute_magnitude_h
        r = assess(orbit, h, mjd(start), mjd(END), args.draws, rng)
        lin, mc = r["sigma_m_10yr_lin"], r["sigma_m_10yr_mc"]
        if abs(mc - lin) / lin > 0.05 and lin < 30:
            lin_bad.append((b.designation, lin, mc))
        rows.append(dict(rank=b["rank"], designation=b.designation, U=b.orbit_condition_code,
                         H=h, **r))
    res = pd.DataFrame(rows)
    res["verdict"] = res.apply(verdict, axis=1)

    print(f"  assessed {len(res)}; not in NEODyS {misses}")
    print(f"LINEAR vs MONTE CARLO  sigma_M after ten years: {len(res) - len(lin_bad)} of"
          f" {len(res)} agree to 5%"
          + (" (OK)" if not lin_bad else " -- disagreeing: "
             + ", ".join(f"{d} {a:.3g} vs {m:.3g}" for d, a, m in lin_bad[:5])))
    if lin_bad:
        print("  MISMATCH -- the covariance is not being propagated consistently")
        return

    good = res.dropna(subset=["U"])
    rho = good[["sigma_m_10yr_lin", "U"]].corr(method="spearman").iloc[0, 1]
    print(f"  10-year mean-anomaly drift against the catalog's U, Spearman: {rho:.3f}"
          f" (n = {len(good)})")

    print(f"\nRECOVERY BEFORE {END}  (V <= {V_LIMIT}, elongation >= {MIN_ELONGATION_DEG:.0f} deg)")
    for v in ("recoverable", "search", "lost", "no window"):
        s = res[res.verdict == v]
        print(f"  {v:12s} {len(s):4,}   median H {s.H.median() if len(s) else np.nan:5.1f}"
              f"   median p_observable {s.p_observable.median() if len(s) else np.nan:.2f}")
    print("\n  by condition code:")
    for lab, sel in (("U 0-4", res.U <= 4), ("U 5-6", res.U.between(5, 6)), ("U 7-9", res.U >= 7)):
        s = res[sel.fillna(False)]
        if len(s):
            print(f"    {lab}  n={len(s):3}  " + "  ".join(
                f"{v} {int((s.verdict == v).sum())}" for v in
                ("recoverable", "search", "lost", "no window")))
    print("\n  per body:")
    cols = ["rank", "designation", "U", "H", "v_best", "window_date", "search_deg",
            "p_observable", "sigma_m_10yr_lin", "verdict"]
    print(res[cols].to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
