# -*- coding: utf-8 -*-
"""Can a body that spins this fast hold the regolith the model digs?

economicspace mines "run-of-mine regolith at bulk grade" from every body it
ranks, and it reads `rotation_period_h` only to size the battery for the
night. A body that spins faster than about two hours cannot keep loose
material on its equator: centrifugal acceleration there beats its gravity. It
is a monolith, or held together by cohesion, and either way it has no loose
layer for a bucket-drum to scoop.

Two limits, both for a sphere of the body's own density:

    lift period     sqrt(3 pi / (G rho)). Faster than this, a loose grain on
                    the equator is not bound. Size-independent.
    cohesionless    Holsapple (2007): the fastest a strengthless rubble pile
                    with a 32 degree friction angle can spin before it yields
                    (Drucker-Prager). Faster than this it needs cohesion,
                    and `min_cohesion_pa` says how much.

The cohesion model is adapted from Karmanplus/prospector's
`prospector/enrichment/models/cohesion.py` (Apache-2.0), reduced to the
sphere, where Holsapple's three shape integrals have a closed form. That
reduction found a bug upstream. prospector's `_shape_integral` gives the
x-axis integrand the exponents (1.5, 0.5, 1.5), where Holsapple's are
(1.5, 0.5, 0.5), so for a sphere it returns A_x = 0.4 instead of 2/3 and the
three integrals sum to 1.733 instead of 2. This probe uses 2/3 for all three.
It checks that it reproduces upstream's figures exactly when A_x is forced to
0.4, which shows the difference is that one exponent and nothing else.

The top of the profitability ranking has no measured periods, so for each
such body the probe reports the share of NEOs within 0.75 mag of its H whose
measured period is faster than its lift period. The window and its width come
from prospector's `spin.py`; the empirical share replaces its Johnson S_B fit.

This script READS ONLY. It reads the catalog and the profitability catalog,
does not import master, does not build a stage, and writes nothing to
`asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_spin.py
    py ../starred-repos/probe_spin.py --top 40
"""

import argparse

import numpy as np
import pandas as pd

CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"

# From prospector/constants.py and cohesion.py (0.4.3).
G = 6.674e-11                        # m^3 kg^-1 s^-2
FRICTION_ANGLE_DEG = 32.0
SECONDS_PER_HOUR = 3600.0
# Holsapple's shape integrals for a sphere: each is the integral of
# (u + 1)^(-5/2) from 0 to infinity, which is 2/3, and the three sum to 2.
A_SPHERE = 2.0 / 3.0
# What prospector 0.4.3's `_shape_integral` returns for the x-axis of a sphere:
# the integral of (u + 1)^(-7/2), from its exponents (1.5, 0.5, 1.5).
A_X_UPSTREAM = 0.4
MAG_WINDOW = 0.75                    # prospector spin.py MAGNITUDE_WINDOW

# prospector 0.4.3's own `min_cohesion_pa` and `critical_period_h(k=0)`, run
# from its source on 2026-10-29: (radius m, period h, density kg/m^3,
# cohesion Pa, cohesionless critical period h).
UPSTREAM = [
    (0.5,    1 / 60, 2000.0, 1.176425195e+00, 3.495457719),
    (50.0,   0.5,    2000.0, 1.277704262e+01, 3.495457719),
    (150.0,  2.0,    2400.0, 5.177945298e+00, 3.190901736),
    (500.0,  3.0,    1300.0, 1.182154628e+01, 4.335581705),
    (1000.0, 8.0,    2000.0, 0.0,             3.495457719),
]


def friction_coefficient(phi_deg=FRICTION_ANGLE_DEG):
    """prospector's Drucker-Prager slope s for a friction angle."""
    phi = np.radians(phi_deg)
    return (2 * np.sin(phi)) / (np.sqrt(3) * (3 - np.sin(phi)))


def stresses(r, omega, rho, a_x=A_SPHERE):
    """Principal stresses at the centre of a spinning sphere (Pa).

    prospector's `_principal_stresses` with the shape integrals in closed
    form. `a_x` exists only so the upstream check can reproduce its x-axis.
    """
    gravity = 2 * np.pi * rho ** 2 * G
    sx = (rho * omega ** 2 - gravity * a_x) * r ** 2 / 5
    sy = (rho * omega ** 2 - gravity * A_SPHERE) * r ** 2 / 5
    sz = (-gravity * A_SPHERE) * r ** 2 / 5
    return sx, sy, sz


def min_cohesion_pa(r, period_h, rho, a_x=A_SPHERE):
    """Cohesion (Pa) a sphere needs to survive its spin; never negative. Vectorised."""
    omega = 2 * np.pi / (np.asarray(period_h, float) * SECONDS_PER_HOUR)
    sx, sy, sz = stresses(np.asarray(r, float), omega, np.asarray(rho, float), a_x)
    shear = ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) / 6
    return np.maximum(np.sqrt(shear) + friction_coefficient() * (sx + sy + sz), 0.0)


def cohesionless_period_h(rho, a_x=A_SPHERE):
    """The fastest period (h) a strengthless rubble pile survives. Size-independent.

    For a sphere with a_x = 2/3 the yield condition is linear in omega^2, so
    this is closed form: omega^2 = 4 pi s rho G / (1/sqrt(3) + 2 s). With
    upstream's a_x it is solved numerically, by bisection on min_cohesion_pa.
    """
    rho = np.asarray(rho, float)
    s = friction_coefficient()
    if a_x == A_SPHERE:
        omega2 = 4 * np.pi * s * rho * G / (1 / np.sqrt(3) + 2 * s)
        return 2 * np.pi / np.sqrt(omega2) / SECONDS_PER_HOUR
    lo, hi = np.full(rho.shape, 0.5), np.full(rho.shape, 20.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        needs = min_cohesion_pa(100.0, mid, rho, a_x) > 0
        lo, hi = np.where(needs, mid, lo), np.where(needs, hi, mid)
    return 0.5 * (lo + hi)


def lift_period_h(rho):
    """Faster than this, a loose grain on the equator of a sphere is unbound."""
    return np.sqrt(3 * np.pi / (G * np.asarray(rho, float))) / SECONDS_PER_HOUR


def check_upstream():
    """Upstream's figures with its a_x, and the closed form against the bisection."""
    print("UPSTREAM  (with prospector's a_x = 0.4, this must reproduce its own figures)")
    worst_k, worst_p = 0.0, 0.0
    for r, p, rho, k_up, pc_up in UPSTREAM:
        k = float(min_cohesion_pa(r, p, rho, A_X_UPSTREAM))
        pc = float(cohesionless_period_h(rho, A_X_UPSTREAM))
        worst_k = max(worst_k, abs(k - k_up) / max(k_up, 1e-9))
        worst_p = max(worst_p, abs(pc - pc_up))
        print(f"  r {r:7.1f} m  P {p:7.4f} h  rho {rho:5.0f}   cohesion {k:12.6e} (upstream"
              f" {k_up:12.6e})   cohesionless period {pc:.6f} h (upstream {pc_up:.6f})")
    ok = worst_k < 1e-8 and worst_p < 1e-6
    print(f"  worst relative cohesion gap {worst_k:.1e}, worst period gap {worst_p:.1e} h "
          f"({'OK' if ok else 'MISMATCH -- this is not prospector model'})")

    bisected = float(cohesionless_period_h(2000.0, A_SPHERE + 1e-15))
    closed = float(cohesionless_period_h(2000.0))
    flat_ok = abs(bisected - closed) < 1e-6
    print(f"  corrected model, rho 2000: closed form {closed:.6f} h, bisection {bisected:.6f} h"
          f" ({'OK' if flat_ok else 'MISMATCH -- the closed form is wrong'})")
    print(f"  lift period, rho 2000: {float(lift_period_h(2000.0)):.4f} h")
    return ok and flat_ok


def load(path):
    cols = ["designation", "is_neo", "absolute_magnitude_h", "diameter_km",
            "rotation_period_h", "comp_density_est_gcm3", "comp_group"]
    c = pd.read_csv(path, usecols=cols, low_memory=False)
    c = c[c["is_neo"].astype(str).str.lower().isin(["true", "1"])].copy()
    c["designation"] = c["designation"].astype(str).str.strip()
    c["rho"] = pd.to_numeric(c["comp_density_est_gcm3"], errors="coerce") * 1000.0
    c["r_m"] = pd.to_numeric(c["diameter_km"], errors="coerce") * 500.0
    return c


def main():
    ap = argparse.ArgumentParser(description="spin limits against the regolith assumption")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    if not check_upstream():
        return

    neo = load(args.catalog)
    m = neo.dropna(subset=["rotation_period_h", "rho", "r_m"]).copy()
    m = m[m.rotation_period_h > 0]
    m["p_lift"] = lift_period_h(m.rho)
    m["p_cohesionless"] = cohesionless_period_h(m.rho)
    m["k_pa"] = min_cohesion_pa(m.r_m, m.rotation_period_h, m.rho)
    m["k_pa_upstream"] = min_cohesion_pa(m.r_m, m.rotation_period_h, m.rho, A_X_UPSTREAM)
    faster_lift = m.rotation_period_h < m.p_lift
    faster_dp = m.rotation_period_h < m.p_cohesionless

    print(f"\nNEOs: {len(neo):,}   with a measured period, density and size: {len(m):,}")
    print(f"  lift period         median {m.p_lift.median():.3f} h   faster than it:"
          f" {int(faster_lift.sum()):,} ({faster_lift.mean() * 100:.1f}%)")
    print(f"  cohesionless period median {m.p_cohesionless.median():.3f} h   faster than it:"
          f" {int(faster_dp.sum()):,} ({faster_dp.mean() * 100:.1f}%)")
    print(f"  upstream's cohesionless period would call {int((m.rotation_period_h < cohesionless_period_h(m.rho.to_numpy(), A_X_UPSTREAM)).sum()):,}"
          " fast, against the corrected model's above")
    k = m.k_pa[faster_dp]
    print(f"  cohesion the fast ones need: median {k.median():.2f} Pa  p90 {k.quantile(0.9):.1f} Pa"
          f"  max {k.max():.0f} Pa;  over 100 Pa: {int((k > 100).sum()):,}")

    print("\n  by H (measured periods):")
    hb = pd.cut(m.absolute_magnitude_h, [0, 18, 20, 22, 25, 28, 40], right=False)
    for b, s in m.groupby(hb, observed=False):
        if len(s):
            fl = (s.rotation_period_h < s.p_lift).mean()
            print(f"    H {str(b):12s} n={len(s):5,}  median period {s.rotation_period_h.median():8.3f} h"
                  f"   faster than lift {fl * 100:5.1f}%")

    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"])
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    top = prof.head(1000).merge(neo, on="designation", how="inner")
    mags = m.absolute_magnitude_h.to_numpy(float)
    fast = faster_lift.to_numpy()
    share, n_win = [], []
    for h in top.absolute_magnitude_h.to_numpy(float):
        w = np.abs(mags - h) < MAG_WINDOW
        n_win.append(int(w.sum()))
        share.append(fast[w].mean() if w.sum() >= 10 else np.nan)
    top["n_window"] = n_win
    top["p_faster_than_lift"] = share
    top["measured_fast"] = top.rotation_period_h < lift_period_h(top.rho)

    print("\nPROFITABILITY RANKING  (NEOs among the top rows)")
    for n in (25, 100, 1000):
        t = top[top["rank"] <= n]
        known = t.p_faster_than_lift.notna()
        print(f"  top {n:5,}: {len(t):4,} NEOs   measured period {int(t.rotation_period_h.notna().sum()):3,}"
              f" (faster than lift {int(t.measured_fast.sum()):2,})"
              f"   size-class share faster than lift: median {t.p_faster_than_lift.median():.3f}"
              f" over {int(known.sum()):,}; unknown (window < 10) {int((~known).sum()):,}")
    print(f"\n  the top {args.top}:")
    cols = ["rank", "designation", "absolute_magnitude_h", "diameter_km",
            "rotation_period_h", "n_window", "p_faster_than_lift"]
    print(top[top["rank"] <= args.top][cols].to_string(
        index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
