# -*- coding: utf-8 -*-
"""How well is delta-v known for the bodies at the top of the ranking?

SPACE-MAP.md found the ranking blind to orbit quality, with poorly determined
orbits enriched 3.1x at its top, and left the response open: a cutoff, an arc
floor, or a confidence weighting. SECOND-PASS.md F9 found astroquery's NEODyS
client, which returns each NEO's full orbit covariance, and proposed pushing
that covariance through the delta-v estimator to get an error bar instead of a
filter. This probe does that for the NEOs at the top of the profitability
ranking.

For each body it fetches NEODyS's Keplerian elements at an epoch near the
present with their 6x6 covariance, draws from the Gaussian on (a, e, i), and
prices every draw with the shipped estimator (`dv_out_shipped` from
probe_nhats.py, both apsides). It reports the 5th to 95th percentile of
`dv_out` per body, against the orbit condition code U, and against how far
apart the bodies' nominal delta-v figures are.

Two checks run first. The shipped copy must reproduce the five docstring
figures (from probe_nhats.py). And for every body fetched, the square roots
of the parsed covariance's diagonal must equal the RMS line NEODyS prints
beside it, which shows the 21 upper-triangle values were placed correctly.

Nothing is copied from astroquery (BSD-3). The endpoint and the line tags
(KEP, MJD, COV) come from reading its `solarsystem/neodys/core.py`; the
parsing here is written against the file format.

A Gaussian on Keplerian elements is a linearisation. For an orbit fitted to a
day or two of observations the real uncertainty region is curved and far from
Gaussian, so for those bodies the interval says the orbit is poorly known and
roughly how poorly, and no more.

This script READS ONLY. It fetches one public read-only service, one request
per body with a pause between, does not import master, does not build a stage,
and writes nothing to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_covariance.py
    py ../starred-repos/probe_covariance.py --top 200 --draws 4000
"""

import argparse
import time
import urllib.error

import numpy as np
import pandas as pd

import snapshots
from probe_nhats import check_reference, dv_out_shipped

NEODYS_URL = "https://newton.spacedys.com/~neodys2/epoch/{}.ke1"
CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"
PAUSE_S = 0.3


def packed(designation):
    """NEODyS's file name: a number as is, a provisional designation without its space."""
    return str(designation).strip().replace(" ", "")


def parse(text):
    """KEP (6), the epoch, the 6x6 covariance and NEODyS's RMS line, or None."""
    kep, mjd, cov, rms = None, None, [], None
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "KEP":
            kep = np.array([float(x) for x in parts[1:7]])
        elif parts[0] == "MJD":
            mjd = float(parts[1])
        elif parts[0] == "COV":
            cov.extend(float(x) for x in parts[1:])
        elif parts[:2] == ["!", "RMS"]:
            rms = np.array([float(x) for x in parts[2:8]])
    if kep is None or len(cov) != 21:
        return None
    full = np.zeros((6, 6))
    full[np.triu_indices(6)] = cov
    full = full + full.T - np.diag(np.diag(full))
    return {"kep": kep, "mjd": mjd, "cov": full, "rms": rms}


def fetch(designation):
    """One body's NEODyS file, or (None, reason)."""
    url = NEODYS_URL.format(packed(designation))
    try:
        return parse(snapshots.fetch(url, timeout=60).decode("utf-8", errors="replace")), None
    except urllib.error.HTTPError as err:
        return None, f"HTTP {err.code}"
    except (urllib.error.URLError, TimeoutError) as err:
        return None, type(err).__name__


def rms_gap(orbit):
    """Worst relative gap between sqrt(diag(COV)) and the RMS line."""
    if orbit["rms"] is None:
        return np.nan
    sd = np.sqrt(np.diag(orbit["cov"]))
    return float(np.max(np.abs(sd - orbit["rms"]) / np.maximum(orbit["rms"], 1e-300)))


def dv_interval(orbit, draws, rng):
    """Shipped dv_out at the nominal elements and over Gaussian draws on (a, e, i)."""
    mean = orbit["kep"][:3]
    cov = orbit["cov"][:3, :3]
    x = rng.multivariate_normal(mean, cov, size=draws, method="eigh")
    ok = (x[:, 0] > 0) & (x[:, 1] >= 0) & (x[:, 1] < 1) & (x[:, 2] >= 0) & (x[:, 2] <= 180)
    dv, _ = dv_out_shipped(x[ok, 0], x[ok, 1], x[ok, 2])
    nominal, _ = dv_out_shipped(mean[:1], mean[1:2], mean[2:3])
    return float(nominal[0]), np.percentile(dv, [5, 50, 95]), 1 - ok.mean()


def main():
    ap = argparse.ArgumentParser(description="NEODyS covariance through the delta-v estimator")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--top", type=int, default=100, help="ranking rows to take NEOs from")
    ap.add_argument("--draws", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    if not check_reference():
        return

    cat = pd.read_csv(args.catalog, low_memory=False,
                      usecols=["designation", "is_neo", "orbit_condition_code",
                               "observation_arc_days", "absolute_magnitude_h",
                               "semi_major_axis_au", "eccentricity", "inclination_deg"])
    cat["designation"] = cat["designation"].astype(str).str.strip()
    cat = cat[cat["is_neo"].astype(str).str.lower().isin(["true", "1"])]
    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"],
                       nrows=args.top)
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    top = prof.merge(cat, on="designation", how="inner")
    print(f"\ntop {args.top} rows of the ranking: {len(top)} NEOs; fetching NEODyS")
    print(f"  {snapshots.describe()}")

    rng = np.random.default_rng(args.seed)
    rows, misses, worst_rms = [], {}, 0.0
    for _, b in top.iterrows():
        orbit, why = fetch(b.designation)
        if not snapshots.replay_folder():
            time.sleep(PAUSE_S)
        if orbit is None:
            misses[why or "unparsed"] = misses.get(why or "unparsed", 0) + 1
            continue
        gap = rms_gap(orbit)
        worst_rms = max(worst_rms, gap if np.isfinite(gap) else 0.0)
        nominal, (p05, p50, p95), rejected = dv_interval(orbit, args.draws, rng)
        rows.append({"rank": b["rank"], "designation": b.designation,
                     "U": b.orbit_condition_code, "arc_d": b.observation_arc_days,
                     "H": b.absolute_magnitude_h, "dv_nominal": nominal,
                     "dv_p05": p05, "dv_p50": p50, "dv_p95": p95,
                     "width90": p95 - p05, "rejected": rejected,
                     "sd_a": float(np.sqrt(orbit["cov"][0, 0])),
                     "rms_gap": gap})
    r = pd.DataFrame(rows)
    print(f"  fetched {len(r)}; missing {sum(misses.values())}"
          + (f" ({', '.join(f'{k}: {v}' for k, v in misses.items())})" if misses else ""))
    ok = worst_rms < 1e-3
    print(f"  worst |sqrt(diag COV) - RMS| / RMS over all fetched: {worst_rms:.1e} "
          f"({'OK' if ok else 'MISMATCH -- the covariance is misparsed'})")
    if not ok or r.empty:
        return

    print(f"\nDELTA-V UNCERTAINTY  (shipped dv_out, 5th-95th percentile, {args.draws:,} draws each)")
    print(f"  90% interval width: median {r.width90.median():.3f} km/s"
          f"   p90 {r.width90.quantile(0.9):.3f}   max {r.width90.max():.3f}")
    for lim in (0.1, 0.5, 1.0):
        print(f"  wider than {lim:3.1f} km/s: {int((r.width90 > lim).sum()):3,} of {len(r)}")
    print(f"  draws rejected as unbound or unphysical: median {r.rejected.median() * 100:.2f}%"
          f"   max {r.rejected.max() * 100:.1f}%")
    q1, q3 = r.dv_nominal.quantile([0.25, 0.75])
    print(f"  spread of nominal dv across these bodies: IQR {q3 - q1:.3f} km/s"
          f" ({q1:.3f} to {q3:.3f})")

    print("\n  by orbit condition code:")
    for lab, sel in (("U 0-2", r.U <= 2), ("U 3-4", r.U.between(3, 4)),
                     ("U 5-6", r.U.between(5, 6)), ("U 7-9", r.U >= 7)):
        s = r[sel.fillna(False)]
        if len(s):
            print(f"    {lab}  n={len(s):3}  width90 median {s.width90.median():8.3f} km/s"
                  f"   max {s.width90.max():8.3f}")

    print("\n  per body:")
    cols = ["rank", "designation", "U", "arc_d", "H", "dv_nominal", "dv_p05", "dv_p95",
            "width90", "rejected"]
    print(r[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
