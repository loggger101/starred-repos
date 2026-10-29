# -*- coding: utf-8 -*-
"""A debiased prior on whether a NEO is dark, from where in the main belt it came.

F11 measured the builder's no-albedo rule against NEOs with a measured class
and found a guessed C right 29% of the time. It also had to admit that bodies
with a measured class are biased toward bright S-types, so its prior is too
generous to S. A debiased prior needs a model of the whole NEO population,
corrected for what surveys can see.

Granvik et al. (2018) built that model. Their low-resolution model file gives,
for each cell of (a, e, i, H), the probability that a NEO there came from each
of seven source regions (the nu6 resonance, the 3:1, 5:2 and 2:1 resonances
with Jupiter, the Hungaria and Phocaea groups, and the Jupiter-family comets),
and the debiased probability that its geometric albedo is at most 0.1, or
between 0.1 and 0.3. An albedo at most 0.1 is the builder's own line between C
and S, so the model gives every NEO a survey-corrected P(the builder's split
would call it C).

This probe reads that file and asks:

    calibration  among NEOs with a MEASURED albedo, does the model's
                 P(albedo <= 0.1) match how often they really are that dark
    the guess    for the 38,616 NEOs typed from an assumed albedo, how often
                 the model disagrees with the class the builder gave, and how
                 many of them the model expects to be dark in total
    the top      P(dark) and the most likely source region for the NEOs at the
                 top of the profitability ranking

The model's grid is laid out like Toliou et al.'s (F12): the same axes and
steps, partly filled. The nearest-cell lookup is probe_perihelion.Grid,
checked against prospector's full scan before measuring, and each cell's seven
source probabilities must sum to 1.

The file is downloaded once (33 MB) to ~/.cache/starred-repos, checked
against the sha256 below, parsed once to a .npy there, and never committed.

This script READS ONLY. It reads the catalog and the profitability catalog
and writes nothing to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_source_regions.py
"""

import argparse
import gzip
import hashlib
import os
import urllib.request

import numpy as np
import pandas as pd

from probe_perihelion import Grid, check_lookup
from probe_taxonomy_prior import CLASSES, collapse

URL = ("https://www.mv.helsinki.fi/home/mgranvik/data/Granvik+_2018_Icarus/"
       "Granvik+_2018_Icarus_model_lowres_v20181213_albedo_cp.out.gz")
NAME = "Granvik+_2018_Icarus_model_lowres_v20181213_albedo_cp.out.gz"
# The download of 2026-10-29: 33,342,324 bytes.
SHA256 = "81e26f962f2573edb098ffae2d7bf254588ba2e6382bd3157e9ad10e057de42d"
DEFAULT_CACHE = os.path.join(os.path.expanduser("~"), ".cache", "starred-repos")
CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"
N_CELLS = 773280                       # the file's own header
SOURCES = ["nu6", "5/2", "2/1", "HUN", "3/1", "PHO", "JFC"]
# Columns: a, e, i, H, N, sigma N, seven source probabilities, their seven
# sigmas, P(pv <= 0.1), P(0.1 < pv <= 0.3), their two sigmas, collision
# probability, MOID.
COL_SRC = slice(6, 13)
COL_DARK, COL_MID = 20, 21
H_LO, H_HI = 15.0, 25.0                # the grid's H runs 15.125 to 24.875


def load_model(cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    gz = os.path.join(cache_dir, NAME)
    if not os.path.exists(gz):
        print(f"  downloading {URL} (33 MB) to {cache_dir}")
        urllib.request.urlretrieve(URL, gz + ".part")
        os.replace(gz + ".part", gz)
    with open(gz, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    if sha != SHA256:
        raise SystemExit(f"{gz}: sha256 {sha} is not the recorded {SHA256}; re-record it deliberately.")
    npy = os.path.join(cache_dir, "granvik_2018_lowres.npy")
    if os.path.exists(npy) and os.path.getmtime(npy) >= os.path.getmtime(gz):
        return np.load(npy)
    rows = []
    with gzip.open(gz, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if line.startswith("!") or len(parts) != 26:
                continue
            rows.append([float(x) for x in parts])
    table = np.asarray(rows)
    np.save(npy, table)
    return table


def check_model(table):
    """Row count against the header, and each cell's source probabilities summing to 1."""
    sums = table[:, COL_SRC].sum(axis=1)
    alb = table[:, COL_DARK] + table[:, COL_MID]
    ok_n = len(table) == N_CELLS
    ok_sum = bool(np.all(np.abs(sums - 1) < 2e-3))
    ok_alb = bool(np.all((alb >= -1e-6) & (alb <= 1 + 1e-6)))
    print(f"MODEL  {len(table):,} cells (header says {N_CELLS:,}); source probabilities sum to 1"
          f" within {np.abs(sums - 1).max():.1e}; albedo categories within [0, 1]: {ok_alb}"
          f" ({'OK' if ok_n and ok_sum and ok_alb else 'MISMATCH -- the file is misread'})")
    return ok_n and ok_sum and ok_alb


def annotate(grid, df):
    rows, n_fb = grid.lookup(df.semi_major_axis_au, df.eccentricity, df.inclination_deg,
                             df.absolute_magnitude_h)
    t = grid.table[rows]
    df["p_dark"] = t[:, COL_DARK]
    df["p_bright"] = np.clip(1 - t[:, COL_DARK] - t[:, COL_MID], 0, 1)
    src = t[:, COL_SRC]
    df["source"] = np.asarray(SOURCES)[src.argmax(axis=1)]
    df["p_source"] = src.max(axis=1)
    df["p_jfc"] = src[:, SOURCES.index("JFC")]
    return n_fb


def load(path):
    cols = ["designation", "is_neo", "semi_major_axis_au", "eccentricity", "inclination_deg",
            "absolute_magnitude_h", "albedo", "spectral_type_source", "comp_group"]
    c = pd.read_csv(path, usecols=cols, low_memory=False)
    c = c[c["is_neo"].astype(str).str.lower().isin(["true", "1"])].copy()
    c["designation"] = c["designation"].astype(str).str.strip()
    return c.dropna(subset=["semi_major_axis_au", "eccentricity", "inclination_deg",
                            "absolute_magnitude_h"])


def main():
    ap = argparse.ArgumentParser(description="Granvik et al. (2018) debiased albedo prior")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    table = load_model(args.cache)
    if not check_model(table):
        return
    grid = Grid(table)
    neo = load(args.catalog)
    if not check_lookup(grid, neo):
        return
    n_fb = annotate(grid, neo)
    in_h = neo.absolute_magnitude_h.between(H_LO, H_HI)
    print(f"\nNEOs {len(neo):,}; empty-cell fallbacks {n_fb:,}; H inside the grid {int(in_h.sum()):,}")

    m = neo[neo.albedo.between(0, 1, inclusive="neither") & in_h].copy()
    m["dark"] = m.albedo <= 0.1
    print(f"\nCALIBRATION  NEOs with a measured albedo and H inside the grid: {len(m):,}")
    print(f"  dark (albedo <= 0.1): observed {m.dark.mean():.3f}, model mean P {m.p_dark.mean():.3f}")
    for lo, hi in ((0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 1.01)):
        s = m[(m.p_dark >= lo) & (m.p_dark < hi)]
        if len(s) >= 20:
            print(f"    model P {lo:.1f}-{min(hi, 1):.1f}  n={len(s):5,}  mean P {s.p_dark.mean():.3f}"
                  f"  observed dark {s.dark.mean():.3f}")
    brier = np.mean((m.p_dark - m.dark) ** 2)
    base = np.mean((m.dark.mean() - m.dark) ** 2)
    print(f"  Brier score {brier:.4f} (a constant at the observed rate scores {base:.4f})")

    g = neo[neo.spectral_type_source.eq("albedo_assumed")].copy()
    g["assigned"] = collapse(g.comp_group)
    gi = g.absolute_magnitude_h.between(H_LO, H_HI)
    print(f"\nTHE GUESS  NEOs typed from an assumed albedo {len(g):,} (H inside the grid {int(gi.sum()):,})")
    for c, s in g.groupby("assigned"):
        print(f"  assigned {c:10s} n={len(s):6,}  model P(dark) median {s.p_dark.median():.3f}"
              f"  mean {s.p_dark.mean():.3f};  model disagrees (P(dark) on the other side of 0.5)"
              f" on {int(((s.p_dark > 0.5) if c != 'C-complex' else (s.p_dark < 0.5)).sum()):,}")
    print(f"  the builder calls {int((g.assigned == 'C-complex').sum()):,} of them dark (C);"
          f" the model expects {g.p_dark.sum():,.0f} ({g.p_dark.mean() * 100:.1f}%)")
    print(f"  inside the grid's H only: builder {int((g[gi].assigned == 'C-complex').sum()):,},"
          f" model {g[gi].p_dark.sum():,.0f}")

    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"], nrows=1000)
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    top = prof.merge(neo, on="designation", how="inner")
    top["assigned"] = collapse(top.comp_group)
    print("\nPROFITABILITY RANKING  (NEOs among the top rows)")
    for n in (25, 100, 1000):
        t = top[top["rank"] <= n]
        c = t[t.assigned == "C-complex"]
        print(f"  top {n:5,}: {len(t):4} NEOs; H inside the grid {int(t.absolute_magnitude_h.between(H_LO, H_HI).sum()):4};"
              f" typed C {len(c):4} with mean model P(dark) {c.p_dark.mean():.3f};"
              f" sources {t.source.value_counts().to_dict()}")
    print(f"\n  the top {args.top}:")
    cols = ["rank", "designation", "comp_group", "spectral_type_source", "absolute_magnitude_h",
            "p_dark", "p_bright", "source", "p_source"]
    print(top[top["rank"] <= args.top][cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
