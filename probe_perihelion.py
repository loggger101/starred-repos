# -*- coding: utf-8 -*-
"""How close to the Sun has a near-Earth body probably been, and does it matter?

Water is worth more in space than anything else this pipeline sells, and a
C-type NEO is credited with water from its taxonomy alone (FINDINGS.md F2:
90.2% of NEO taxonomies are inferred from an albedo, mostly an assumed one).
A body's present orbit is not its history. NEOs random-walk under planetary
encounters, and one that once passed inside about 0.2 AU has had its surface
heated and cracked enough to drive hydrated minerals off.

Toliou, Granvik and Tsirvoulis (2021) integrated a synthetic NEO population
and tabulated, on a grid of (a, e, i, H), the probability that a body's
minimum perihelion ever fell below each of 26 thresholds from 0.05 to 1.30 AU.
This probe reads that grid for every NEO in the catalog and for the top of the
profitability ranking, and reports two numbers per body:

    q_min_exp     the expected minimum perihelion, prospector's reduction of
                  the cumulative curve to its mean
    p_q_lt_0p2    P(q_min < 0.2 AU), straight off the table's fourth column

Both describe the POPULATION on this body's orbit today, not this body. Two
bodies on one orbit can have had different pasts, so read either as a prior.

The reader and `expected_q_min` are adapted from Karmanplus/prospector's
`prospector/enrichment/models/perihelion.py` (Apache-2.0); CITATIONS.md
section 3 carries its notice and what changed. Before measuring, the probe
reproduces prospector 0.4.3's own figures on five bodies and checks its fast
grid lookup against prospector's brute-force nearest cell on a sample.

The table is 248 MB gzipped and is downloaded once to a cache OUTSIDE this
checkout (default ~/.cache/starred-repos, since this repo lives in OneDrive),
checked against the sha256 recorded below, and parsed once to a .npy there.

This script READS ONLY. It reads the catalog and the profitability catalog,
does not import master, does not build a stage, and writes nothing to
`asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_perihelion.py
    py ../starred-repos/probe_perihelion.py --cache D:/somewhere --top 100
"""

import argparse
import gzip
import hashlib
import os
import urllib.request

import numpy as np
import pandas as pd

TOLIOU_URL = ("https://www.mv.helsinki.fi/home/mgranvik/data/"
              "Toliou+_2021_MNRAS/Toliou+_2021_MNRAS.dat.gz")
# The download of 2026-10-29: 247,823,013 bytes.
TOLIOU_SHA256 = "8b895c663b461054a86ecd43ee17dfc0b84e3278c96fd753a29bf19d68a7ad86"
DEFAULT_CACHE = os.path.join(os.path.expanduser("~"), ".cache", "starred-repos")
CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"

# From prospector/enrichment/models/perihelion.py: the ranges that make the
# four axes comparable in the nearest-cell distance, and the 26 thresholds of
# columns 5-30, P(q_min <= q_s) for q_s = 0.05 ... 1.30 AU.
AXIS_RANGES = np.array([4.2, 1.0, 180.0, 10.0])      # a (AU), e, i (deg), H (mag)
THRESHOLDS_AU = np.arange(0.05, 1.35, 0.05)
N_COLS = 30
COL_Q_LT_0P2 = 4 + 3                                  # the 0.20 AU column

# prospector 0.4.3's own expected_min_perihelion_au on these bodies, run from
# its source against the same download on 2026-10-29.
UPSTREAM = pd.DataFrame({
    "body":    ["Bennu", "Itokawa", "Ryugu", "Eros", "Phaethon"],
    "a":       [1.126, 1.324, 1.190, 1.458, 1.271],
    "e":       [0.204, 0.280, 0.190, 0.223, 0.890],
    "i":       [6.0,   1.6,   5.88,  10.8,  22.3],
    "H":       [20.2,  19.2,  19.3,  10.4,  14.3],
    "q_min_upstream": [0.569230, 0.564210, 0.576365, 0.738385, 0.047345],
})


def fetch(cache_dir):
    """The gzipped table, downloaded once and checked against TOLIOU_SHA256."""
    os.makedirs(cache_dir, exist_ok=True)
    gz = os.path.join(cache_dir, "Toliou+_2021_MNRAS.dat.gz")
    if not os.path.exists(gz):
        print(f"  downloading {TOLIOU_URL} (248 MB) to {cache_dir}")
        urllib.request.urlretrieve(TOLIOU_URL, gz + ".part")
        os.replace(gz + ".part", gz)
    h = hashlib.sha256()
    with open(gz, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    if h.hexdigest() != TOLIOU_SHA256:
        raise SystemExit(f"{gz}: sha256 {h.hexdigest()} is not the recorded "
                         f"{TOLIOU_SHA256}. The upstream file changed; re-record "
                         "it deliberately, since every figure in F12 rests on it.")
    return gz


def load_table(cache_dir):
    """The grid as an (n, 30) array: a, e, i, H, then 26 cumulative probabilities.

    Adapted from prospector's `load_toliou_table`. The download is a tar
    archive inside gzip, so the tar header and padding sit in the text; a line
    is kept only if it carries 30 numbers. Parsed once and saved beside the
    download, as prospector does, because parsing takes about a minute.
    """
    gz = fetch(cache_dir)
    npy = os.path.join(cache_dir, "toliou_2021_q.npy")
    if os.path.exists(npy) and os.path.getmtime(npy) >= os.path.getmtime(gz):
        return np.load(npy)
    rows = []
    with gzip.open(gz, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < N_COLS:
                continue
            try:
                rows.append([float(x) for x in parts[:N_COLS]])
            except ValueError:
                continue
    table = np.asarray(rows)
    np.save(npy, table)
    return table


class Grid:
    """Nearest grid cell, as prospector picks it, without a scan per body.

    Prospector takes the argmin of a normalised Euclidean distance over every
    row. That metric is a sum over axes, so on a Cartesian grid the nearest
    cell is the nearest value on each axis independently. The table is
    regular on all four axes (39 x 25 x 45 x 40) but only 26% filled, so a
    snapped cell that does not exist falls back to prospector's full scan.
    """

    def __init__(self, table):
        self.table = table
        self.axes = [np.unique(table[:, k]) for k in range(4)]
        shape = tuple(len(ax) for ax in self.axes)
        self.index = np.full(shape, -1, dtype=np.int64)
        idx = tuple(np.searchsorted(self.axes[k], table[:, k]) for k in range(4))
        self.index[idx] = np.arange(len(table))

    def brute(self, x):
        """Prospector's own rule: argmin of the normalised distance over all rows."""
        out = np.empty(len(x), dtype=np.int64)
        for j, row in enumerate(x):
            d = np.sum(((self.table[:, :4] - row) / AXIS_RANGES) ** 2, axis=1)
            out[j] = np.argmin(d)
        return out

    def lookup(self, a, e, i, h):
        """Row index of the nearest cell for each body. Returns (rows, n_fallback)."""
        x = np.column_stack([a, e, i, h]).astype(float)
        snapped = []
        for k, ax in enumerate(self.axes):
            j = np.clip(np.searchsorted(ax, x[:, k]), 1, len(ax) - 1)
            nearer_low = (x[:, k] - ax[j - 1]) <= (ax[j] - x[:, k])
            snapped.append(np.where(nearer_low, j - 1, j))
        rows = self.index[tuple(snapped)]
        missing = rows < 0
        if missing.any():
            rows[missing] = self.brute(x[missing])
        return rows, int(missing.sum())


def expected_q_min(cum):
    """Prospector's reduction: difference the cumulative curve, average the bin midpoints."""
    per_bin = np.diff(cum, prepend=0, axis=-1)
    return (per_bin * (THRESHOLDS_AU - 0.025)).sum(axis=-1)


def check_upstream(grid):
    """Reproduce prospector's figures, then the fast lookup against its full scan."""
    u = UPSTREAM.copy()
    rows, _ = grid.lookup(u.a, u.e, u.i, u.H)
    u["q_min_adapted"] = expected_q_min(grid.table[rows, 4:N_COLS])
    u["delta"] = u["q_min_adapted"] - u["q_min_upstream"]
    print("UPSTREAM  (must reproduce prospector 0.4.3's expected_min_perihelion_au)")
    print(u.to_string(index=False, float_format=lambda v: f"{v:9.6f}"))
    worst = u["delta"].abs().max()
    ok = worst < 1e-6
    print(f"\n  worst |adapted - upstream| = {worst:.1e} AU "
          f"({'OK' if ok else 'MISMATCH -- this is not prospector lookup'})")
    print("  Eros (H = 10.4) and Phaethon (H = 14.3) lie below the grid's H range,"
          " 15.125-24.875, and take its edge, as upstream does")
    return ok


def check_lookup(grid, neo, n=1500, seed=1):
    """The fast lookup must pick a cell exactly as near as prospector's full scan."""
    s = neo.sample(n=min(n, len(neo)), random_state=seed)
    x = s[["semi_major_axis_au", "eccentricity", "inclination_deg",
           "absolute_magnitude_h"]].to_numpy(float)
    fast, _ = grid.lookup(*x.T)
    slow = grid.brute(x)

    def dist(rows):
        return np.sum(((grid.table[rows, :4] - x) / AXIS_RANGES) ** 2, axis=1)

    gap = np.abs(dist(fast) - dist(slow)).max()
    ok = gap < 1e-12
    print(f"  fast lookup vs full scan on {len(s):,} catalog NEOs: worst distance gap "
          f"{gap:.1e} ({'OK' if ok else 'MISMATCH -- the snapped cell is not the nearest'})")
    return ok


def load_neos(path):
    """Every NEO with usable elements and H."""
    cols = ["designation", "is_neo", "semi_major_axis_au", "eccentricity",
            "inclination_deg", "perihelion_au", "absolute_magnitude_h",
            "spectral_type", "spectral_type_source", "comp_group", "comp_ice_fraction"]
    c = pd.read_csv(path, usecols=cols, low_memory=False)
    c = c[c["is_neo"].astype(str).str.lower().isin(["true", "1"])]
    c["designation"] = c["designation"].astype(str).str.strip()
    need = ["semi_major_axis_au", "eccentricity", "inclination_deg", "absolute_magnitude_h"]
    return c.dropna(subset=need).copy()


def annotate(grid, df):
    rows, n_fb = grid.lookup(df.semi_major_axis_au, df.eccentricity,
                             df.inclination_deg, df.absolute_magnitude_h)
    cum = grid.table[rows, 4:N_COLS]
    df["q_min_exp"] = expected_q_min(cum)
    df["p_q_lt_0p2"] = grid.table[rows, COL_Q_LT_0P2]
    return n_fb


def describe(df, label):
    p = df["p_q_lt_0p2"]
    print(f"  {label:28s} n={len(df):7,}  q_min_exp median {df.q_min_exp.median():.3f} AU"
          f"   P(q_min<0.2) median {p.median():.3f}  mean {p.mean():.3f}"
          f"   > 0.5 on {(p > 0.5).mean() * 100:5.1f}%")


def main():
    """Checks, then every NEO, then the top of the profitability ranking."""
    ap = argparse.ArgumentParser(description="Toliou et al. (2021) perihelion prior")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    ap.add_argument("--top", type=int, default=25, help="rows of the ranking to list")
    args = ap.parse_args()

    table = load_table(args.cache)
    grid = Grid(table)
    print(f"Toliou grid: {len(table):,} cells, axes "
          + " x ".join(str(len(ax)) for ax in grid.axes) + "\n")
    if not check_upstream(grid):
        return
    neo = load_neos(args.catalog)
    if not check_lookup(grid, neo):
        return

    n_fb = annotate(grid, neo)
    h = neo.absolute_magnitude_h
    above_q = (neo.q_min_exp > neo.perihelion_au + 0.05).mean() * 100
    print(f"\nNEOs: {len(neo):,}   empty-cell fallbacks {n_fb:,}"
          f"   H outside the grid {int(((h < 15.0) | (h > 25.0)).sum()):,}")
    print(f"  expected q_min more than 0.05 AU above today's perihelion: {above_q:.1f}%"
          "  (a minimum cannot exceed the present value; this is grid coarseness)")

    print("\nBY COMPOSITION GROUP")
    for g, s in neo.groupby("comp_group"):
        if len(s) >= 200:
            describe(s, g)
    ice = neo[pd.to_numeric(neo.comp_ice_fraction, errors="coerce") > 0]
    describe(ice, "credited with ice")
    describe(ice[ice.spectral_type_source.ne("source")], "  of which taxonomy inferred")

    prof = pd.read_csv(args.profitability, low_memory=False,
                       usecols=["designation", "profit_usd", "cargo_water_kg",
                                "water_liberated_kg", "payload_dominant_phase"])
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    top = prof.merge(neo, on="designation", how="left")
    print(f"\nPROFITABILITY RANKING  ({len(prof):,} rows, in file order, which is by profit)")
    for n in (25, 100, 1000):
        t = top.head(n)
        joined = t["q_min_exp"].notna()
        print(f"  top {n:5,}: NEOs {int(joined.sum()):5,}"
              f"   C-complex {int(t.comp_group.eq('C-complex').sum()):5,}"
              f"   carrying water {int((pd.to_numeric(t.cargo_water_kg, errors='coerce') > 0).sum()):5,}")
        describe(t[joined], f"top {n:,} NEOs")
        wet = t[joined & (pd.to_numeric(t.cargo_water_kg, errors="coerce") > 0)]
        if len(wet):
            describe(wet, f"top {n:,} carrying water")

    print(f"\n  the top {args.top}:")
    cols = ["rank", "designation", "spectral_type", "spectral_type_source",
            "payload_dominant_phase", "perihelion_au", "q_min_exp", "p_q_lt_0p2"]
    print(top.head(args.top)[cols].to_string(
        index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
