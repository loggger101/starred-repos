# -*- coding: utf-8 -*-
"""More measured spins for F13, from the Asteroid Lightcurve Database.

F13 found that a third of the NEOs with a measured rotation period spin faster
than loose regolith can stay on, and that past H 28 all 51 measured bodies do.
Its periods came from the catalog, which has one for 2,440 NEOs. The bodies
that lead the profitability ranking sit near H 32, beyond any of them.

The LCDB (Warner, Harris and Pravec; minplanobs.org) is the standard
compilation of asteroid rotation periods, with a reliability code U for each.
space-datasets' `update-lcdb.py` (F8) is where it was found. This probe reads
its summary table, joins it to the catalog's NEOs, and asks four things:

    coverage     how many NEOs have a reliable LCDB period (U >= 2, no
                 upper- or lower-limit flag), against the catalog's 2,440
    agreement    where both have one, do they agree
    spin         F13's table again with LCDB periods: the share faster than
                 the lift period by H, and how far into faint H it now reaches
    albedo       how many NEOs typed from an assumed albedo have an albedo in
                 the LCDB that the LCDB does not itself flag as assumed (F16)

The fixed-width layout is read from the table's own header line, not copied
from anyone's column specification. Values the LCDB marks with source code A
(assumed from class or family) are not used as measurements.

The table is downloaded once (40 MB) to ~/.cache/starred-repos, checked
against the sha256 below, and never committed. Lift periods come from
probe_spin.py, so the physics is F13's.

This script READS ONLY. It reads the catalog and the profitability catalog
and writes nothing to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_lightcurves.py
"""

import argparse
import hashlib
import io
import os
import re
import urllib.request
import zipfile

import numpy as np
import pandas as pd

from probe_spin import lift_period_h
from probe_taxonomy_prior import CLASSES, collapse, rule_class

LCDB_URL = "https://minplanobs.org/MPInfo/datazips/LCLIST_PUB_CURRENT.zip"
# The download of 2026-10-29 (39,966,552 bytes), "GENERATED: Oct 01, 2023".
LCDB_SHA256 = "1a200680375897fbd62143572b2078137a0d42b170eabb764b1b1eef620640e2"
DEFAULT_CACHE = os.path.join(os.path.expanduser("~"), ".cache", "starred-repos")
CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"

# Rows the parse must reproduce, read off the table by eye on 2026-10-29:
# (designation as the catalog keys it, period h, H).
KNOWN = [("1", 9.07417, 3.53), ("25143", 12.132, 19.48), ("101955", 4.2968, 20.72),
         ("2008 TC3", 0.0269409, 30.9)]
H_BINS = [0, 18, 20, 22, 25, 28, 30, 40]


def fetch(cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, "LCLIST_PUB_CURRENT.zip")
    if not os.path.exists(path):
        print(f"  downloading {LCDB_URL} (40 MB) to {cache_dir}")
        urllib.request.urlretrieve(LCDB_URL, path + ".part")
        os.replace(path + ".part", path)
    with open(path, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    if sha != LCDB_SHA256:
        raise SystemExit(f"{path}: sha256 {sha} is not the recorded {LCDB_SHA256}. The LCDB "
                         "was updated; re-record it deliberately, since F18 rests on it.")
    return path


def header_columns(header, dashes):
    """Column (name, start, end) from the header line, one field per header token.

    Each field runs from its token's start to the next token's start. Repeated
    names (the LCDB reuses S and F for each field's source and flag) take the
    name of the field they precede, so S before CLASS becomes S_CLASS.
    """
    tokens = [(m.group(), m.start()) for m in re.finditer(r"\S+", header)]
    names = []
    for k, (tok, _) in enumerate(tokens):
        if tok in ("S", "F", "P") and k + 1 < len(tokens):
            nxt = next(t for t, _ in tokens[k + 1:] if t not in ("S", "F", "P"))
            names.append(f"{tok}_{nxt}")
        else:
            names.append(tok)
    starts = [s for _, s in tokens] + [len(dashes)]
    return [(n, starts[k], starts[k + 1]) for k, n in enumerate(names)]


def parse_summary(text):
    lines = text.splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith("NUMBER"))
    cols = header_columns(lines[k], lines[k + 1])
    body = "\n".join(lines[k + 2:])
    d = pd.read_fwf(io.StringIO(body), colspecs=[(s, e) for _, s, e in cols],
                    names=[n for n, _, _ in cols], dtype=str)
    d = d.apply(lambda s: s.str.strip())
    num = pd.to_numeric(d["NUMBER"], errors="coerce")
    d["designation"] = np.where(num > 0, num.astype("Int64").astype(str), d["DESIG"])
    for c in ("PERIOD", "H", "ALBEDO", "DIA.", "AMAX"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["U_num"] = pd.to_numeric(d["U"].str[0], errors="coerce")
    return d


def load_lcdb(cache_dir):
    with zipfile.ZipFile(fetch(cache_dir)) as z:
        return parse_summary(z.read("lc_summary_pub.txt").decode("latin-1"))


def check_parse(lc):
    ok = True
    print("PARSE  (known rows of the summary table)")
    for des, period, h in KNOWN:
        row = lc[lc.designation == des]
        got = (float(row.PERIOD.iloc[0]), float(row.H.iloc[0])) if len(row) else (np.nan, np.nan)
        good = len(row) == 1 and abs(got[0] - period) < 1e-9 and abs(got[1] - h) < 1e-9
        ok &= good
        print(f"  {des:10s} period {got[0]!s:>12} h (expect {period}), H {got[1]} (expect {h})"
              f"  {'OK' if good else 'MISMATCH'}")
    if not ok:
        print("  MISMATCH -- the fixed-width layout is misread")
    return ok


def reliable(lc):
    """Periods with U >= 2 and no limit flag (> or <) on the period."""
    flag = lc.get("F_PERIOD", pd.Series("", index=lc.index)).fillna("")
    return lc[(lc.U_num >= 2) & lc.PERIOD.gt(0) & ~flag.isin([">", "<"])]


def load_catalog(path):
    cols = ["designation", "provisional_designation", "is_neo", "absolute_magnitude_h",
            "rotation_period_h", "comp_density_est_gcm3", "spectral_type_source",
            "comp_group"]
    c = pd.read_csv(path, usecols=cols, low_memory=False)
    c = c[c["is_neo"].astype(str).str.lower().isin(["true", "1"])].copy()
    for k in ("designation", "provisional_designation"):
        c[k] = c[k].astype(str).str.strip()
    c["rho"] = pd.to_numeric(c["comp_density_est_gcm3"], errors="coerce") * 1000.0
    return c


def join(neo, lc, columns):
    """LCDB rows joined to the catalog on designation, then on provisional designation."""
    a = neo.merge(lc[["designation"] + columns], on="designation", how="inner")
    rest = neo[~neo.designation.isin(a.designation)]
    b = (rest.merge(lc[["designation"] + columns].rename(columns={"designation": "provisional_designation"}),
                    on="provisional_designation", how="inner"))
    return pd.concat([a, b], ignore_index=True).drop_duplicates("designation")


def main():
    ap = argparse.ArgumentParser(description="LCDB periods against F13")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    args = ap.parse_args()

    lc = load_lcdb(args.cache)
    print(f"LCDB summary: {len(lc):,} rows\n")
    if not check_parse(lc):
        return
    neo = load_catalog(args.catalog)
    rel = reliable(lc)
    j = join(neo, rel, ["PERIOD", "U", "H"]).rename(columns={"H": "H_lcdb"})
    j = j.dropna(subset=["rho"])
    print(f"\nCOVERAGE  NEOs {len(neo):,}; with a catalog period {int(neo.rotation_period_h.notna().sum()):,};"
          f" with a reliable LCDB period (U >= 2, no limit flag) {len(j):,}")
    new = j[j.rotation_period_h.isna()]
    print(f"  LCDB periods for NEOs the catalog has none for: {len(new):,}")

    # Which LCDB quality the catalog's own periods carry. F13 used every one.
    allp = join(neo[neo.rotation_period_h.notna()], lc, ["PERIOD", "U", "U_num", "F_PERIOD"])
    limit = allp.F_PERIOD.fillna("").isin([">", "<"])
    print(f"  the catalog's {int(neo.rotation_period_h.notna().sum()):,} NEO periods in the LCDB:"
          f" {len(allp):,}; by quality: U 3 {int(((allp.U_num == 3) & ~limit).sum()):,},"
          f" U 2 {int(((allp.U_num == 2) & ~limit).sum()):,}, U 1 {int(((allp.U_num == 1) & ~limit).sum()):,},"
          f" a limit (> or <) {int(limit.sum()):,}, no U {int(allp.U_num.isna().sum()):,}")

    both = j.dropna(subset=["rotation_period_h"])
    ratio = both.PERIOD / both.rotation_period_h
    print(f"\nAGREEMENT  where both have a period (n = {len(both):,}): within 1%"
          f" {np.mean(np.abs(ratio - 1) < 0.01) * 100:.1f}%, within 5% {np.mean(np.abs(ratio - 1) < 0.05) * 100:.1f}%;"
          f" a factor 2 apart {np.mean((np.abs(ratio - 2) < 0.05) | (np.abs(ratio - 0.5) < 0.025)) * 100:.1f}%")

    j["p_lift"] = lift_period_h(j.rho)
    j["fast"] = j.PERIOD < j.p_lift
    print("\nSPIN  (LCDB periods; faster than the lift period for the body's catalog density)")
    hb = pd.cut(j.absolute_magnitude_h, H_BINS, right=False)
    for b, s in j.groupby(hb, observed=False):
        if len(s):
            print(f"  H {str(b):10s} n={len(s):5,}  median period {s.PERIOD.median():9.4f} h"
                  f"   faster than lift {s.fast.mean() * 100:5.1f}%   (catalog had {int(s.rotation_period_h.notna().sum()):4,})")
    print(f"  all     n={len(j):5,}  faster than lift {j.fast.mean() * 100:.1f}%")
    faint = j[j.absolute_magnitude_h >= 30]
    if len(faint):
        print(f"  faintest: H {faint.absolute_magnitude_h.max():.2f}; at H >= 30 the slowest period is"
              f" {faint.PERIOD.max() * 60:.1f} min")

    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"], nrows=1000)
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    t = prof.merge(j[["designation", "PERIOD", "fast"]], on="designation", how="inner")
    print(f"\n  top 1,000 rows with a reliable LCDB period: {len(t)}"
          f" (faster than lift {int(t.fast.sum())}); in the top 25: {int((t['rank'] <= 25).sum())}")

    alb = lc[lc.ALBEDO.gt(0) & lc.S_ALBEDO.fillna("A").ne("A")]
    guessed = neo[neo.spectral_type_source.eq("albedo_assumed")]
    ga = join(guessed, alb, ["ALBEDO", "S_ALBEDO"])
    print(f"\nALBEDO  LCDB albedos not flagged assumed: {len(alb):,}; for NEOs typed from an assumed"
          f" albedo: {len(ga):,}")
    if len(ga):
        implied = np.asarray(CLASSES)[rule_class(ga.ALBEDO.to_numpy(float))]
        flip = implied != collapse(ga.comp_group)
        print(f"  their class would change on {int(flip.sum())} of {len(ga)};"
              f" source codes {ga.S_ALBEDO.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
