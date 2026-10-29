# -*- coding: utf-8 -*-
"""Do the NEOs typed from an assumed albedo actually have a measured one somewhere?

F11 found that a class guessed from an ASSUMED albedo is right about half the
time for a NEO, against about 70% for one guessed from a MEASURED albedo. The
builder assumes an albedo for 38,616 of the 42,775 NEOs, because none of its
sources gave one. If some of them do have a measured albedo somewhere, the
cheapest improvement to F2 is to use it.

Two independent places to look:

    NEOWISE V2.0    the table the builder already reads (IRSA TAP,
                    neowisesbpropv2). Re-read here directly, all 183k rows of
                    every orbit type, and joined on number and provisional
                    designation against both the catalog's designation and its
                    provisional_designation: V2.0 predates many numberings
                    and lists those bodies only provisionally. A match among
                    the albedo_assumed NEOs would be a merge the builder
                    missed.
    SsODNet datacloud  every published diameter and albedo SsODNet knows for a
                    body, one request per body, for a sample: every
                    albedo_assumed NEO in the top 1,000 rows of the ranking,
                    and 300 more drawn at random. SsODNet's summary card, which
                    the builder reads, carries one best value; the datacloud
                    carries all of them.

For every albedo found, the probe applies the builder's split (C below 0.10, S
otherwise) and counts how many NEOs would change class.

Checks before measuring: the mirrored builder rules must reproduce the class
the release assigned on every inferred row (probe_taxonomy_prior.check_builder),
and the NEOWISE join must find at least 95% of the NEOs the catalog matched
to NEOWISE (`source_neowise`), which shows the join keys work before an
absence of matches is believed.

Network reads go through snapshots.py. This script READS ONLY and writes
nothing to `asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_albedo_coverage.py
    py ../starred-repos/probe_albedo_coverage.py --random 600 --seed 2
"""

import argparse
import io
import json
import time
import urllib.error
import urllib.parse

import numpy as np
import pandas as pd

import snapshots
from probe_taxonomy_prior import (
    C_BELOW_ALBEDO, CLASSES, assumed_albedo, check_builder, collapse, rule_class)

CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
# Every row, whatever its orbit type: NEOWISE typed some of today's NEOs as
# main-belt when it fitted them.
NEOWISE_QUERY = "SELECT asteroid_number, prov_desig, v_albedo, type FROM neowisesbpropv2"
DATACLOUD = ("https://ssp.imcce.fr/webservices/ssodnet/api/datacloud.php"
             "?-name={}&-resource=diamalbedo&-mime=json")
PAUSE_S = 0.25


def neowise_rows():
    url = TAP + "?" + urllib.parse.urlencode({"QUERY": NEOWISE_QUERY, "FORMAT": "csv"})
    d = pd.read_csv(io.BytesIO(snapshots.fetch(url, timeout=600)))
    d["v_albedo"] = pd.to_numeric(d["v_albedo"], errors="coerce")
    num = pd.to_numeric(d["asteroid_number"], errors="coerce").astype("Int64").astype(str)
    keyed = pd.concat([
        pd.DataFrame({"designation": num, "nw_albedo": d["v_albedo"], "nw_type": d["type"]}),
        pd.DataFrame({"designation": d["prov_desig"].astype(str).str.strip(),
                      "nw_albedo": d["v_albedo"], "nw_type": d["type"]}),
    ])
    keyed = keyed[~keyed["designation"].isin(["<NA>", "nan", ""])]
    # Every designation NEOWISE lists, albedo or not (for the join check), and
    # the ones with an albedo (for the measurement).
    listed = set(keyed["designation"])
    with_albedo = keyed.dropna(subset=["nw_albedo"]).groupby("designation", as_index=False).agg(
        nw_albedo=("nw_albedo", "median"), nw_fits=("nw_albedo", "size"),
        nw_type=("nw_type", "first"))
    return listed, with_albedo


def in_neowise(df, keys):
    """Whether each row joins NEOWISE on its designation or its provisional one."""
    return df.designation.isin(keys) | df.provisional_designation.isin(keys)


def neowise_matches(df, nw):
    """Rows of `df` with a NEOWISE albedo, keyed either way, one row per body."""
    by_prov = (df.drop(columns="designation")
               .rename(columns={"provisional_designation": "designation"})
               .merge(nw, on="designation", how="inner"))
    by_prov["designation"] = df.set_index("provisional_designation").loc[
        by_prov["designation"], "designation"].to_numpy()
    both = pd.concat([df.drop(columns="provisional_designation")
                      .merge(nw, on="designation", how="inner"), by_prov], ignore_index=True)
    return both.drop_duplicates("designation")


def same_body(identity, designation):
    """Whether a datacloud identity is the body asked for.

    SsODNet resolves names loosely and can answer with several bodies: asked
    for 2019 XG1 it also returns (3104) Durer, which carries the alias
    1955 XG1, with eight albedos. Only the body whose name, number or alias
    is the designation itself counts.
    """
    want = str(designation).strip()
    names = {str(identity.get("name", "")), str(identity.get("number", ""))}
    names |= {str(a) for a in identity.get("aliases") or []}
    return want in names


def parse_datacloud(text, designation):
    """Every albedo for this body in a datacloud answer, as [(albedo, method, shortbib)].

    Returns (albedos, n_other_bodies), the second counting answers that were
    for some other body and were dropped.
    """
    try:
        body = json.loads(text)
    except ValueError:
        return None
    out, others = [], 0
    for obj in body.get("data") or []:
        if not same_body(obj.get("identity") or {}, designation):
            others += 1
            continue
        for row in (obj.get("datacloud") or {}).get("diamalbedo") or []:
            a = row.get("albedo")
            if isinstance(a, (int, float)) and 0 < a < 1:
                out.append((float(a), row.get("method"), row.get("shortbib")))
    return out, others


def datacloud(designation):
    url = DATACLOUD.format(urllib.parse.quote(str(designation)))
    try:
        return parse_datacloud(snapshots.fetch(url, timeout=60).decode("utf-8", "replace"),
                               designation)
    except (urllib.error.URLError, TimeoutError):
        return None


def load(path):
    cols = ["designation", "is_neo", "semi_major_axis_au", "albedo", "albedo_provider",
            "spectral_type_source", "comp_group", "absolute_magnitude_h", "source_neowise",
            "provisional_designation"]
    c = pd.read_csv(path, usecols=cols, low_memory=False)
    c["designation"] = c["designation"].astype(str).str.strip()
    c["provisional_designation"] = c["provisional_designation"].astype(str).str.strip()
    c["neo"] = c["is_neo"].astype(str).str.lower().isin(["true", "1"])
    return c


def flips(df, albedo_col):
    """Assigned class (the builder's) against the class the measured albedo implies."""
    assigned = collapse(df["comp_group"])
    implied = np.asarray(CLASSES)[rule_class(df[albedo_col].to_numpy(float))]
    return pd.crosstab(pd.Series(assigned, name="assigned"), pd.Series(implied, name="implied"))


def main():
    ap = argparse.ArgumentParser(description="measured albedos for albedo-guessed NEOs")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--random", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    cat = load(args.catalog)
    if not check_builder(cat[cat["semi_major_axis_au"] <= 4.6]):
        return
    neo = cat[cat.neo]
    guessed = neo[neo.spectral_type_source.eq("albedo_assumed")].copy()
    print(f"\nNEOs {len(neo):,}; typed from an assumed albedo {len(guessed):,}")
    print(snapshots.describe())

    listed, nw = neowise_rows()
    # The builder prefers JPL's copy of a NEOWISE albedo, so `albedo_provider`
    # rarely says NEOWISE; `source_neowise` marks every body NEOWISE matched.
    from_nw = neo[neo["source_neowise"].astype(str).str.lower().isin(["true", "1"])]
    hit = in_neowise(from_nw, listed)
    found = hit.mean() if len(from_nw) else np.nan
    print(f"\nNEOWISE V2.0  {len(listed):,} designations listed, {len(nw):,} with an albedo,"
          " all orbit types")
    print(f"  join check: NEOs the catalog matched to NEOWISE, found again:"
          f" {int(hit.sum()):,} / {len(from_nw):,}"
          f" ({found * 100:.1f}%) {'(OK)' if found >= 0.95 else '-- MISMATCH: the join keys do not work'}")
    if not found >= 0.95:
        return
    m = neowise_matches(guessed, nw)
    print(f"  albedo_assumed NEOs with a NEOWISE V2.0 albedo: {len(m):,}")
    if len(m):
        print(flips(m, "nw_albedo").to_string())

    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"], nrows=1000)
    prof["designation"] = prof["designation"].astype(str).str.strip()
    top = guessed[guessed.designation.isin(prof.designation)]
    rest = guessed[~guessed.designation.isin(prof.designation)]
    sample = pd.concat([top.assign(group="top 1,000"),
                        rest.sample(n=min(args.random, len(rest)), random_state=args.seed)
                        .assign(group="random")])
    print(f"\nSsODNet DATACLOUD  {len(top)} albedo_assumed NEOs in the top 1,000 rows,"
          f" plus {len(sample) - len(top)} at random")
    rows, failed, other_bodies = [], 0, 0
    for _, b in sample.iterrows():
        got = datacloud(b.designation)
        if not snapshots.replay_folder():
            time.sleep(PAUSE_S)
        if got is None:
            failed += 1
            continue
        got, others = got
        other_bodies += others > 0
        alb = [a for a, _, _ in got]
        rows.append({"designation": b.designation, "group": b.group, "H": b.absolute_magnitude_h,
                     "comp_group": b.comp_group, "n_albedo": len(alb),
                     "dc_albedo": float(np.median(alb)) if alb else np.nan,
                     "methods": ";".join(sorted({str(mt) for _, mt, _ in got})),
                     "refs": ";".join(sorted({str(r) for _, _, r in got}))})
    r = pd.DataFrame(rows)
    print(f"  answered {len(r)}; failed {failed}; answers that also carried some other body,"
          f" dropped: {other_bodies}")
    for g, s in r.groupby("group"):
        hit = s[s.n_albedo > 0]
        print(f"  {g:10s} n={len(s):4}  with any albedo {len(hit):4} ({len(hit) / len(s) * 100:.1f}%)")
    hit = r[r.n_albedo > 0]
    if len(hit):
        print("\n  class the datacloud albedo implies, against the builder's:")
        print(flips(hit, "dc_albedo").to_string())
        print("\n  bodies with a datacloud albedo:")
        print(hit[["designation", "group", "H", "comp_group", "dc_albedo", "n_albedo", "methods",
                   "refs"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\n  (the builder's split: C below albedo {C_BELOW_ALBEDO}, S otherwise)")


if __name__ == "__main__":
    main()
