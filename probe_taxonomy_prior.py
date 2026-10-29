# -*- coding: utf-8 -*-
"""How likely is a guessed taxonomy to be right?

FINDINGS.md F2: 90.2% of NEO taxonomies are inferred, almost all from an
albedo that was itself ASSUMED from the orbit, and the model treats the guess
as a fact. F2 proposed a categorical prior over the class, conditioned on
albedo and orbit and calibrated on the bodies with a measured class.
THIRD-PASS.md observed that this is a multinomial logistic regression, which
statsmodels fits as `MNLogit` and which needs no sampler.

The builder (AsteroidCatalog `asteroid_catalog/enrich.py`) has two rules:

    albedo          a measured albedo below 0.10 is C, anything else S
    albedo_assumed  no albedo: the orbit bin's assumed albedo, then the same
                    split. For a NEO that means S inside 2.5 AU, C beyond

This probe takes the 170k bodies within 4.6 AU whose class came from a
source, holds 20% out, and on the held-out rows measures each rule against two
fitted models:

    orbit model     class ~ orbit bin, e, i, H, NEO flag     (the inputs the
                    albedo_assumed rule has)
    albedo model    the same plus log albedo                 (the inputs the
                    albedo rule has)

On NEOs the orbit model does not beat base rates, so the prior it applies to
the albedo_assumed NEOs and the top of the profitability ranking is the
empirical one instead: among the 2,943 NEOs with a source class, how the true
class falls for each class the rule would give, binned by H. A bin with fewer
than 30 bodies takes the pooled row for its rule.

Bodies with a measured class are biased toward bright S-types, which
spectroscopic surveys reach first, so this prior is too generous to S and too
hard on C. It is the only calibration data there is.

The fit is a numpy Newton-Raphson multinomial logit, not statsmodels, so the
probe keeps the repo's two dependencies. It is the textbook model, written
from the likelihood; nothing is copied. Before measuring it must recover known
coefficients from synthetic data, and the builder's albedo rule must reproduce
the 84.0% its own docstring records.

Classes are `comp_group` collapsed to four: C-complex, S-complex, X-complex,
and Other (D, V, L, K, A, Q and the rest). Bodies past 4.6 AU are left out:
the builder types them D from the orbit, which is a different rule.

This script READS ONLY. It reads the catalog and the profitability catalog,
does not import master, does not build a stage, and writes nothing to
`asteroid_pipeline/`.

Usage:
    py ../starred-repos/probe_taxonomy_prior.py
    py ../starred-repos/probe_taxonomy_prior.py --top 40 --seed 2
"""

import argparse

import numpy as np
import pandas as pd

CATALOG = "asteroid_pipeline/asteroid_catalog.csv"
PROFITABILITY = "asteroid_pipeline/profitability_catalog.csv"
CLASSES = ["C-complex", "S-complex", "X-complex", "Other"]

# Mirrored from AsteroidCatalog asteroid_catalog/derive.py (catalog 92f7739):
# the albedo assumed for a body with no albedo and no class, by semi-major
# axis, with a separate table that NEOs take in preference. (lower, upper, p).
ALBEDO_BY_A = [(0.0, 1.3, 0.1870), (1.3, 2.0, 0.4135), (2.0, 2.5, 0.1900),
               (2.5, 2.82, 0.0860), (2.82, 3.27, 0.0660), (3.27, 3.7, 0.0570),
               (3.7, 4.6, 0.0550), (4.6, 5.5, 0.0700), (5.5, 1e9, 0.0880)]
ALBEDO_BY_A_NEO = [(0.0, 1.3, 0.1870), (1.3, 2.0, 0.1700), (2.0, 2.5, 0.1370),
                   (2.5, 2.82, 0.0620), (2.82, 1e9, 0.0370)]
# enrich.py `_class_from_albedo`, and OUTER_SOLAR_SYSTEM_AU, past which an
# untyped body is D from its orbit.
C_BELOW_ALBEDO = 0.10
OUTER_AU = 4.6
# enrich.py's docstring: the split "is right for 84.0%" of the 65,159 bodies
# of the 2026-09-23 release with a source class and a measured albedo.
BUILDER_SPLIT_RIGHT = 0.840

BINS = [b[:2] for b in ALBEDO_BY_A if b[0] < OUTER_AU]

# H bins for the NEO prior. The last is open-ended and mostly empty of bodies
# with a measured class, so its probabilities are an extrapolation.
H_BINS = [0.0, 18.0, 20.0, 22.0, 25.0, 40.0]
MIN_CELL = 30


# ── the model ────────────────────────────────────────────────────────────────
def softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    p = np.exp(z)
    return p / p.sum(axis=1, keepdims=True)


def fit_mnlogit(X, y, k, ridge=1e-4, iters=50, tol=1e-9):
    """Multinomial logit by Newton-Raphson, first class the reference.

    X is (n, d) with an intercept column; y holds class indices 0..k-1.
    Returns W of shape (d, k) with W[:, 0] = 0. A small ridge keeps the
    Hessian invertible when a class is rare in some bin.
    """
    n, d = X.shape
    Y = np.eye(k)[y]
    w = np.zeros(d * (k - 1))
    for _ in range(iters):
        W = np.column_stack([np.zeros(d), w.reshape(k - 1, d).T])
        P = softmax(X @ W)
        grad = (X.T @ (Y - P))[:, 1:].T.ravel() - ridge * n * w
        H = np.zeros((d * (k - 1), d * (k - 1)))
        for a in range(1, k):
            for b in range(1, k):
                wts = P[:, a] * ((a == b) - P[:, b])
                H[(a - 1) * d:a * d, (b - 1) * d:b * d] = -(X * wts[:, None]).T @ X
        H -= ridge * n * np.eye(d * (k - 1))
        step = np.linalg.solve(H, grad)
        w -= step
        if np.abs(step).max() < tol:
            break
    return np.column_stack([np.zeros(d), w.reshape(k - 1, d).T])


def check_recovery(seed=0):
    """The fit must recover known coefficients from data drawn from them."""
    rng = np.random.default_rng(seed)
    n, d, k = 60_000, 4, 4
    X = np.column_stack([np.ones(n), rng.normal(size=(n, d - 1))])
    W = np.column_stack([np.zeros(d), rng.normal(scale=0.8, size=(d, k - 1))])
    P = softmax(X @ W)
    y = (rng.random(n)[:, None] > np.cumsum(P, axis=1)).sum(axis=1)
    W_hat = fit_mnlogit(X, y, k, ridge=0.0)
    err = np.abs(W_hat - W).max()
    ok = err < 0.06
    print(f"SYNTHETIC RECOVERY  n={n:,}: worst |W_hat - W| = {err:.4f} "
          f"({'OK' if ok else 'MISMATCH -- the fit does not recover its own model'})")
    return ok


# ── features and rules ──────────────────────────────────────────────────────
def assumed_albedo(a, neo):
    """The builder's assumed albedo for a body with neither albedo nor class."""
    out = np.full(len(a), np.nan)
    for table, mask in ((ALBEDO_BY_A, ~neo), (ALBEDO_BY_A_NEO, neo)):
        for lo, hi, p in table:
            sel = mask & (a >= lo) & (a < hi)
            out[sel] = p
    return out


def rule_class(albedo):
    """enrich.py's split: C below 0.10, S otherwise, as a CLASSES index."""
    return np.where(albedo < C_BELOW_ALBEDO, 0, 1)


def features(df, with_albedo):
    """Intercept, orbit-bin dummies (first bin the reference), e, i, H, NEO, NEO x a."""
    a = df["semi_major_axis_au"].to_numpy(float)
    cols = [np.ones(len(df))]
    for lo, hi in BINS[1:]:
        cols.append(((a >= lo) & (a < hi)).astype(float))
    e = df["eccentricity"].to_numpy(float)
    i = df["inclination_deg"].to_numpy(float)
    h = df["absolute_magnitude_h"].to_numpy(float)
    neo = df["neo"].to_numpy(float)
    cols += [(e - 0.15) / 0.1, (i - 8.0) / 6.0, (h - 15.0) / 2.0, neo, neo * (a - 2.0)]
    if with_albedo:
        la = np.log10(df["albedo"].to_numpy(float))
        cols += [(la + 1.0) / 0.3, ((la + 1.0) / 0.3) ** 2]
    return np.column_stack(cols)


def collapse(group):
    g = group.astype(str)
    return np.where(g.isin(CLASSES[:3]), g, "Other")


def load(path):
    cols = ["designation", "is_neo", "semi_major_axis_au", "eccentricity",
            "inclination_deg", "absolute_magnitude_h", "albedo", "spectral_type",
            "spectral_type_source", "comp_group", "comp_metal_fraction",
            "comp_carbon_fraction", "comp_ice_fraction"]
    c = pd.read_csv(path, usecols=cols, low_memory=False)
    c["designation"] = c["designation"].astype(str).str.strip()
    c["neo"] = c["is_neo"].astype(str).str.lower().isin(["true", "1"])
    need = ["semi_major_axis_au", "eccentricity", "inclination_deg", "absolute_magnitude_h"]
    c = c.dropna(subset=need)
    return c[c["semi_major_axis_au"] <= OUTER_AU].copy()


# ── scoring ──────────────────────────────────────────────────────────────────
def log_loss(P, y):
    return -np.mean(np.log(np.clip(P[np.arange(len(y)), y], 1e-12, None)))


def score(label, P, y, rule=None):
    """Accuracy and log loss on held-out rows, beside the rule's accuracy."""
    acc = (P.argmax(axis=1) == y).mean()
    base = np.bincount(y, minlength=len(CLASSES)) / len(y)
    line = (f"  {label:30s} n={len(y):7,}  model acc {acc * 100:5.1f}%"
            f"  log loss {log_loss(P, y):.3f} (base rates {log_loss(np.tile(base, (len(y), 1)), y):.3f})")
    if rule is not None:
        line += f"  rule acc {(rule == y).mean() * 100:5.1f}%"
    print(line)


def calibration(P, y, label):
    """Predicted probability of the true class against how often it was right."""
    p_top = P.max(axis=1)
    right = P.argmax(axis=1) == y
    print(f"  calibration, {label}: predicted top-class probability vs observed")
    for lo, hi in ((0.0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01)):
        s = (p_top >= lo) & (p_top < hi)
        if s.sum() >= 30:
            print(f"    {lo:.2f}-{min(hi, 1):.2f}  n={int(s.sum()):7,}"
                  f"  predicted {p_top[s].mean():.3f}  observed {right[s].mean():.3f}")


def main():
    ap = argparse.ArgumentParser(description="calibrated prior over the guessed taxonomy")
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--profitability", default=PROFITABILITY)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    if not check_recovery():
        return
    cat = load(args.catalog)
    src = cat[cat["spectral_type_source"].eq("source")].copy()
    src["y"] = pd.Categorical(collapse(src["comp_group"]), categories=CLASSES).codes
    print(f"\nbodies with a source class, a <= {OUTER_AU} AU: {len(src):,}"
          f"  ({int(src.neo.sum()):,} NEOs)")
    print("  " + "  ".join(f"{c} {int((src.y == j).sum()):,}" for j, c in enumerate(CLASSES)))

    # The mirrored rules must give the class the builder actually assigned.
    print("\nBUILDER CHECK  the mirrored rules against the classes the release assigned")
    ok = True
    for source, guess in (
            ("albedo", lambda d: rule_class(d.albedo.to_numpy(float))),
            ("albedo_assumed", lambda d: rule_class(assumed_albedo(
                d.semi_major_axis_au.to_numpy(float), d.neo.to_numpy(bool))))):
        d = cat[cat["spectral_type_source"].eq(source)]
        same = np.asarray(CLASSES)[guess(d)] == collapse(d["comp_group"])
        print(f"  {source:15s} {int(same.sum()):9,} / {len(d):9,}  ({same.mean() * 100:.3f}%)")
        ok &= same.mean() > 0.999
    if not ok:
        print("  MISMATCH -- these are not the builder's rules")
        return

    alb = src[(src["albedo"] > 0) & (src["albedo"] < 1)]
    c_or_not = (rule_class(alb.albedo.to_numpy(float)) == 0) == (alb.y.to_numpy() == 0)
    print(f"  albedo split, C versus not-C, on {len(alb):,} source-typed bodies with an albedo:"
          f" {c_or_not.mean() * 100:.1f}%. The builder's docstring records"
          f" {BUILDER_SPLIT_RIGHT * 100:.1f}% on the 2026-09-23 release, beyond 4.6 AU"
          " included and its scoring unrecorded")

    rng = np.random.default_rng(args.seed)
    test = rng.random(len(src)) < 0.2
    tr, te = src[~test], src[test]
    k = len(CLASSES)

    W_orbit = fit_mnlogit(features(tr, False), tr.y.to_numpy(), k)
    P_orbit = softmax(features(te, False) @ W_orbit)
    rule_orbit = rule_class(assumed_albedo(te.semi_major_axis_au.to_numpy(float),
                                           te.neo.to_numpy(bool)))

    tra = tr[(tr.albedo > 0) & (tr.albedo < 1)]
    tea_mask = ((te.albedo > 0) & (te.albedo < 1)).to_numpy()
    tea = te[tea_mask]
    W_alb = fit_mnlogit(features(tra, True), tra.y.to_numpy(), k)
    P_alb = softmax(features(tea, True) @ W_alb)

    y_te = te.y.to_numpy()
    neo_te = te.neo.to_numpy(bool)
    print(f"\nHELD OUT  (20%, seed {args.seed}; four classes, so a rule that says only C or S"
          " is wrong on every X and Other)")
    score("orbit model, all", P_orbit, y_te, rule_orbit)
    score("orbit model, NEOs", P_orbit[neo_te], y_te[neo_te], rule_orbit[neo_te])
    score("albedo model, all", P_alb, tea.y.to_numpy(), rule_class(tea.albedo.to_numpy(float)))
    na = tea.neo.to_numpy(bool)
    score("albedo model, NEOs", P_alb[na], tea.y.to_numpy()[na],
          rule_class(tea.albedo.to_numpy(float))[na])
    calibration(P_orbit, y_te, "orbit model, all")
    calibration(P_orbit[neo_te], y_te[neo_te], "orbit model, NEOs")

    # For NEOs the fitted model does not beat base rates, so the prior applied
    # to them is the empirical one: among NEOs with a source class, how the
    # true class falls for each class the albedo_assumed rule would give, by H.
    neo_src = src[src.neo].copy()
    neo_src["rule"] = rule_class(assumed_albedo(neo_src.semi_major_axis_au.to_numpy(float),
                                                np.ones(len(neo_src), bool)))
    neo_src["hbin"] = pd.cut(neo_src.absolute_magnitude_h, H_BINS, right=False)
    print(f"\nNEOs WITH A SOURCE CLASS  (n = {len(neo_src):,}): the true class, by the class"
          " the albedo_assumed rule would give")
    table = {}
    for r in (0, 1):
        pooled_rows = neo_src[neo_src.rule == r]
        pooled = np.bincount(pooled_rows.y, minlength=k) / len(pooled_rows)
        for hb, s in pooled_rows.groupby("hbin", observed=False):
            frac = np.bincount(s.y, minlength=k) / max(len(s), 1)
            table[(r, hb)] = frac if len(s) >= MIN_CELL else pooled
            flag = "" if len(s) >= MIN_CELL else "   (too few: the pooled row stands in)"
            print(f"  rule {CLASSES[r][0]}, H {str(hb):13s} n={len(s):5,}  "
                  + "  ".join(f"{CLASSES[j][0]} {frac[j]:.3f}" for j in range(k)) + flag)
        print(f"  rule {CLASSES[r][0]}, pooled         n={len(pooled_rows):5,}  "
              + "  ".join(f"{CLASSES[j][0]} {pooled[j]:.3f}" for j in range(k)))

    target = cat[cat["spectral_type_source"].eq("albedo_assumed") & cat["neo"]].copy()
    target["rule"] = rule_class(assumed_albedo(target.semi_major_axis_au.to_numpy(float),
                                               np.ones(len(target), bool)))
    target["hbin"] = pd.cut(target.absolute_magnitude_h, H_BINS, right=False)
    probs = np.vstack([table.get((r, hb), np.full(k, np.nan))
                       for r, hb in zip(target.rule, target.hbin)])
    for j, c in enumerate(CLASSES):
        target["p_" + c[0]] = probs[:, j]
    target["p_assigned"] = probs[np.arange(len(target)), target.rule.to_numpy()]
    beyond = target.absolute_magnitude_h >= H_BINS[-2]
    print(f"\nALBEDO_ASSUMED NEOs  (n = {len(target):,}; {int(beyond.sum()):,} at H >= {H_BINS[-2]},"
          f" a bin {int(((neo_src.absolute_magnitude_h >= H_BINS[-2])).sum())} source-typed NEOs fill)")
    for r, s in target.groupby("rule"):
        print(f"  assigned {CLASSES[r]:10s} n={len(s):7,}  P(assigned) median {s.p_assigned.median():.3f}"
              f"   mean P(C) {s.p_C.mean():.3f}  P(S) {s.p_S.mean():.3f}"
              f"  P(X) {s.p_X.mean():.3f}  P(Other) {s.p_O.mean():.3f}")

    # Expected composition under that prior, against the point estimate. Each
    # class's fractions are its mean over the source rows of that class.
    print("\n  value-bearing fractions, point estimate vs expectation under the prior:")
    for f in ("metal", "carbon", "ice"):
        col = f"comp_{f}_fraction"
        means = np.array([pd.to_numeric(src.loc[src.y == j, col], errors="coerce").mean()
                          for j in range(k)])
        expect = probs @ means
        point = pd.to_numeric(target[col], errors="coerce").to_numpy(float)
        for r in (0, 1):
            m = target.rule.to_numpy() == r
            print(f"    {f:6s} assigned {CLASSES[r][0]}: point median {np.nanmedian(point[m]):.3f}"
                  f"  expected median {np.nanmedian(expect[m]):.3f}"
                  + ("   class means " + "  ".join(f"{CLASSES[j][0]} {means[j]:.3f}"
                                                   for j in range(k)) if r == 0 else ""))

    prof = pd.read_csv(args.profitability, low_memory=False, usecols=["designation"])
    prof["designation"] = prof["designation"].astype(str).str.strip()
    prof["rank"] = np.arange(1, len(prof) + 1)
    top = prof.head(1000).merge(target, on="designation", how="left")
    print("\nPROFITABILITY RANKING  (albedo_assumed NEOs among the top rows)")
    for n in (25, 100, 1000):
        t = top.head(n)
        t = t[t.p_assigned.notna()]
        print(f"  top {n:5,}: {len(t):4,} albedo_assumed NEOs"
              f"   H >= {H_BINS[-2]}: {int((t.absolute_magnitude_h >= H_BINS[-2]).sum()):4,}"
              f"   P(assigned) median {t.p_assigned.median():.3f}")
    print(f"\n  the top {args.top}:")
    cols = ["rank", "designation", "comp_group", "semi_major_axis_au",
            "absolute_magnitude_h", "p_C", "p_S", "p_X", "p_O", "p_assigned"]
    print(top.head(args.top)[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
