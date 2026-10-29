"""Shared helpers: the repo root on sys.path, quiet calls, fixtures, a synthetic Toliou grid."""

import contextlib
import io
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(ROOT, "tests", "fixtures")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def quiet(f, *args, **kwargs):
    """Call a probe's check without its printout."""
    with contextlib.redirect_stdout(io.StringIO()):
        return f(*args, **kwargs)


def fixture(name):
    return os.path.join(FIXTURES, name)


def builder_rows():
    """A stratified sample of data-2026-10-21b's inferred rows; see tests/README.md."""
    d = pd.read_csv(fixture("builder_rows.csv"))
    d["neo"] = d["is_neo"].astype(str).str.lower().isin(["true", "1"])
    return d


def synthetic_toliou(seed=0, fill=0.4):
    """A small regular grid shaped like Toliou et al.'s, partly filled.

    Axes are coarser than the real table's but regular, the cumulative
    columns are monotone and end at 1, and only `fill` of the cells exist,
    so the empty-cell fallback is exercised.
    """
    rng = np.random.default_rng(seed)
    a = np.round(np.arange(0.55, 2.0, 0.1), 3)
    e = np.round(np.arange(0.02, 0.9, 0.08), 3)
    i = np.arange(2.0, 40.0, 4.0)
    h = np.arange(15.125, 25.0, 1.0)
    cells = np.array(np.meshgrid(a, e, i, h, indexing="ij")).reshape(4, -1).T
    cells = cells[rng.random(len(cells)) < fill]
    steps = rng.random((len(cells), 26))
    cum = np.cumsum(steps, axis=1)
    cum /= cum[:, -1:]
    return np.column_stack([cells, cum])


def synthetic_neos(n=300, seed=3):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"semi_major_axis_au": rng.uniform(0.5, 2.1, n),
                         "eccentricity": rng.uniform(0.0, 0.9, n),
                         "inclination_deg": rng.uniform(0, 45, n),
                         "absolute_magnitude_h": rng.uniform(14, 26, n)})
