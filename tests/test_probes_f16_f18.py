"""Checks behind F16 (albedo coverage), F17 (source regions) and F18 (lightcurves), and their planted defects."""

import json
import os
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from tests._util import fixture, quiet

import probe_albedo_coverage as pa
import probe_lightcurves as pl
import probe_source_regions as ps

HAVE_GRANVIK = os.path.exists(os.path.join(ps.DEFAULT_CACHE, ps.NAME))


def datacloud_answer():
    """Two bodies for one name, the way SsODNet answers 2019 XG1."""
    return json.dumps({"data": [
        {"identity": {"name": "2019 XG1", "number": None, "aliases": ["2007 JN4", "K19X01G"]},
         "datacloud": {"diamalbedo": []}},
        {"identity": {"name": "Durer", "number": 3104, "aliases": [3104, "1955 XG1"]},
         "datacloud": {"diamalbedo": [{"albedo": 0.30, "method": "NEATM", "shortbib": "Masiero+2011"},
                                      {"albedo": 0.28, "method": "TPM", "shortbib": "Hung+2022"}]}},
    ]})


class AlbedoCoverage(unittest.TestCase):

    def test_datacloud_keeps_only_the_body_asked_for(self):
        got, others = pa.parse_datacloud(datacloud_answer(), "2019 XG1")
        self.assertEqual(got, [])
        self.assertEqual(others, 1)

    def test_datacloud_finds_the_body_by_alias(self):
        got, _ = pa.parse_datacloud(datacloud_answer(), "1955 XG1")
        self.assertEqual(sorted(a for a, _, _ in got), [0.28, 0.30])

    def test_planted_a_loose_identity_check_takes_the_wrong_body(self):
        with mock.patch.object(pa, "same_body", lambda identity, designation: True):
            got, _ = pa.parse_datacloud(datacloud_answer(), "2019 XG1")
        self.assertEqual(len(got), 2)

    def test_neowise_join_uses_both_designations(self):
        nw = pd.DataFrame({"designation": ["2005 AB1", "77"], "nw_albedo": [0.05, 0.2],
                           "nw_fits": [1, 1], "nw_type": ["neos", "neos"]})
        df = pd.DataFrame({"designation": ["613939", "77", "999"],
                           "provisional_designation": ["2005 AB1", "nan", "2001 ZZ9"],
                           "comp_group": ["S-complex", "C-complex", "S-complex"]})
        hit = pa.in_neowise(df, set(nw.designation))
        self.assertEqual(list(hit), [True, True, False])
        m = pa.neowise_matches(df, nw)
        self.assertEqual(sorted(m.designation), ["613939", "77"])


class Lightcurves(unittest.TestCase):

    def setUp(self):
        with open(fixture("lcdb_sample.txt"), encoding="latin-1") as fh:
            self.text = fh.read()

    def test_parse_reproduces_known_rows(self):
        self.assertTrue(quiet(pl.check_parse, pl.parse_summary(self.text)))

    def test_repeated_header_names_take_the_next_field(self):
        names = [n for n, _, _ in pl.header_columns(
            "NUMBER  S CLASS  S F DIA.  F PERIOD", "-" * 36)]
        self.assertEqual(names, ["NUMBER", "S_CLASS", "CLASS", "S_DIA.", "F_DIA.", "DIA.",
                                 "F_PERIOD", "PERIOD"])

    def test_reliable_drops_limits_and_u1(self):
        lc = pd.DataFrame({"U_num": [3, 2, 1, 3], "PERIOD": [2.0, 3.0, 4.0, 5.0],
                           "F_PERIOD": ["", "", "", ">"]})
        self.assertEqual(list(pl.reliable(lc).PERIOD), [2.0, 3.0])

    def test_planted_a_shifted_header_misreads_the_table(self):
        lines = self.text.splitlines()
        k = next(i for i, l in enumerate(lines) if l.startswith("NUMBER"))
        lines[k] = lines[k].replace("NUMBER ", "NUMBER  ", 1)       # every later column one right
        self.assertFalse(quiet(pl.check_parse, pl.parse_summary("\n".join(lines))))


def synthetic_model(n=50, seed=0):
    rng = np.random.default_rng(seed)
    t = np.zeros((n, 26))
    t[:, 0:4] = rng.uniform([0.5, 0.0, 0.0, 15.0], [3.0, 0.9, 40.0, 25.0], (n, 4))
    src = rng.random((n, 7))
    t[:, ps.COL_SRC] = src / src.sum(axis=1, keepdims=True)
    t[:, ps.COL_DARK] = rng.uniform(0, 0.5, n)
    t[:, ps.COL_MID] = rng.uniform(0, 0.5, n)
    return t


class SourceRegions(unittest.TestCase):

    def test_model_check_passes_a_consistent_table(self):
        with mock.patch.object(ps, "N_CELLS", 50):
            self.assertTrue(quiet(ps.check_model, synthetic_model()))

    def test_planted_source_probabilities_that_do_not_sum_to_one(self):
        t = synthetic_model()
        t[3, 6] += 0.01
        with mock.patch.object(ps, "N_CELLS", 50):
            self.assertFalse(quiet(ps.check_model, t))

    def test_planted_a_missing_row(self):
        with mock.patch.object(ps, "N_CELLS", 51):
            self.assertFalse(quiet(ps.check_model, synthetic_model()))

    @unittest.skipUnless(HAVE_GRANVIK, "Granvik model not downloaded")
    def test_real_model(self):
        self.assertTrue(quiet(ps.check_model, ps.load_model(ps.DEFAULT_CACHE)))


if __name__ == "__main__":
    unittest.main()
