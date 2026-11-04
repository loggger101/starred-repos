"""run_all.py is itself a check: it must fail when a probe does, and figures.json must hold."""

import contextlib
import glob
import io
import json
import os
import re
import sys
import tempfile
import unittest
from unittest import mock

from tests._util import ROOT

import run_all

# Figures whose probe output the note trimmed or paraphrased, so the pattern
# cannot read them back from the note. The number itself must still be there.
NOT_VERBATIM = {
    "survey rows (1-in-40 stride)",
    "comp_group agreement with SDSS, %",
    "top 1,000 carrying water, mean P(q_min<0.2)",
    "top-100 NEOs fetched from NEODyS",
}


def note_text(name):
    with open(os.path.join(ROOT, name), encoding="utf-8") as fh:
        return fh.read()


def in_note(figure):
    """Whether a figure's recorded value reads back out of the note that records it."""
    text = note_text(figure["note"])
    if figure["label"] in NOT_VERBATIM:
        return (f"{figure['recorded']:,}" in text) or (str(figure["recorded"]) in text)
    return run_all.extract(figure["pattern"], text) == figure["recorded"]


class Spec(unittest.TestCase):

    def setUp(self):
        self.spec = run_all.load_spec()

    def test_every_probe_is_listed_and_has_a_figure(self):
        scripts = {p["script"] for p in self.spec["probes"].values()}
        on_disk = {os.path.basename(p) for p in glob.glob(os.path.join(ROOT, "probe_*.py"))}
        self.assertEqual(scripts, on_disk)
        self.assertEqual({f["probe"] for f in self.spec["figures"]}, set(self.spec["probes"]))

    def test_every_pattern_has_one_group_and_every_figure_a_date(self):
        for f in self.spec["figures"]:
            self.assertEqual(re.compile(f["pattern"]).groups, 1, f["label"])
            self.assertRegex(f["recorded_on"], r"^\d{4}-\d{2}-\d{2}$")

    def test_every_figure_reads_back_out_of_its_note(self):
        for f in self.spec["figures"]:
            self.assertTrue(in_note(f), f"{f['finding']} {f['label']}")

    def test_planted_a_moved_figure_is_caught(self):
        for f in self.spec["figures"]:
            self.assertFalse(in_note(dict(f, recorded=f["recorded"] + 12345.5)), f["label"])

    def test_planted_a_broken_pattern_is_caught(self):
        verbatim = [f for f in self.spec["figures"] if f["label"] not in NOT_VERBATIM]
        for f in verbatim:
            self.assertFalse(in_note(dict(f, pattern="NO SUCH LINE (\\d+)")), f["label"])

    def test_not_verbatim_lists_only_real_figures(self):
        self.assertLessEqual(NOT_VERBATIM, {f["label"] for f in self.spec["figures"]})


class Extract(unittest.TestCase):

    def test_reads_separators_and_signs(self):
        self.assertEqual(run_all.extract(r"n = ([\d,]+)", "n = 39,398 rows"), 39398.0)
        self.assertEqual(run_all.extract(r"gap ([+\-\d.]+)", "gap +1.296 km/s"), 1.296)

    def test_a_missing_figure_is_none_not_zero(self):
        self.assertIsNone(run_all.extract(r"n = (\d+)", "nothing here"))


class Offline(unittest.TestCase):

    def test_an_api_probe_needs_the_network(self):
        self.assertTrue(run_all.needs_network({"network": True}))

    def test_a_table_only_needs_it_until_it_is_cached(self):
        with tempfile.TemporaryDirectory() as cache:
            probe = {"download": "table.gz"}
            self.assertTrue(run_all.needs_network(probe, cache))
            open(os.path.join(cache, "table.gz"), "wb").close()
            self.assertFalse(run_all.needs_network(probe, cache))

    def test_every_probe_that_downloads_says_so(self):
        probes = run_all.load_spec()["probes"]
        for key, p in probes.items():
            with open(os.path.join(ROOT, p["script"]), encoding="utf-8") as fh:
                # probe_taxonomy's table is committed; a fresh checkout has it.
                if "urlretrieve" in fh.read() and key != "taxonomy":
                    self.assertIn("download", p, key)


def fake_probe(folder, body):
    path = os.path.join(folder, "fake_probe.py")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    return path


class Timeout(unittest.TestCase):

    def test_the_default_when_a_probe_names_none(self):
        self.assertEqual(run_all.timeout_for({}), run_all.DEFAULT_TIMEOUT_S)

    def test_a_probe_can_ask_for_longer(self):
        self.assertEqual(run_all.timeout_for({"timeout": 14400}), 14400)

    def test_the_command_line_wins(self):
        self.assertEqual(run_all.timeout_for({"timeout": 14400}, 60), 60)

    def test_a_probe_over_the_default_says_so(self):
        # The live F20 run outlasts the default, and a timeout loses every figure.
        self.assertGreater(run_all.timeout_for(run_all.load_spec()["probes"]["nhats_legs"]),
                           run_all.DEFAULT_TIMEOUT_S)


class RunProbe(unittest.TestCase):

    def run_fake(self, body, timeout=60):
        with tempfile.TemporaryDirectory() as folder:
            return run_all.run_probe(folder, fake_probe(folder, body), [], timeout)

    def test_a_clean_probe_is_ok(self):
        out, _, status = self.run_fake("print('n = 5')")
        self.assertEqual(status, "ok")
        self.assertIn("n = 5", out)

    def test_a_printed_mismatch_fails_even_on_exit_zero(self):
        self.assertEqual(self.run_fake("print('worst 1.0 (MISMATCH -- no)')")[2], "self-check failed")

    def test_a_crash_fails(self):
        self.assertEqual(self.run_fake("raise SystemExit(3)")[2], "exit 3")

    def test_a_hung_probe_times_out(self):
        self.assertEqual(self.run_fake("import time; time.sleep(30)", timeout=1)[2], "timeout")


class Main(unittest.TestCase):
    """The exit status is what a caller reads, so it must carry each failure."""

    def run_main(self, body, pattern=r"n = (\d+)"):
        with tempfile.TemporaryDirectory() as folder:
            repo = os.path.join(folder, "repo", "asteroid_pipeline")
            os.makedirs(repo)
            open(os.path.join(repo, "asteroid_catalog.csv"), "w").close()
            spec = {"probes": {"fake": {"script": fake_probe(folder, body), "args": []}},
                    "figures": [{"probe": "fake", "finding": "F0", "label": "a figure",
                                 "pattern": pattern, "recorded": 5, "recorded_on": "2026-11-04"}]}
            argv = ["run_all.py", "--repo", os.path.dirname(repo), "--out",
                    os.path.join(folder, "runs")]
            with mock.patch.object(run_all, "load_spec", lambda: spec), \
                    mock.patch.object(sys, "argv", argv), \
                    contextlib.redirect_stdout(io.StringIO()):
                try:
                    run_all.main()
                    code = 0
                except SystemExit as err:
                    code = err.code
            names = os.listdir(os.path.join(folder, "runs"))
            with open(os.path.join(folder, "runs", names[0]), encoding="utf-8") as fh:
                return code, json.load(fh)

    def test_a_clean_run_exits_zero_and_records_the_figure(self):
        code, run = self.run_main("print('n = 7')")
        self.assertEqual(code, 0)
        self.assertEqual(run["figures"][0]["now"], 7.0)

    def test_a_missing_figure_fails_the_run(self):
        code, run = self.run_main("print('nothing')")
        self.assertEqual(code, 1)
        self.assertIsNone(run["figures"][0]["now"])

    def test_a_failed_self_check_fails_the_run(self):
        code, run = self.run_main("print('MISMATCH'); print('n = 5')")
        self.assertEqual(code, 1)
        self.assertEqual(run["probes"]["fake"]["status"], "self-check failed")


if __name__ == "__main__":
    unittest.main()
