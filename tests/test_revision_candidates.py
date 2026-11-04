"""The revision candidates register: each row well formed, in the right file, and pinned."""

import csv
import datetime
import glob
import os
import re
import subprocess
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN = os.path.join(ROOT, "revision_candidates.csv")
CLOSED = os.path.join(ROOT, "revision_candidates_closed.csv")

COLUMNS = ["id", "status", "kind", "finding", "found_date", "target_repo", "target_file",
           "target_row", "field", "current_value", "proposed", "note", "evidence",
           "checked_against", "checked_date"]
LIVE = {"open", "blocked"}
SETTLED = {"applied", "declined", "superseded", "withdrawn"}
KINDS = {"value", "method", "wording", "citation"}


def read(path):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))
    return rows[0], [dict(zip(rows[0], r)) for r in rows[1:]]


def slug(heading):
    """GitHub's anchor rule, as check.yml applies it."""
    return re.sub(r"[^\w\- ]", "", heading.strip().lower()).replace(" ", "-")


def known_anchors():
    """Every note's headings, {file: {anchor}}."""
    anchors = {}
    for path in glob.glob(os.path.join(ROOT, "*.md")):
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        anchors[os.path.basename(path)] = {slug(h) for h in re.findall(r"^#+\s+(.*)$", text, re.M)}
    return anchors


def known_findings():
    """The finding ids in README's table."""
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as fh:
        return set(re.findall(r"^\| (F\d+) \|", fh.read(), re.M))


def is_date(text):
    try:
        return datetime.date.fromisoformat(text)
    except ValueError:
        return None


def problems(open_rows, closed_rows, findings, anchors, today=None):
    """What is wrong with the two files, as (code, id) pairs."""
    today = today or datetime.date.today()
    found, seen = [], set()
    for rows, live in ((open_rows, True), (closed_rows, False)):
        for r in rows:
            cid = r["id"]
            if any(not r[c].strip() for c in COLUMNS):
                found.append(("empty", cid))
            if not re.fullmatch(r"sr-\d{3}", cid) or cid in seen:
                found.append(("id", cid))
            seen.add(cid)
            if r["status"] not in LIVE | SETTLED:
                found.append(("status", cid))
            elif (r["status"] in LIVE) != live:
                found.append(("wrong-file", cid))
            if r["kind"] not in KINDS:
                found.append(("kind", cid))
            if any(f != "-" and f not in findings for f in r["finding"].split(";")):
                found.append(("finding", cid))
            note, _, anchor = r["note"].partition("#")
            if anchor not in anchors.get(note, ()):
                found.append(("note", cid))
            if not re.fullmatch(re.escape(r["target_repo"]) + r"@[0-9a-f]{40}", r["checked_against"]):
                found.append(("pin", cid))
            first, last = is_date(r["found_date"]), is_date(r["checked_date"])
            if not first or not last or first > last or last > today:
                found.append(("date", cid))
    return found


class Register(unittest.TestCase):

    def setUp(self):
        self.open_header, self.open = read(OPEN)
        self.closed_header, self.closed = read(CLOSED)

    def test_both_files_have_the_same_columns(self):
        self.assertEqual(self.open_header, COLUMNS)
        self.assertEqual(self.closed_header, COLUMNS)

    def test_every_row_is_well_formed(self):
        self.assertEqual(problems(self.open, self.closed, known_findings(), known_anchors()), [])

    def test_each_pin_is_a_commit_of_its_repo(self):
        for r in self.open + self.closed:
            repo, _, sha = r["checked_against"].partition("@")
            checkout = os.path.join(os.path.dirname(ROOT), repo)
            if not os.path.isdir(os.path.join(checkout, ".git")):
                self.skipTest(f"no {repo} checkout beside this repo")
            done = subprocess.run(["git", "-C", checkout, "cat-file", "-e", sha + "^{commit}"],
                                  capture_output=True)
            self.assertEqual(done.returncode, 0, f"{r['id']}: {sha} is not in {repo}")


class Planted(unittest.TestCase):
    """One defect for each rule, on a row that passes."""

    def row(self, **change):
        good = dict.fromkeys(COLUMNS, "x")
        good.update(id="sr-001", status="open", kind="value", finding="F1;-", found_date="2026-09-03",
                    target_repo="economicspace", note="FINDINGS.md#f1",
                    checked_against="economicspace@" + "a" * 40, checked_date="2026-11-04")
        return dict(good, **change)

    def codes(self, open_rows=(), closed_rows=()):
        return {code for code, _ in problems(open_rows, closed_rows, {"F1"}, {"FINDINGS.md": {"f1"}},
                                             today=datetime.date(2026, 11, 4))}

    def test_the_good_row_passes(self):
        self.assertEqual(self.codes([self.row()]), set())

    def test_an_empty_field(self):
        self.assertIn("empty", self.codes([self.row(evidence=" ")]))

    def test_general_researchs_id_is_refused(self):
        self.assertIn("id", self.codes([self.row(id="rc-001")]))

    def test_an_id_used_in_both_files(self):
        self.assertIn("id", self.codes([self.row()], [self.row(status="applied")]))

    def test_an_unknown_status(self):
        self.assertIn("status", self.codes([self.row(status="done")]))

    def test_a_settled_row_left_in_the_open_file(self):
        self.assertIn("wrong-file", self.codes([self.row(status="applied")]))

    def test_a_live_row_in_the_closed_file(self):
        self.assertIn("wrong-file", self.codes([], [self.row()]))

    def test_an_unknown_kind(self):
        self.assertIn("kind", self.codes([self.row(kind="idea")]))

    def test_a_finding_the_readme_does_not_list(self):
        self.assertIn("finding", self.codes([self.row(finding="F99")]))

    def test_a_note_with_no_such_heading(self):
        self.assertIn("note", self.codes([self.row(note="FINDINGS.md#f2")]))

    def test_a_note_with_no_such_file(self):
        self.assertIn("note", self.codes([self.row(note="NOPE.md#f1")]))

    def test_a_short_commit(self):
        self.assertIn("pin", self.codes([self.row(checked_against="economicspace@54e30b0")]))

    def test_a_commit_pinned_to_another_repo(self):
        self.assertIn("pin", self.codes([self.row(checked_against="samebytes@" + "a" * 40)]))

    def test_checked_before_it_was_found(self):
        self.assertIn("date", self.codes([self.row(checked_date="2026-09-02")]))

    def test_checked_in_the_future(self):
        self.assertIn("date", self.codes([self.row(checked_date="2026-11-05")]))


if __name__ == "__main__":
    unittest.main()
