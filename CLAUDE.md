# CLAUDE.md: working in this repo

An audit of 17 starred repositories against
[economicspace](https://github.com/loggger101/economicspace), moved out of its
`research/starred-repos/` on 2026-10-20. There is no build and no package:
markdown notes, read-only probe scripts, three helper modules (`orbital.py`,
`snapshots.py`, `run_all.py`), one archival data file and a stdlib `unittest`
suite. [README.md](README.md) says what each file is and how to run the
probes. CI (`.github/workflows/check.yml`) checks only what needs no
economicspace checkout: the scripts compile, the tests pass (every probe's
self-check, and each failing on a planted defect), `sdsstax_ast.tab`'s sha256,
LF line endings, and the notes' relative links.

| file | is the authority for |
|---|---|
| [`README.md`](README.md) | what is here, how to run it, and how far the figures have drifted |
| [`CITATIONS.md`](CITATIONS.md) | every reference, licence and attribution obligation |
| `FINDINGS.md`, `SECOND-PASS.md`, `SPACE-MAP.md`, `SOURCES.md`, `THIRD-PASS.md`, `FOURTH-PASS.md` | the audit as it was measured, with dates |
| [`figures.json`](figures.json) | each recorded headline figure, the release date it belongs to, and the pattern `run_all.py` reads it back with |
| this file | how to edit safely |

## Rules

- The notes are a measurement record. Don't edit a recorded figure to match a
  re-run, because the catalog release under the probes changes and each figure
  belongs to the release it was taken on. Add a new measurement as a new,
  dated section, and note a figure that has moved in README.
- When the notes say "this repo" they mean economicspace, where they were
  written. `modules/`, `verify.py`, `MODEL.md` and the like are its files.
- The probes read only. They run with an economicspace checkout as the working
  directory, never write into its `asteroid_pipeline/`, never run a Stage, and
  never import its modules, because importing `modules/calc.py` runs its
  pip-install block. They mirror its constants instead, with a comment naming
  the source.
- Any probe that prices delta-v reproduces `asteroid_transfer_dv_km_s`'s five
  docstring figures before it measures anything, as the plane-change and NHATS
  probes do. That check is what shows it is pricing the shipped model.
- Copy nothing from the eight copyleft stars (nyx, space-map, CamPyRoS,
  Celestia, and from the third pass ASSIST, REBOUND, REBOUNDx, celmech);
  reimplement from the paper. A GPL tool may be run as a separate step whose
  output a probe reads, never imported by one. Code adapted from a permissive
  project gets its notice, and a note of what changed, in `CITATIONS.md` in
  the same commit.
- `sdsstax_ast.tab` is byte for byte the archive's file (CRLF, `-text` in
  `.gitattributes`). Never rewrite, renormalise or re-download over it.
- Every self-check a probe runs gets a test in `tests/` that runs without a
  catalog, and a planted defect in `tests/test_planted.py` (or beside the
  probe's own tests) that makes it fail. A check that cannot fail proves
  nothing.
- A new headline figure goes into `figures.json` with the date it was taken
  and a pattern that reads it out of the probe's output; check the pattern
  against the recorded output before committing.
- Network reads go through `snapshots.fetch`, so every answer is recorded with
  its sha256 and can be replayed. Downloaded tables go to
  `~/.cache/starred-repos` (outside OneDrive), are pinned by sha256 in the
  probe, and are never committed. Only small excerpts go in
  `tests/fixtures/`, with their provenance in `tests/README.md`.
- Bash heredocs on this machine eat backslashes and `python -` hangs: write a
  script containing `\n` or regexes with the Write tool, then run it.

## Environment

- The checkout lives in OneDrive and is shared with a second machine. Confirm
  a phantom-modified file with `git diff` before acting on it. A
  `*-DESKTOP-PJS73RO.*` file is a sync conflict copy: diff it, delete it, and
  never commit it.
- Text files are LF on every checkout (`.gitattributes`). Python on Windows
  can write CRLF, so check `git diff` after a scripted edit.
