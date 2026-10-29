# CLAUDE.md: working in this repo

An audit of 17 starred repositories against
[economicspace](https://github.com/loggger101/economicspace), moved out of its
`research/starred-repos/` on 2026-10-20. There is no build, no test suite and
no package: markdown notes, nine read-only probe scripts, one helper module and
one archival data file. [README.md](README.md) says what each file is and how
to run the probes. CI (`.github/workflows/check.yml`) checks only what needs
no economicspace checkout: the scripts compile, `sdsstax_ast.tab`'s sha256,
LF line endings, and the notes' relative links.

| file | is the authority for |
|---|---|
| [`README.md`](README.md) | what is here, how to run it, and how far the figures have drifted |
| [`CITATIONS.md`](CITATIONS.md) | every reference, licence and attribution obligation |
| `FINDINGS.md`, `SECOND-PASS.md`, `SPACE-MAP.md`, `SOURCES.md` | the audit as it was measured, with dates |
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

## Environment

- The checkout lives in OneDrive and is shared with a second machine. Confirm
  a phantom-modified file with `git diff` before acting on it. A
  `*-DESKTOP-PJS73RO.*` file is a sync conflict copy: diff it, delete it, and
  never commit it.
- Text files are LF on every checkout (`.gitattributes`). Python on Windows
  can write CRLF, so check `git diff` after a scripted edit.
