# starred-repos

An audit of 17 starred repositories against the
[economicspace](https://github.com/loggger101/economicspace) asteroid-mining
pipeline, and the code and data worth taking from them. It lived in that repo
as `research/starred-repos/` until 2026-10-20 and moved here with its history;
where these notes say "this repo" or "the pipeline", they mean economicspace.

Nothing here is wired into the pipeline. No module, config or CSV was touched
and no `pipeline_version` moved. Every script reads only, and none of them runs
Stage 1, 2 or 3, so none can invalidate a `.verify` baseline.

## What to read

| file | what it is |
|---|---|
| [`FINDINGS.md`](FINDINGS.md) | Start here. Five findings, each a measurement on economicspace's own catalog, and what became of each |
| [`SOURCES.md`](SOURCES.md) | per-repository inventory: licence, the specific files worth reading, and a verdict. Several stars are worth less here than they look |
| [`CITATIONS.md`](CITATIONS.md) | every reference and attribution obligation in this repo, including what was adapted from whom. The pipeline keeps its own in [economicspace's](https://github.com/loggger101/economicspace/blob/main/CITATIONS.md) |
| [`SPACE-MAP.md`](SPACE-MAP.md) | a re-audit of the one repo the first pass got wrong. It found the orbit-quality gap |
| [`SECOND-PASS.md`](SECOND-PASS.md) | a re-read of all seventeen for what the first pass did not open. Four more findings, including the only check on the delta-v model that economicspace did not build itself |
| [`THIRD-PASS.md`](THIRD-PASS.md) | the 194 stars added after the first audit, read on 2026-10-28. Eighteen touch the domain; prospector, a second asteroid-mining model under Apache-2.0, had code worth adapting. F10 to F14, measured on 2026-10-29, test the shipped delta-v estimator, the guessed taxonomies, the water and spin of the bodies at the top of the ranking, and how well their orbits pin down delta-v |
| [`CLAUDE.md`](CLAUDE.md) | rules for editing here. The notes are a dated record, the probes read only, and nothing comes from a copyleft star |

## What to run

Every script takes the root of an economicspace checkout as its working
directory. It reads that checkout's `asteroid_pipeline/asteroid_catalog.csv`,
and `probe_taxonomy.py` also reads its `catalog_taxonomy.json`; economicspace's
Stage 1 installs both. `probe_perihelion.py`, `probe_taxonomy_prior.py`,
`probe_spin.py` and `probe_covariance.py` also read its
`profitability_catalog.csv`, which a full run writes. The only dependencies
are numpy and pandas ([`requirements.txt`](requirements.txt)).
`probe_nhats.py`, `probe_prospector.py` and `probe_covariance.py` need network
access to public read-only APIs, and `probe_perihelion.py` downloads a 248 MB
table once, to
`~/.cache/starred-repos` rather than this checkout. With this repo cloned
beside economicspace:

```bash
py -m pip install -r ../starred-repos/requirements.txt
py ../starred-repos/probe_plane_change.py
py ../starred-repos/probe_taxonomy.py
py ../starred-repos/probe_lambert.py --n 400 --seed 1
py ../starred-repos/probe_nhats.py --asterank
py ../starred-repos/probe_prospector.py
py ../starred-repos/probe_taxonomy_prior.py
py ../starred-repos/probe_perihelion.py
py ../starred-repos/probe_spin.py
py ../starred-repos/probe_covariance.py
```

| file | what it measures | runtime |
|---|---|---|
| `probe_plane_change.py` | where the plane change is paid (F1). It first reproduces the five figures in `asteroid_transfer_dv_km_s`'s docstring, which shows it is measuring the shipped model | seconds |
| `probe_taxonomy.py` | taxonomy provenance and cross-survey disagreement (F2) | ~30 s |
| `probe_lambert.py` | the closed-form estimator against a real optimised transfer (F4) | ~10 min at n=400 |
| `probe_nhats.py` | the shipped estimator against two oracles other people built: JPL's NHATS and the Shoemaker-Helin closed form (F6). Reproduces the docstring first, and prices both apsides, which the plane-change probe deliberately does not | ~30 s |
| `probe_prospector.py` | prospector's closed form against the shipped one, with NHATS ranking both (F10). Reproduces the docstring and prospector's own output first, then isolates the plane change, the 1 AU crossing and the departure pricing. Imports `probe_nhats.py`; adapted code, see [`CITATIONS.md`](CITATIONS.md) | ~30 s |
| `probe_taxonomy_prior.py` | how often a guessed taxonomy is right (F11): the builder's two rules against fitted multinomial logits on a hold-out, then an empirical prior for NEOs applied to the albedo_assumed NEOs and the top of the ranking. Checks its mirror of the builder's rules against the release, and its fit against synthetic data, first | ~40 s |
| `probe_perihelion.py` | the Toliou et al. (2021) minimum-perihelion prior for every NEO and the top of the ranking (F12). Reproduces prospector's own figures and checks its grid lookup against prospector's full scan first; adapted code, see [`CITATIONS.md`](CITATIONS.md) | ~2.5 min, plus a one-time download and parse |
| `probe_spin.py` | whether the NEOs, and the top of the ranking, spin too fast to hold loose regolith (F13). Prospector's Holsapple cohesion model reduced to a sphere, with its shape-integral bug corrected; reproduces upstream's figures with the bug, and checks the closed form against a bisection, first. Adapted code, see [`CITATIONS.md`](CITATIONS.md) | ~1 min |
| `probe_covariance.py` | NEODyS orbit covariance pushed through the shipped delta-v estimator for the NEOs in the top 100 rows (F14). Reproduces the docstring, and checks each parsed covariance against NEODyS's RMS line, first. One request per body | ~1.5 min |
| `orbital.py` | Kepler solver, elements to state vectors, and a validated Izzo Lambert solver. `probe_lambert.py` imports it; it is not a probe. This is the only copy of the solvers and the one F4 was measured with. economicspace's root `orbital.py` keeps just `elements_to_state`, for the dashboard's orbit diagram, and the two copies of that function are the same text | - |

The figures in the notes were taken on whatever catalog release was installed
at the time, and a re-run against a later release moves them. On 2026-10-20
the plane-change survey drew 39,223 rows where FINDINGS.md's has 38,892, and
the taxonomy probe's `comp_group` agreement read 68.1% against the recorded
66.1%. The reference-body checks, which pin the probes to the shipped model,
still reproduce to under 0.01 km/s.

On 2026-10-26, against catalog release `data-2026-10-21b`, the plane-change
survey drew 39,398 rows and the `comp_group` agreement still read 68.1%. The
reference bodies reproduced to 0.0081 km/s at worst.

On 2026-10-29, on the same release, `probe_prospector.py` joined 7,143 of
NHATS's 7,159 bodies where F6 joined 7,006 of 7,033, and the shipped
estimator's rank agreement with NHATS read 0.5970 against F6's 0.5976.

`py` is the Windows launcher; on Linux or macOS run the same commands with
`python3`.

## Checks

[`.github/workflows/check.yml`](.github/workflows/check.yml) runs on every
push to main and every pull request. It cannot run the probes, which need an
economicspace checkout and its installed catalog. It checks that every script
compiles, that `sdsstax_ast.tab` still has the sha256 recorded in
[`CITATIONS.md`](CITATIONS.md), that every other text file is LF, and that the
notes' relative links and anchors resolve.

## Data

`sdsstax_ast.tab`, 7.5 MB, is the SDSS-based Asteroid Taxonomy V1.1 asteroid
table from the PDS Small Bodies Node. Unlike the Stage 1 catalog it is static
and does not change under a re-fetch, so it is safe to commit.
`probe_taxonomy.py` downloads it if it is missing. Using it obliges you to cite
it; the citation is in [`CITATIONS.md`](CITATIONS.md).

## The short version

The two Δv findings point in opposite directions and partially cancel. F1 says
the model overcharges the plane change; F4 says it undercharges the transfer
geometry by more. Fixing either one alone makes the model worse. If you read
only one thing here, read F4's reversal table.

[F6](SECOND-PASS.md#f6-two-independent-delta-v-oracles-exist-and-neither-was-used)
checks the same thing with the first instruments economicspace did not build
itself. The shipped estimator agrees with the other published closed form at
0.949 and with JPL's optimised trajectories at 0.598. Patched conics as a
family disagree with real trajectories, and the pipeline sits where the
published closed form does.
