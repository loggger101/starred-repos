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
| [`THIRD-PASS.md`](THIRD-PASS.md) | the 194 stars added after the first audit, read on 2026-10-28. Eighteen touch the domain; prospector, a second asteroid-mining model under Apache-2.0, had code worth adapting, and F10 measures it against the shipped delta-v estimator |
| [`FOURTH-PASS.md`](FOURTH-PASS.md) | the measurements that followed, F11 onward, each with its own probe: guessed taxonomies, water, spin, orbit uncertainty, and what came after |
| [`CLAUDE.md`](CLAUDE.md) | rules for editing here. The notes are a dated record, the probes read only, and nothing comes from a copyleft star |

## Findings

Every finding in the notes, one line each. The section linked is the record;
this table only points at it.

| finding | where | what it found | status |
|---|---|---|---|
| F1 | [FINDINGS](FINDINGS.md#f1-two-orbital-elements-are-fetched-stored-and-read-by-nothing) | the whole plane change is charged at departure; splitting it saves a median 0.47 km/s | not applied: F4 measures the fix making the estimator worse |
| F2 | [FINDINGS](FINDINGS.md#f2-composition-rests-on-a-taxonomy-that-is-guessed-for-93-of-the-bodies-that-matter) | 90% of NEO taxonomies are guessed, and two surveys disagree on a third of the bodies both cover | open; F11 measures the guess |
| F3 | [FINDINGS](FINDINGS.md#what-happened-to-each-finding) | pandas 3's dtype inference was unprobed | closed in economicspace's `platform_check.py` |
| F4 | [FINDINGS](FINDINGS.md#f4-the-estimator-is-optimistic-by-11-and-f1-is-the-smaller-of-two-errors-that-cancel) | the closed form sits ~11% below a Lambert optimum, and F1's error partly cancels it | open; a documented limitation in MODEL.md. F19 rules out the oracle's two-body error as the cause |
| F5 | [FINDINGS](FINDINGS.md#f5-mean_anomaly_deg-is-stored-without-its-epoch-so-it-cannot-be-used) | `mean_anomaly_deg` was stored without its epoch | applied in catalog 1.2.0 |
| - | [SPACE-MAP](SPACE-MAP.md#1-the-finding-the-ranking-could-not-see-orbit-quality) | the ranking cannot see orbit quality; poor orbits are enriched 3.1x at its top | open; F14 and F15 say what a poor orbit costs: a lost target, not a wrong delta-v |
| F6 | [SECOND-PASS](SECOND-PASS.md#f6-two-independent-delta-v-oracles-exist-and-neither-was-used) | the estimator ranks with JPL's NHATS at 0.60 and with Shoemaker-Helin at 0.95; 14 one-way legs exceed JPL's round trip | measured |
| F7 | [SECOND-PASS](SECOND-PASS.md#f7-brahe-is-not-thin-any-more-and-it-closes-space-mapmds-open-item) | brahe carries a SPICE reader and a Horizons SPK client | capability; F19 bounds F4 another way |
| F8 | [SECOND-PASS](SECOND-PASS.md#f8-space-datasets-was-read-at-2-per-cent) | space-datasets carries cited sources for the softest assumptions | inventory |
| F9 | [SECOND-PASS](SECOND-PASS.md#f9-astroquery-carries-an-orbit-covariance-client) | NEODyS serves full orbit covariance | used by F14 and F15 |
| F10 | [THIRD-PASS](THIRD-PASS.md#f10-prospectors-estimator-ranks-closer-to-jpl-and-the-plane-change-is-not-why) | moving the plane change to arrival is never cheaper on NHATS NEAs; prospector's better rank comes from its departure pricing | measured |
| F11 | [FOURTH-PASS](FOURTH-PASS.md#f11-a-guessed-taxonomy-is-right-about-half-the-time-and-a-guessed-c-less) | a guessed class is right about half the time; a guessed C 29% | measured |
| F12 | [FOURTH-PASS](FOURTH-PASS.md#f12-the-top-of-the-ranking-probably-never-passed-close-to-the-sun) | the water-carrying NEOs at the top almost never passed inside 0.2 AU | measured |
| F13 | [FOURTH-PASS](FOURTH-PASS.md#f13-the-bodies-that-lead-the-ranking-spin-too-fast-to-hold-loose-regolith) | a third of NEOs with a measured spin, and all past H 28, cannot hold loose regolith | measured; prospector's cohesion model corrected; F18 puts the share at 25.6% on reliable periods; economicspace calc 1.34.0 adds an opt-in limit |
| F14 | [FOURTH-PASS](FOURTH-PASS.md#f14-orbit-uncertainty-barely-moves-delta-v-it-moves-where-the-body-is) | orbit uncertainty barely moves delta-v for the top bodies; it moves where they are | measured |
| F15 | [FOURTH-PASS](FOURTH-PASS.md#f15-most-of-the-top-of-the-ranking-cannot-be-found-again-before-2036) | 42 of the 69 top-100 NEOs are never observable before 2036, and every body with U 7 to 9 is unobservable or lost | measured |
| F16 | [FOURTH-PASS](FOURTH-PASS.md#f16-the-guessed-neos-have-no-measured-albedo-anywhere-the-catalog-missed) | no albedo-guessed NEO has a measured albedo in NEOWISE, SsODNet's datacloud or the LCDB | measured; SsODNet name-matching trap recorded |
| F17 | [FOURTH-PASS](FOURTH-PASS.md#f17-a-debiased-population-model-agrees-the-guessed-cs-are-mostly-not-dark) | Granvik et al.'s debiased model disagrees with 2,959 of the 3,460 guessed Cs | measured |
| F18 | [FOURTH-PASS](FOURTH-PASS.md#f18-the-lightcurve-database-adds-no-spins-and-trims-f13s-headline) | the LCDB adds no NEO periods; on its reliable ones F13's share is 25.6%, and the faint bins hold | measured |
| F19 | [FOURTH-PASS](FOURTH-PASS.md#f19-f4s-oracle-is-not-the-problem-its-two-body-error-is-00035-kms) | F4's oracle carries a median 0.0035 km/s of two-body error against a 1.3 km/s gap | measured; closes SPACE-MAP section 4 |

## What to run

Every script takes the root of an economicspace checkout as its working
directory. It reads that checkout's `asteroid_pipeline/asteroid_catalog.csv`,
and `probe_taxonomy.py` also reads its `catalog_taxonomy.json`; economicspace's
Stage 1 installs both. The probes that look at the top of the ranking also
read its `profitability_catalog.csv`, which a full run writes. The only
dependencies are numpy and pandas ([`requirements.txt`](requirements.txt)).

Network and downloads:
- **Public read-only APIs.** These probes call them:
  - `probe_nhats.py` and `probe_prospector.py`: NHATS and Asterank
  - `probe_covariance.py` and `probe_recovery.py`: NEODyS
  - `probe_albedo_coverage.py`: IRSA and SsODNet
  - `probe_horizons.py`: JPL Horizons
- **Recorded responses.** Every response goes through
  [`snapshots.py`](snapshots.py), which saves it under
  `~/.cache/starred-repos/snapshots/<date>/` with its sha256. Setting
  `STARRED_REPOS_REPLAY=<that folder>` (or passing `run_all.py --replay`)
  re-runs a probe from the recording instead of the network.
- **Static tables.** Three probes download one once each, to
  `~/.cache/starred-repos` rather than this checkout, and pin it by sha256:
  - `probe_perihelion.py`: Toliou, 248 MB
  - `probe_source_regions.py`: Granvik, 33 MB
  - `probe_lightcurves.py`: the LCDB, 40 MB

With this repo cloned beside economicspace:

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
py ../starred-repos/probe_recovery.py
py ../starred-repos/probe_albedo_coverage.py
py ../starred-repos/probe_source_regions.py
py ../starred-repos/probe_lightcurves.py
py ../starred-repos/probe_horizons.py
```

| file | what it measures | runtime |
|---|---|---|
| `probe_plane_change.py` | where the plane change is paid (F1). It first reproduces the five figures in `asteroid_transfer_dv_km_s`'s docstring, which shows it is measuring the shipped model | seconds |
| `probe_taxonomy.py` | taxonomy provenance and cross-survey disagreement (F2) | ~30 s |
| `probe_lambert.py` | the closed-form estimator against a real optimised transfer (F4). Since 2026-10-29 it reproduces the docstring and checks its Lambert solver on four round trips before measuring | ~3 to 10 min at n=400 (164 s on 2026-10-29) |
| `probe_nhats.py` | the shipped estimator against two oracles other people built: JPL's NHATS and the Shoemaker-Helin closed form (F6). Reproduces the docstring first, and prices both apsides, which the plane-change probe deliberately does not | ~30 s |
| `probe_prospector.py` | prospector's closed form against the shipped one, with NHATS ranking both (F10). Reproduces the docstring and prospector's own output first, then isolates the plane change, the 1 AU crossing and the departure pricing. Imports `probe_nhats.py`; adapted code, see [`CITATIONS.md`](CITATIONS.md) | ~30 s |
| `probe_taxonomy_prior.py` | how often a guessed taxonomy is right (F11): the builder's two rules against fitted multinomial logits on a hold-out, then an empirical prior for NEOs applied to the albedo_assumed NEOs and the top of the ranking. Checks its mirror of the builder's rules against the release, and its fit against synthetic data, first | ~40 s |
| `probe_perihelion.py` | the Toliou et al. (2021) minimum-perihelion prior for every NEO and the top of the ranking (F12). Reproduces prospector's own figures and checks its grid lookup against prospector's full scan first; adapted code, see [`CITATIONS.md`](CITATIONS.md) | ~2.5 min, plus a one-time download and parse |
| `probe_spin.py` | whether the NEOs, and the top of the ranking, spin too fast to hold loose regolith (F13). Prospector's Holsapple cohesion model reduced to a sphere, with its shape-integral bug corrected; reproduces upstream's figures with the bug, and checks the closed form against a bisection, first. Adapted code, see [`CITATIONS.md`](CITATIONS.md) | ~1 min |
| `probe_covariance.py` | NEODyS orbit covariance pushed through the shipped delta-v estimator for the NEOs in the top 100 rows (F14). Reproduces the docstring, and checks each parsed covariance against NEODyS's RMS line, first. One request per body | ~1.5 min |
| `probe_recovery.py` | whether the NEOs in the top 100 rows can be observed again before 2036, and how large a patch of sky would have to be searched (F15). Checks its Earth against the almanac and its covariance propagation against draws first | 35 s replayed; live adds one NEODyS request per body |
| `probe_albedo_coverage.py` | whether the albedo-guessed NEOs have a measured albedo in NEOWISE V2.0 or SsODNet's datacloud (F16). Mirrors the builder's rules, and proves its NEOWISE join on bodies the catalog already matched, first | about 40 s replayed; live adds about 660 SsODNet requests |
| `probe_source_regions.py` | Granvik et al.'s (2018) debiased source-region and albedo model against the guessed classes (F17). Checks the model file and the grid lookup first | ~2 min, plus a one-time download and parse |
| `probe_lightcurves.py` | the Asteroid Lightcurve Database against the catalog's periods and F13 (F18). Reproduces four known rows of the table first | ~20 s, plus a one-time download |
| `probe_horizons.py` | the two-body error in F4's Lambert oracle, measured by running it on JPL Horizons' trajectories (F19). Reproduces the docstring and solver checks, and checks its interpolation and frame, first | ~3 min live |
| `run_all.py` | not a probe: runs the probes and compares their headline figures with [`figures.json`](figures.json); see [Checks](#checks) | the sum of the probes |
| `snapshots.py` | not a probe: records and replays every network response | - |
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

Also on 2026-10-29, the first `run_all.py` run found three more figures moved.
- **F4's median gap** (`probe_lambert.py` at n=400, seed 1): +1.367 km/s
  against the recorded +1.296. The sample is drawn from today's catalog, so it
  is not F4's 400 bodies.
- **F1's survey**: 39,398 rows (9.901 km/s at the median; saving 0.4687),
  against 38,892 rows (9.904; 0.4692).
- **F2's agreement**: still 68.1%, against the recorded 66.1%.

Every F10 to F19 figure reproduced exactly.

`py` is the Windows launcher; on Linux or macOS run the same commands with
`python3`.

## Checks

[`.github/workflows/check.yml`](.github/workflows/check.yml) runs on every
push to main and every pull request. It cannot run the probes, which need an
economicspace checkout and its installed catalog. It checks that every script
compiles, that `sdsstax_ast.tab` still has the sha256 recorded in
[`CITATIONS.md`](CITATIONS.md), that every other text file is LF, and that the
notes' relative links and anchors resolve.

It also runs [`tests/`](tests/README.md): every probe's self-check (the
docstring figures, prospector's own output, a synthetic recovery, a parser
against a fixture), and every one of those checks failing on a planted
defect. They need no catalog. Run them locally with
`python -m unittest discover -s tests -t . -v`.

[`run_all.py`](run_all.py) is the other half. Run from here with
economicspace beside it, it runs every probe against the installed release,
reads each headline figure back out with the patterns in
[`figures.json`](figures.json), and prints it beside the recorded value with
the release it was taken on. It writes the run to a dated JSON file under
`~/.cache/starred-repos/runs/` and exits non-zero if a self-check fails or a
figure goes missing. `--slow` adds the Lambert probe, `--offline` leaves out
the probes that call public APIs, and `--replay <folder>` runs them from a
recorded snapshot instead (see [`snapshots.py`](snapshots.py)). Replaying the
2026-10-29 recording reproduced the live NHATS and NEODyS outputs line for
line.

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
