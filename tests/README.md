# tests

Every probe checks itself before it measures. It reproduces a figure from the
shipped model or from upstream code, or recovers a known answer, and it stops
if it cannot. Those checks need no economicspace checkout, so they live here
too and CI runs them on every push.

```bash
python -m unittest discover -s tests -t . -v
```

| file | what it holds |
|---|---|
| `test_checks.py` | each probe's self-check for F1 to F15, plus the small pure functions behind them (parsers, lookups, closed forms) |
| `test_planted.py` | each of those self-checks failing on a planted defect: a recorded figure moved, a constant changed, an exponent restored, a parse order reversed. A check that cannot fail proves nothing |
| `test_probes_f16_f18.py` | F16's identity filter and NEOWISE join, F17's model check, F18's table parser, each with its planted defect |
| `test_probes_f15_f19.py` | F19's Horizons parser, Hermite interpolation and its fourth-power check, with a planted linear interpolant |
| `test_snapshots.py` | `snapshots.py` recording, replaying, and refusing a tampered or unrecorded response |
| `_util.py` | the repo root on `sys.path`, quiet calls, fixture loaders, a synthetic Toliou-shaped grid |

Tests that need a local download skip when it is absent: the Toliou table
(`probe_perihelion.py`) and the Granvik model (`probe_source_regions.py`),
both cached under `~/.cache/starred-repos`. Run locally with those caches
present and nothing skips.

## Fixtures

All are small text samples, committed LF (`.gitattributes`); the CSV is
exempt from the repo's `*.csv` ignore rule.

`fixtures/builder_rows.csv` is 1,416 rows of catalog release
`data-2026-10-21b` (AsteroidCatalog, MIT): up to 60 rows of each
`spectral_type_source` in {`albedo`, `albedo_assumed`} for each orbit bin and
NEO flag the builder distinguishes, within 4.6 AU, sampled with seed 1. It
lets `check_builder` run without the catalog. On the full release the same
check covers all 1,381,412 inferred rows.

`fixtures/neodys_2000SG344.ke1` is NEODyS's Keplerian file for 2000 SG344 as
served on 2026-10-29 (University of Pisa and SpaceDyS), line endings converted
to LF. It is one of the fourteen bodies F6 flags, and it exercises the
covariance parser, the RMS check and F15's covariance propagation.

`fixtures/lcdb_sample.txt` is the header and four rows (Ceres, Itokawa,
Bennu, 2008 TC3) of the Asteroid Lightcurve Database's `lc_summary_pub.txt`,
release of 2023-10-01 (Warner, Harris and Pravec; see CITATIONS.md). It is
what `check_parse` reproduces.

`fixtures/horizons_vectors.txt` is JPL Horizons' answer for 2000 SG344's
heliocentric ecliptic state vectors, 2027-01-01 to 2027-03-01 every two days,
as served on 2026-10-29, line endings converted to LF.
