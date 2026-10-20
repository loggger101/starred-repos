# julie-dujardin/space-map, audited on its own

A focused second look at the one star the first pass filed under "thin". That
verdict was wrong. The code is unusable here for licensing reasons, and the
project is a renderer rather than a model, but what it knows about JPL's API is
worth more than most of the permissively-licensed repositories combined.

| | |
|---|---|
| licence | AGPL-3.0 |
| live at | spacemap.co |
| audited | 2026-09-03 |
| shape | `data/` Python ingest pipeline, `frontend/` Svelte + WebGL renderer |

Nothing here is copied, because AGPL would make this pipeline AGPL. Everything
below is an API fact, a data-source pointer, or an idea, and none of those
carries a licence obligation. The one thing taken is a list of field names that
JPL publishes.

---

## 1. The finding: the ranking could not see orbit quality

`data/src/space_map_data/download/providers/objects/sbdb.py` does not hardcode
its field list. It calls the discovery endpoint:

```
https://ssd-api.jpl.nasa.gov/sbdb_query.api?info=field
```

SBDB exposes 79 queryable fields. This pipeline requested 23.

Among the 56 it did not request was `condition_code`, the MPC orbit-uncertainty
parameter U, running 0 (well determined) to 9 (barely constrained). Δv comes
from `a`, `e` and `i`, and the whole economic ranking comes from Δv. So the
model has been quoting cost/revenue ratios to six figures for bodies whose
orbits are provisional, with no way to tell them apart from bodies observed for
70 years.

Measured on the profitability catalog on disk, joined against the full SBDB
table (1,562,105 bodies with a condition code):

| set | n | U >= 5 | vs population | binomial p |
|---|---|---|---|---|
| top 30 by cost/revenue | 30 | 43.3% | 3.1x | 8.3e-05 |
| top 50 | 50 | 30.0% | 2.2x | 2.5e-03 |
| all matched | 157 | 18.5% | 1.3x | 0.065 |
| population | 1.56 M | 13.9% | - | - |

The enrichment is concentrated at the top and disappears down the ranking,
which is what separates a winner's curse from a merely bad population.
Quartiles of the same 157 rows: Q1 35.9%, then 5.1%, 15.4%, 17.5%. Selecting
the extreme of a ranking derived from noisy elements preferentially selects
bodies whose fit errors flatter them.

The bodies at the top, with their orbit arcs:

| designation | ratio | U | arc (days) | obs |
|---|---|---|---|---|
| 2017 MC1 (rank 1) | 18.26 | 7 | 32 | 28 |
| 2015 BN515 (rank 2) | 21.27 | 8 | 12 | 17 |
| 2024 SS4 | 32.67 | 7 | 5 | 33 |
| 2024 XY11 | 36.59 | 8 | 4 | 24 |
| 2015 KJ292 | 37.52 | 9 | 4 | 18 |

Applied in catalog `1.2.0`, as data only: `orbit_condition_code`,
`observation_arc_days`, `n_observations`, `orbit_fit_rms_arcsec`,
`earth_moid_au`, `orbit_class`, `orbit_solution_date`. All seven are 100%
populated.

No filter was applied, deliberately. A U cutoff, an arc-length floor
and a confidence-weighted ranking are all defensible and they answer different
questions. Choosing one changes what the model says, so it is a modelling
decision. Until one of them is chosen, nobody should read the top of the
ranking as a target list.

Coverage was checked before adding anything, and fields with negligible
coverage were deliberately left out: `H_sigma` 0.00%, `G` 0.10%, `BV` 0.85%,
`UB` 0.81%, `IR` 0.00%, `GM` 0.01%, `extent` 0.02%. `BV` and `UB` are colour
indices and would have been relevant to the taxonomy uncertainty in
[FINDINGS.md](FINDINGS.md) F2, but at under 1% they cannot carry it.

---

## 2. Schema discovery beats a hardcoded list

The seven columns matter less than how they came to be missing.

This pipeline's field list is a literal, so it can only ever contain what
somebody thought of. `epoch` was missing for nine releases
([FINDINGS.md](FINDINGS.md) F5); `condition_code` was missing since the
beginning. Both are one call away from being visible:

```bash
curl -s "https://ssd-api.jpl.nasa.gov/sbdb_query.api?info=field" | py -m json.tool
```

Suggested, not applied: a check that fetches the field list and reports
anything SBDB offers that `_JPL_FIELDS` does not request. It costs one HTTP
call, it is the natural home for the coverage numbers above, and it turns "we
never thought of that field" into a line of output. It has the same shape as
`schema_check()` and `stamp_check()`, one level further upstream: those two
verify what arrived, and this one verifies what was *asked for*.

Stage 1's fetch now lives in
[AsteroidCatalog](https://github.com/loggger101/AsteroidCatalog), so the check
belongs there, beside its gates and out of the fetch's hot path.

---

## 3. An incremental mirror, which is the reproducibility answer

Their `SBDBDownloader` keeps `sbdb.sqlite` keyed by `spkid`, pages the full
table once, and from then on fetches only records whose `soln_date` moved. It
reconciles against the remote spkid list to catch deletions (merged
designations) and spkid migrations when an unnumbered body becomes numbered.

That last one is a real hazard. The catalog joins on `designation`, which
*changes* when a body is numbered, while `spkid` is stable and the catalog
already carries it as `spk_id`.

The rest of the case for a mirror has been answered already. Stage 1 now
downloads a pinned, sha256-checked catalog release that AsteroidCatalog builds,
which is the pinned artefact this section argued for. Their operational facts
belong to that repo: page at 5,000 rows with a pause between pages; `soln_date`
filters are date-only, so overlap the window by a couple of days and rely on
upserts; past roughly 200,000 changed records, a full resync beats delta
pagination.

---

## 4. It points at the unquantified error in our own Lambert oracle

`data/scripts/elements_benchmark.py` propagates exported Keplerian elements and
scores them against SPICE truth from `spiceypy.spkezr`, using JPL's
`sb441-n373.bsp`.

That bears directly on [FINDINGS.md](FINDINGS.md) F4. Our
[`probe_lambert.py`](probe_lambert.py) propagates two-body Keplerian orbits
over a time-of-flight grid running to 8 years. Real asteroid orbits precess
under planetary perturbations, so the oracle carries an error nobody has
measured, and that includes the +1.30 km/s headline.

The kernels exist and can be fetched:

| kernel | bodies | size |
|---|---|---|
| `sb441-n16.bsp` | 16 | 646 MB |
| `sb441-n373s.bsp` | 373, short span | 982 MB |
| `sb441-n373.bsp` | 373, full span | 15.2 GB |

at `https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/`.

Reading them needs `spiceypy`, so this is work for an offline validation
harness and never a pipeline dependency. It would bound the two-body assumption
without removing it, and that bound is currently unknown in both directions.

---

## 5. What is not useful here

- The frontend (Svelte and WebGL): spacekit does the same job under MIT.
- Chebyshev ephemeris compression: a renderer needs positions at arbitrary
  times, and this pipeline's Δv model is date-free and needs none.
- `generate_small_body_colors.py`: colour is its output, where taxonomy is
  F2's input, so it adds no classification evidence.
- The probe and spacecraft machinery, which has nothing to do with asteroid
  economics.

---

## Verdict

The most valuable star of the seventeen, and the one whose code is least
usable. Its worth is entirely in what it knows about JPL's API: that the field
list is discoverable, that `condition_code` exists, and that a small-body
catalog is better kept as an incrementally-synced mirror than re-downloaded.

The first pass rated it thin because it judged the code, which AGPL puts out
of reach.
