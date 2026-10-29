# Third pass: the stars added since the first audit

[`FINDINGS.md`](FINDINGS.md) and [`SECOND-PASS.md`](SECOND-PASS.md) read the
same seventeen repositories twice. This pass reads the ones starred after
them. On 2026-10-28 `loggger101` had 213 stars, 194 of them added after the
first audit on 2026-09-03.

Most of the 194 have nothing to say about an asteroid-mining model. Five are
the user's own repositories (economicspace, spacecost, AsteroidCatalog,
asteroid-belt-gradient, General_Research), which are not audit targets. About
170 are agent skills, web and UI frameworks, linters, curated lists and
general-purpose libraries (numpy, scipy, pandas, pytest and the like); the
first pass set agent tooling aside, and this one sets the rest of that group
aside with it. The 18 that touch the domain were cloned and read against the
pipeline, and they are the subject of this note.

One of them is the most useful star since space-datasets.
Karmanplus/prospector is a second asteroid-mining model, Apache-2.0, and parts
of it can be adapted here.

Nothing here is wired into the pipeline. The pass itself, on 2026-10-28, was a
reading plus the live checks dated below. One measurement followed on
2026-10-29, with a probe:
[F10](#f10-prospectors-estimator-ranks-closer-to-jpl-and-the-plane-change-is-not-why).

| repository | licence | verdict |
|---|---|---|
| [Karmanplus/prospector](#karmanplusprospector-apache-20) | Apache-2.0 | tier 1. The closest prior art this project has, with code that can be adapted |
| [esa/pykep](#esapykep-mpl-20) | MPL-2.0 | tier 1, as `SOURCES.md` predicted. Not on PyPI for Windows |
| [matthewholman/assist](#matthewholmanassist-gpl-30) | GPL-3.0 | tier 1 as an instrument, never as code. It can bound F4's two-body error on near-Earth bodies |
| [statsmodels/statsmodels](#tier-2-worth-one-specific-narrow-thing) | BSD-3 | tier 2. F2's prior is a multinomial logit, which needs no sampler |
| [nyx-space/anise](#tier-2-worth-one-specific-narrow-thing) | MPL-2.0 | tier 2. A second SPICE reader with a Windows wheel |
| [relf/EGObox](#tier-2-worth-one-specific-narrow-thing) | Apache-2.0 | tier 2. Sample-efficient calibration, where pygmo2 was proposed |
| [visgl/deck.gl](#tier-2-worth-one-specific-narrow-thing) | MIT | tier 2. Millions of points in Streamlit, through pydeck |
| [AndrewAnnex/SpiceyPy](#tier-2-worth-one-specific-narrow-thing) | MIT | tier 2. The reference to check brahe and ANISE against |
| the other ten | | [tier 3](#tier-3-thin-for-this-project-and-why) |

---

## Karmanplus/prospector, Apache-2.0

Prospector screens the small-body catalogue against a low-thrust vehicle's
delta-v budget, ranks what survives by what would make it worth mining, and
hands the best candidates to progressively more expensive trajectory solvers.
It is the same problem economicspace solves, approached from the trajectory
end rather than the economics end. Version 0.4.3, last commit 2026-09-15.

Asterank was the only prior art `SECOND-PASS.md` found, and its economics come
back degenerate through the public API. Prospector's physics is validated:
[`docs/validation.md`](https://github.com/karmanplus/prospector/blob/main/docs/validation.md)
flies Dawn, Psyche, Hayabusa2 and DART with their real engines and masses. The
search is never told where the missions went, and it finds Dawn's Mars flyby
four days and 42 km from the flown one and Hayabusa2's Earth swing-by three
days and 560 km from it. Dawn's propellant comes out at 224 kg against 247 kg
flown, and Hayabusa2's at 22 kg against 24 kg.

### It prices the plane change where F1 says to

`prospector/solvers/edelbaum.py`, 151 lines of numpy, carries two closed forms.
`intercept_dv` is the impulsive one, and it is built the way
`asteroid_transfer_options_km_s` is: a Hohmann half-ellipse from 1 AU out to an
apsis of the target, and one burn to match the target there. Somebody else
arrived at the same construction independently.

It differs in two places, and both are live questions here.

1. It charges the inclination on the match burn, at the target's local speed.
   economicspace charges it on departure, in the `v_inf` law of cosines in
   `_transfer_legs_for_apsis`. In `probe_plane_change.py`'s terms, the shipped
   model is placement A and prospector is placement B. F1's fix is placement
   C, the best split, which has A and B as its endpoints, and F4 measured it
   making the estimator worse in every inclination band.
2. It tries a third meeting point, the 1 AU crossing, for orbits with
   q <= 1 <= Q. Its own documentation says the crossing usually loses, because
   the target carries a large radial velocity there.

That makes it a third closed-form oracle beside Shoemaker-Helin, and the only
one that takes the other side of F1. Put both estimators against NHATS, as
`probe_nhats.py` does for the shipped one, and JPL referees F1 and F4 directly.
Two cautions carry over from F6. Prospector departs from Earth's heliocentric
orbit, with no escape from LEO, so only ranks compare and not absolute figures.
And a rank against NHATS measures agreement with a round trip, not the one-way
leg.

`spiral_dv` and `lowthrust_dv` are the low-thrust half: an Edelbaum spiral with
an eccentricity term, which runs high, bracketed against `intercept_dv`, which
runs low. That bracket bears on `dv_penalty_factor`, the scalar the pipeline
uses for "the rocket equation ignores thrust; trajectories don't".

### The low-thrust solve that would calibrate `dv_penalty_factor`

`SOURCES.md` recommended pykep's Sims-Flanagan leg for calibrating
`dv_penalty_factor` on a few hundred sampled bodies. `prospector/solvers/simsflanagan.py`
is that, built and validated. It wraps pykep 3's `sf_pl2pl` with analytic
gradients for its own constraints (arrival deadline, per-segment thrust
ceilings from the solar arrays, departure cone), which took a reference solve
from 39.6 s to 3.3 s. Its documentation measures a finite-difference gradient
as roughly 50% wrong on the departure-epoch column. Anyone repeating this with
pykep alone will hit that.

### Where inclination is bought

`docs/physics.md`, section "Where inclination is bought", is the clearest
statement of F1's question I have seen. A tilt can be bought in the departure
asymptote, at about 1.9 degrees per km/s of out-of-ecliptic v-infinity; in the
geocentric spiral, nearly free in propellant and paid in time; or in the
heliocentric cruise, at about 0.52 km/s per degree at 1 AU. economicspace has
only the first.

### Enrichment: evidence on water that does not come from albedo

F2's problem is that composition is guessed from albedo for most NEOs.
Prospector's `enrichment/` package carries three independent lines of evidence,
and two of them are about water, the resource that matters most in space.

| file | what it is | relevance here |
|---|---|---|
| `enrichment/models/perihelion.py` | reads the Toliou, Granvik and Tsirvoulis (2021) grid: for a body's (a, e, i, H), the probability that its population passed inside each of 26 perihelion thresholds. Below about 0.2 AU, heating strips hydrated minerals | a published prior on whether a C-type NEO still holds its water, independent of the spectrum. The table is live, 247,823,013 bytes gzipped (checked 2026-10-28) |
| `enrichment/models/hydration.py` | a small network trained on 630 RELAB meteorite spectra, mapping albedo, the 0.7 um band depth and four Mahlke taxonomy scores to a petrologic hydration class | needs spectra, so it covers few bodies. The weights ship as a plain `.npz` evaluated in numpy; credit RELAB if used |
| `enrichment/models/cohesion.py` | Holsapple (2007) rotational-disruption strength, the cohesion a rubble pile needs to survive its measured spin | the anchoring question F8 raised against the lightcurve database |
| `enrichment/sources/astorb.py` | Lowell Observatory's AstorbDB GraphQL API: every survey determination of albedo, taxonomy and lightcurve per body, with tumbling and quality flags | a third taxonomy opinion for F2. It returned `502 Bad Gateway` when checked on 2026-10-28 |

`enrichment/models/sizing.py` uses one median albedo per spectral complex
(Marsset et al.), the same flat assumption F2 criticises, so it adds nothing.

### What can be taken

Apache-2.0 permits adapting any of it. The obligations are to carry the
`LICENSE` and `NOTICE` files with what is taken and to mark what was changed.
Under this repo's rules that notice goes into `CITATIONS.md` section 3 in the
same commit as the code. Prospector itself is installed with pixi and conda,
because pykep is not on PyPI for Windows, so take the numpy files rather than
the package. `edelbaum.py` and `perihelion.py` import only numpy, pandas and a
constants module.

---

## esa/pykep, MPL-2.0

`SOURCES.md` called this "the star you are missing". It was starred on
2026-09-23, and the recommendation stands: `lambert_problem`, Keplerian
planets defined by orbital elements, and `trajopt.sf_pl2pl`, with
`pl2pl_N_impulses`, `primer_vector` and `mga_1dsm` beside them. It ships
`de440s.bsp` in `pykep/data`.

Two facts the first pass did not have. Prospector's README says version 3 is
the first under MPL-2.0 and earlier releases were GPL, so pin `>=3`; `COPYING`
on the default branch is MPL-2.0 (last commit 2026-10-09). And
`pip download --only-binary=:all: pykep` finds no distribution at all for
Python 3.14 on Windows (checked 2026-10-28). Prospector's README says PyPI
carries Linux wheels only, with no source fallback, and building needs heyoka,
LLVM and Boost. On this machine it means conda-forge or WSL.

MPL-2.0 is file-level copyleft. Calling pykep from an MIT probe is fine. A
pykep source file copied in stays MPL, so read it, don't vendor it.

---

## matthewholman/assist, GPL-3.0

ASSIST integrates test particles in the field of the Sun, Moon, planets and 16
massive asteroids, with positions from DE441 and its asteroid perturber file,
plus harmonics, relativity and Marsden non-gravitational terms. The README
says every term matches JPL's small-body integrator to machine precision. It
is an extension of REBOUND (hannorein/rebound, also GPL-3.0, also starred).

It is the instrument for `SECOND-PASS.md`'s item 5 and `SPACE-MAP.md` section
4: bounding the error that two-body propagation puts into F4's Lambert oracle.
The route F7 recommended reads truth from the `sb441` kernels through brahe,
and those cover 373 bodies at most, chosen as perturbers for their mass, which
makes them main-belt bodies almost without exception. ASSIST can make truth for any body from its elements, including the
near-Earth ones that planetary encounters move most. The two are
complementary. Use brahe and the kernels to check ASSIST on bodies both cover,
and ASSIST to reach the rest.

The licence decides the shape. Nothing is copied. A probe that imports a GPL
library is a combined work under the FSF's reading, and this repository is
MIT, so keep the GPL step separate: a script run on its own writes propagated
states to a CSV, and an MIT probe reads the CSV. `pip download` finds a Windows
wheel for REBOUND (5.2.2) and none for ASSIST, so ASSIST needs a C compiler or
WSL. It also needs `de440.bsp` and `sb441-n16.bsp`, about 760 MB together.

---

## Tier 2: worth one specific narrow thing

| repo | licence | the one thing | verdict |
|---|---|---|---|
| statsmodels/statsmodels | BSD-3 | `MNLogit` | F2 proposes "a categorical prior over `spectral_type` conditioned on albedo and orbital class, calibrated on the 44,061 bodies with two independent classifications". That model is a multinomial logistic regression, and statsmodels fits it by maximum likelihood with no sampler, giving calibrated class probabilities per body. Fit it before reaching for PyMC, which earns its place only once the Nesvorny families (F8) give the prior a hierarchy |
| nyx-space/anise | MPL-2.0 | a pure-Rust SPICE replacement with Python bindings | `anise 0.10.6` has a Python 3.14 Windows wheel on PyPI (checked 2026-10-28). A second reader for F7's kernels if brahe's falls short. Same author as nyx, but MPL rather than AGPL |
| AndrewAnnex/SpiceyPy | MIT | `spkezr` over CSPICE | what `SPACE-MAP.md` section 4 named, and superseded for this job by brahe (F7). Keep it as the reference implementation brahe and ANISE are checked against |
| relf/EGObox | Apache-2.0 | efficient global optimisation with Gaussian-process surrogates | `SOURCES.md` gave pygmo2 one real use, calibration: which `IN_SPACE_UTILITY_BY_DESTINATION` vector flips a destination. Each evaluation is a pipeline run, so a sample-efficient surrogate optimiser suits that better than pygmo's population methods. `egobox 1.0.0` has a Windows wheel. It gives no admissible bound either, so it cannot replace the exhaustive ladder |
| visgl/deck.gl | MIT | WebGL point layers that draw millions of rows | `ui.py` runs on Streamlit, which ships pydeck as `st.pydeck_chart`. With positions solved in numpy at one epoch, a Cartesian scatter layer would draw the whole catalog, where spacekit (tier 2 in `SOURCES.md`) needs a custom shader. It is a display change only, and it moves no number |

---

## Tier 3: thin for this project, and why

| repo | licence | verdict |
|---|---|---|
| hannorein/rebound | GPL-3.0 | needed only as ASSIST's engine |
| dtamayo/reboundx | GPL-3.0 | extra forces for REBOUND, Yarkovsky among them. Yarkovsky drift matters to NEO ephemerides over decades and does not matter to a mission-cost ranking |
| shadden/celmech | GPL-3.0 | secular and resonant perturbation theory. Nothing the ranking uses |
| CS-SI/Orekit | Apache-2.0 | capable of interplanetary propagation with JPL ephemerides, but it is Java, and brahe and ANISE cover the same ground in Python |
| nyx-space/hifitime | MPL-2.0 | time scales to the nanosecond. The probes take TDB Julian dates from SBDB, and nothing they price is sensitive to a leap second |
| astropy/astropy | BSD-3 | already present as astroquery's dependency. The probes mirror their constants from economicspace on purpose, and astropy's would differ in the last digits |
| jobovy/galpy | BSD-3 | galactic dynamics, the wrong scale entirely |
| HIPS/autograd | MIT | gradients of numpy code, but `calc.py` is scalar `math` Python. A sensitivity study would need a rewrite to use it, and a finite difference costs only a re-run |
| sympy/sympy | BSD-3 | could check the docstring algebra symbolically. Z3 already covers the proof the open branch-and-bound item needs |
| ranaroussi/yfinance | Apache-2.0 | platinum and palladium futures, but economicspace installs pinned, cited prices from a spacevalue release, and Yahoo's terms limit the data to personal use. Not a source a price release could cite |

Toblerity/rtree, Cyan4973/xxHash, iliekturtles/uom, statrs-dev/statrs,
sparsemat/sprs, rust-num/num-traits and GSA/data.gov were looked at and carry
nothing specific to small bodies.

---

## F10. Prospector's estimator ranks closer to JPL, and the plane change is not why

Measured on 2026-10-29 with [`probe_prospector.py`](probe_prospector.py),
against catalog release `data-2026-10-21b` and NHATS as fetched that day. The
probe prices four one-way estimators over every NHATS body in the catalog:

| estimator | what it is |
|---|---|
| `shipped` | placement A at both apsides, departure from LEO: economicspace's model |
| `prospector` | `intercept_dv` as upstream wrote it: placement B, three meeting points, departure priced heliocentrically with no escape from LEO |
| `B_leo` | placement B at both apsides, departure from LEO, on the shipped constants. Against `shipped`, only the plane change differs |
| `B_leo_cross` | `B_leo` plus the 1 AU crossing. Against `B_leo`, only the crossing differs |

Three checks ran before anything was measured. The shipped copy reproduces the
five docstring figures to 0.0081 km/s at worst. The adapted `intercept_dv`
reproduces prospector 0.4.3's own output on the same five bodies to
2.68e-07 km/s. At zero inclination `B_leo` equals `shipped` to 1.92e-13 km/s,
so placement is the only difference between them. Each of the last two checks
was shown to fail on a planted defect: a reference figure moved by 1e-4 km/s,
and `B_leo` given a different Earth speed from `shipped`.

`py ../starred-repos/probe_prospector.py`:

```
NHATS  : 7,159 NEAs from JPL
joined : 7,143  (99.8% of NHATS)

  estimator      median      rho    i < 5 deg     5-15 deg   i > 15 deg   > round trip
                                      n=3,754      n=3,346         n=43
  shipped         4.663   0.5970       0.6043       0.5748       0.5611             14
  prospector      4.976   0.7593       0.9066       0.6888       0.4425              -
  B_leo           5.800   0.3157       0.3553       0.1538      -0.3475            614
  B_leo_cross     5.096   0.6065       0.6317       0.5347      -0.0129            101

  meeting point chosen by B_leo_cross:
    crossing     3,875  (54.2%)
    aphelion     3,100  (43.4%)
    perihelion     168  (2.4%)
  where the crossing wins it saves median 1.317 km/s, max 5.372

  B_leo minus shipped : median +1.145  mean +1.376  p10 +0.151  p90 +2.978 km/s
  B_leo cheaper than shipped on 0 of 7,143 bodies

DIAGNOSIS
  departure priced heliocentrically, apsides only:
    placement A  0.7399
    placement B  0.7591   (prospector without the crossing)
  shipped from LEO, match burn weighted by w:
    w = 0.25   0.7394
    w = 0.50   0.7150
    w = 0.75   0.6588
    w = 1.00   0.5970
    w = 2.00   0.4108
```

`rho` is the Spearman rank agreement with NHATS, and the last column counts
one-way legs that cost more than JPL's round trip. The shipped `rho` was
0.5976 in F6 and is 0.5970 now, on a later catalog and NHATS release.

### On the NEAs JPL can reach, the plane change belongs at departure

`B_leo` is cheaper than `shipped` on none of the 7,143 bodies. Moving the whole
plane change to the match burn adds a median 1.145 km/s, puts 614 one-way legs
above JPL's round trip where the shipped model puts 14, and drops the rank
agreement from 0.5970 to 0.3157.

The likely reason is legible. An NHATS body is on an Earth-like orbit, so at
the meeting point it moves at nearly Earth's speed and a tilt costs about as
much there as at departure. The departure burn is fired deep in Earth's
gravity well from LEO, where the Oberth effect discounts it, and the match burn
gets no such discount. F1's survey found placement B cheaper than A at the
median (9.610 against 9.904 km/s), but that was a 1-in-40 stride of the whole
catalog, which is mostly main belt, where the meeting is slow. For accessible
NEAs, F1's endpoint B has no saving to offer. Whether its best split C saves
anything here was not measured; C cannot be worse than A, and on this evidence
it will sit close to A.

### Prospector's advantage is the heliocentric departure

Prospector ranks at 0.7593 against the shipped 0.5970, and it differs in three
ways. The diagnosis takes them apart. Placement B with a LEO departure ranks
worse, at 0.3157. The 1 AU crossing matters only under a LEO departure, where
it pulls B back to 0.6065, and with a heliocentric departure it adds 0.0002
(0.7591 to 0.7593). Pricing the departure heliocentrically lifts both
placements, A to 0.7399 and B to 0.7591, so that pricing is what carries the
gain.

NHATS being a round trip does not explain it. A round trip pays roughly the
asteroid-side burn twice, and weighting the match burn by 2 drops the shipped
estimator to 0.4108. The match burn is already weighted too heavily.
Down-weighting it to a quarter lifts the shipped model to 0.7394, the same as
a heliocentric departure does. Both say the same thing: across these bodies,
the closed form's apsis match varies more relative to the departure than JPL's
optimised trajectories do. JPL meets each body wherever and whenever is
cheapest, and the closed form only ever meets it at an apsis. That is the
geometry F4 says the estimator undercharges, so F4's undercharge is the natural
suspect. This probe cannot confirm it.

### Cautions

These are ranks against a round trip, so they say nothing about the size of
any estimator's error. Only F4's oracle sizes it. The i > 15 deg band holds 43
bodies and its figures are noise. NHATS is a selected population, the NEAs
that a crewed mission could reach, and the top of the profitability ranking is
mostly not in it. And the weight `w` is a diagnostic. Fitting it would add a
second scalar fudge beside `dv_penalty_factor` to cover a geometry error, and
nothing here proposes that.

No module, config or CSV was touched, and no number a run produces moved.

---

## What I would do, in order

1. ~~Adapt prospector's `intercept_dv` into a probe, with its notice in
   `CITATIONS.md`. Reproduce `asteroid_transfer_dv_km_s`'s five docstring
   figures first, as every pricing probe here does, then put the shipped
   estimator, prospector's and NHATS side by side.~~ Done in F10. On
   accessible NEAs it retires placement B, and it points the shipped
   estimator's rank error at the match leg.
2. Bound F4's oracle with ASSIST on a sample of near-Earth bodies, as a
   separate GPL step writing a CSV. This replaces `SECOND-PASS.md`'s item 5,
   whose kernels reach only the main belt.
3. Fit F2's prior with `MNLogit` before building anything with PyMC.
4. Read the Toliou table against the C-complex NEOs near the top of the
   ranking. Water value there rests on the taxonomy alone, and the table is a
   second, independent reason to believe or doubt it.
5. Calibrate `dv_penalty_factor` against Sims-Flanagan, following
   prospector's solver. Last, because it needs pykep, and on Windows that means
   conda or WSL.

None of these moves a number a run produces.

## Licence position

prospector, EGObox and Orekit are Apache-2.0; statsmodels, astropy and galpy
BSD-3; SpiceyPy, deck.gl and autograd MIT; pykep, ANISE and hifitime MPL-2.0.
ASSIST, REBOUND, REBOUNDx and celmech are GPL-3.0, which brings the copyleft
stars to eight. Nothing was copied from any of them. The Toliou table is
published research data from the University of Helsinki, and RELAB is NASA PDS
data. [`CITATIONS.md`](CITATIONS.md) is the authority.
