# Fourth pass: measuring with what the third pass found

[`THIRD-PASS.md`](THIRD-PASS.md) read the 194 stars added after the first
audit and measured one of them, prospector, against the shipped delta-v
estimator (F10). This note holds the measurements that followed, each with
its own probe. All were taken on 2026-10-29. F11 to F14 were moved here from
THIRD-PASS.md unchanged when it grew too long to read.

Network answers behind F14 to F19 (NEODyS, IRSA, SsODNet, Horizons) were
recorded through [`snapshots.py`](snapshots.py) on that day, so the figures
can be reproduced after the services move on; the static tables (Toliou,
Granvik, the LCDB) are pinned by sha256 in their probes.

Every figure belongs to the catalog release named in its section.
[`figures.json`](figures.json) records the headline ones, and
[`run_all.py`](run_all.py) re-runs every probe and prints them beside the
recorded values.

| finding | what it measures | probe |
|---|---|---|
| [F11](#f11-a-guessed-taxonomy-is-right-about-half-the-time-and-a-guessed-c-less) | how often a guessed taxonomy is right | [`probe_taxonomy_prior.py`](probe_taxonomy_prior.py) |
| [F12](#f12-the-top-of-the-ranking-probably-never-passed-close-to-the-sun) | whether the water at the top of the ranking survived close passes to the Sun | [`probe_perihelion.py`](probe_perihelion.py) |
| [F13](#f13-the-bodies-that-lead-the-ranking-spin-too-fast-to-hold-loose-regolith) | whether those bodies spin too fast to hold loose regolith | [`probe_spin.py`](probe_spin.py) |
| [F14](#f14-orbit-uncertainty-barely-moves-delta-v-it-moves-where-the-body-is) | how much orbit uncertainty moves delta-v | [`probe_covariance.py`](probe_covariance.py) |
| [F15](#f15-most-of-the-top-of-the-ranking-cannot-be-found-again-before-2036) | whether the top-ranked bodies can be observed again before 2036, and where to look | [`probe_recovery.py`](probe_recovery.py) |
| [F16](#f16-the-guessed-neos-have-no-measured-albedo-anywhere-the-catalog-missed) | whether the albedo-guessed NEOs have a measured albedo somewhere | [`probe_albedo_coverage.py`](probe_albedo_coverage.py) |
| [F17](#f17-a-debiased-population-model-agrees-the-guessed-cs-are-mostly-not-dark) | what a debiased NEO population model says about the guessed classes | [`probe_source_regions.py`](probe_source_regions.py) |
| [F18](#f18-the-lightcurve-database-adds-no-spins-and-trims-f13s-headline) | what the lightcurve database adds to F13 | [`probe_lightcurves.py`](probe_lightcurves.py) |
| [F19](#f19-f4s-oracle-is-not-the-problem-its-two-body-error-is-00035-kms) | how much of F4's gap is the oracle's own two-body error | [`probe_horizons.py`](probe_horizons.py) |

---

## F11. A guessed taxonomy is right about half the time, and a guessed C less

Measured on 2026-10-29 with
[`probe_taxonomy_prior.py`](probe_taxonomy_prior.py), against catalog release
`data-2026-10-21b`.

The builder (AsteroidCatalog's `enrich.py`) types a body with no class in two
ways. With a measured albedo, below 0.10 is C and anything else S. With no
albedo, it takes the albedo its orbit bin assumes and applies the same split,
which for a NEO means S inside 2.5 AU and C beyond. That second rule types
1,298,012 bodies, and 38,614 of the 42,775 NEOs.

The probe mirrors both rules from the builder's source, and before measuring
it checks the mirror against the release itself. It reproduces the class of
all 83,404 albedo-typed rows and all 1,298,008 albedo_assumed rows within
4.6 AU. Its numpy multinomial logit recovers known coefficients from 60,000
synthetic rows to 0.036. Each check failed on a planted defect (the split moved
to 0.11, and the labels rotated by one class).

Classes are `comp_group` collapsed to four: C-complex, S-complex, X-complex
and Other. The calibration set is the 169,508 bodies within 4.6 AU whose class
came from a source, 2,943 of them NEOs.

### Fitted models barely beat the rules

On a 20% hold-out (four classes, so a rule that only says C or S is wrong on
every X and Other):

```
orbit model, all     n= 33,766  model acc 51.9%  log loss 1.120 (base rates 1.191)  rule acc 45.4%
orbit model, NEOs    n=    577  model acc 50.6%  log loss 1.294 (base rates 1.215)  rule acc 47.8%
albedo model, all    n= 12,720  model acc 70.8%  log loss 0.808 (base rates 1.147)  rule acc 69.7%
albedo model, NEOs   n=    167  model acc 68.3%  log loss 0.988 (base rates 1.106)  rule acc 70.7%
```

A measured albedo carries most of the information, and the builder's split
already captures nearly all of it. Without an albedo, the orbit tells you
little. On NEOs the fitted orbit model scores worse than plain base rates on
log loss, and its calibration inverts (rows it gives 0.45 are right 0.55 of
the time, rows it gives 0.56 are right 0.48). The builder's own split, scored
as C against everything else, is right for 86.2% of the 64,148 source-typed
bodies with an albedo, where its docstring records 84.0% on an earlier release
with an unrecorded scoring rule.

So `MNLogit` is not the prior for NEOs. The probe uses the empirical one: among
NEOs with a source class, how the true class falls for each class the
albedo_assumed rule would give.

```
rule C, H [0.0, 18.0)   n=  156  C 0.308  S 0.455  X 0.090  O 0.147
rule C, H [18.0, 20.0)  n=  140  C 0.300  S 0.436  X 0.100  O 0.164
rule C, H [20.0, 22.0)  n=   56  C 0.286  S 0.375  X 0.196  O 0.143
rule C, pooled          n=  384  C 0.294  S 0.453  X 0.104  O 0.148
rule S, H [0.0, 18.0)   n=  523  C 0.157  S 0.545  X 0.059  O 0.239
rule S, H [18.0, 20.0)  n=  838  C 0.193  S 0.529  X 0.103  O 0.175
rule S, H [20.0, 22.0)  n=  683  C 0.193  S 0.517  X 0.092  O 0.198
rule S, H [22.0, 25.0)  n=  333  C 0.159  S 0.492  X 0.147  O 0.201
rule S, H [25.0, 40.0)  n=  182  C 0.137  S 0.544  X 0.187  O 0.132
rule S, pooled          n=2,559  C 0.177  S 0.525  X 0.103  O 0.195
```

Where the rule says S, the body is S-complex about half the time. Where it
says C, it is C-complex 29% of the time and S-complex 45%. The rule C cells
past H 22 hold 29 and 3 bodies, and the pooled row stands in for them.

### What it does to the value-bearing fractions

Each class's fractions are the mean over its source-typed bodies. Applied to
the 38,614 albedo_assumed NEOs:

```
metal  assigned C: point median 0.010  expected median 0.074
metal  assigned S: point median 0.060  expected median 0.093
carbon assigned C: point median 0.040  expected median 0.022
carbon assigned S: point median 0.010  expected median 0.019
ice    assigned C: point median 0.150  expected median 0.055
ice    assigned S: point median 0.000  expected median 0.038
```

The 3,460 NEOs typed C from an assumed albedo carry about three times the ice
the prior supports. The 35,154 typed S carry none, and the prior gives them
some. Metal goes up in both, because 10 to 15% of either group is X-complex,
whose mean metal fraction is 0.373.

### At the top of the ranking

The profitability catalog on disk is from pipeline 1.26.0, run on 2026-09-30,
and every one of its 604,422 rows loses money. Its order is by profit, so its
top is the cheapest losing missions. The first 16 rows are NEOs of H 31.5 to
32.5, about a metre across at the assumed albedo, typed S from an assumed
albedo. They fall in the H >= 25 bin, which holds 185 NEOs with a measured
class, two of them at H 30 or fainter. Their P(S) of 0.544 is an
extrapolation.

```
top    25:   21 albedo_assumed NEOs   H >= 25.0:   19   P(assigned) median 0.544
top   100:   61 albedo_assumed NEOs   H >= 25.0:   25   P(assigned) median 0.294
top 1,000:  364 albedo_assumed NEOs   H >= 25.0:   98   P(assigned) median 0.294
```

Past the first sixteen, the rows that carry water are mostly C-types from an
assumed albedo, and the prior gives each of those 0.29.

### Cautions

Bodies with a measured class are not a random sample. Spectroscopic surveys
reach bright, large and S-type bodies first, so the empirical prior is too
generous to S and too hard on C. A debiased prior needs a survey model, which
nothing here has. The direction of F2 survives that caution: a guessed class
is a coin flip, and the model reports it as a fact. This probe does not say how
much the ranking moves, because pricing composition as a distribution would
mean running the value model per class, which is economicspace's job.

No module, config or CSV was touched, and no number a run produces moved.

---

## F12. The top of the ranking probably never passed close to the Sun

Measured on 2026-10-29 with [`probe_perihelion.py`](probe_perihelion.py),
against catalog release `data-2026-10-21b` and the Toliou, Granvik and
Tsirvoulis (2021) grid, downloaded that day (247,823,013 bytes, sha256
`8b895c66...a7ad86`, recorded in the probe).

The grid gives, for a NEO's (a, e, i, H), the probability that its population
ever passed inside each of 26 perihelion thresholds. Below about 0.2 AU,
heating strips hydrated minerals, so it bears on the water the model credits
to C-types. Two numbers per body: prospector's expected minimum perihelion,
and P(q_min < 0.2 AU) from the table's fourth column.

The reader and the expectation are adapted from prospector. The probe
reproduces prospector 0.4.3's own figures on five bodies to 1.1e-16 AU
(Bennu 0.569230, Itokawa 0.564210, Ryugu 0.576365, Eros 0.738385, Phaethon
0.047345). Its lookup snaps each axis to the nearest grid value instead of
scanning 458,169 rows per body, and it picks a cell exactly as near as
prospector's full scan on 1,500 sampled catalog NEOs. Both checks failed on
planted defects. No NEO's expected minimum came out more than 0.05 AU above
its present perihelion, which a minimum cannot exceed.

```
NEOs: 42,773   empty-cell fallbacks 451   H outside the grid 15,570

BY COMPOSITION GROUP
  C-complex           n=  4,585  q_min_exp median 0.675 AU   P(q_min<0.2) median 0.021  mean 0.072   > 0.5 on 2.9%
  S-complex           n= 37,302  q_min_exp median 0.616 AU   P(q_min<0.2) median 0.005  mean 0.046   > 0.5 on 1.6%
  X-complex           n=    303  q_min_exp median 0.597 AU   P(q_min<0.2) median 0.044  mean 0.111   > 0.5 on 4.3%
  credited with ice   n=  4,913  q_min_exp median 0.667 AU   P(q_min<0.2) median 0.022  mean 0.077   > 0.5 on 3.2%

PROFITABILITY RANKING
  top 25 carrying water        n=   5  P(q_min<0.2) median 0.000  mean 0.001
  top 100 carrying water       n=  49  P(q_min<0.2) median 0.011  mean 0.014
  top 1,000 carrying water     n= 364  P(q_min<0.2) median 0.012  mean 0.017   > 0.5 on 0.0%
```

None of the 364 water-carrying NEOs in the top thousand has better than even
odds of having passed inside 0.2 AU, and their mean is 1.7%. Across all NEOs
credited with ice, 3.2% do. The bodies that win are on orbits whose population
stays well clear of the Sun, which is what a cheap rendezvous would select
for, so the water credit at the top survives this check. F11 is the larger
doubt about that water.

Two cautions. The grid's H runs from 15.125 to 24.875, and 15,570 NEOs lie
outside it, mostly small ones that take its faint edge, including the
metre-scale bodies that lead the ranking. And the table describes the
population on a body's orbit, not the body; two bodies on one orbit can have
had different pasts.

No module, config or CSV was touched, and no number a run produces moved.

---

## F13. The bodies that lead the ranking spin too fast to hold loose regolith

Measured on 2026-10-29 with [`probe_spin.py`](probe_spin.py), against catalog
release `data-2026-10-21b`.

economicspace digs "run-of-mine regolith" from every body it ranks, and reads
`rotation_period_h` only to size the battery for the night. Nothing asks
whether the body can hold regolith. A sphere that spins faster than
sqrt(3 pi / (G rho)), about 2.3 h at 2,000 kg/m^3, cannot keep a loose grain
on its equator. Faster than Holsapple's (2007) cohesionless limit, a rubble
pile with a 32 degree friction angle yields unless cohesion holds it.

### prospector's cohesion model has a wrong exponent

The probe adapts prospector's `cohesion.py`, reduced to the sphere, where
Holsapple's three shape integrals have a closed form. The reduction exposed an
upstream bug. `_shape_integral` gives the x-axis integrand the exponents
(1.5, 0.5, 1.5) where Holsapple has (1.5, 0.5, 0.5). For a sphere that returns
A_x = 0.4 instead of 2/3, and the three integrals sum to 1.733 where they must
sum to 2. With A_x forced to 0.4, the probe's closed form reproduces
prospector 0.4.3's own figures on five cases to 2.3e-10, so that exponent is
the whole difference. With the correct 2/3, a strengthless rubble pile at
2,000 kg/m^3 survives down to 2.804 h, where upstream says 3.495 h; the closed
form agrees with a bisection on the cohesion to 1e-6 h. Both checks failed on
planted defects. On the NEOs with a measured period, upstream's limit would
call 1,115 of them fast where the corrected one calls 882.

### What the measured periods say

```
NEOs: 42,775   with a measured period, density and size: 2,437
  lift period         median 2.009 h   faster than it: 810 (33.2%)
  cohesionless period median 2.414 h   faster than it: 882 (36.2%)
  cohesion the fast ones need: median 23.33 Pa  p90 1286.1 Pa  max 309566 Pa;  over 100 Pa: 281

  by H (measured periods):
    H [0, 18)      n=  573  median period    5.786 h   faster than lift   2.6%
    H [18, 20)     n=  601  median period    7.000 h   faster than lift   4.5%
    H [20, 22)     n=  435  median period    6.270 h   faster than lift  17.0%
    H [22, 25)     n=  452  median period    0.372 h   faster than lift  74.6%
    H [25, 28)     n=  325  median period    0.089 h   faster than lift  94.2%
    H [28, 40)     n=   51  median period    0.024 h   faster than lift 100.0%
```

A third of the NEOs with a measured period cannot hold loose material on
their equator. Past H 22, most cannot; past H 28, every one of the 51 measured
spins in under two minutes at the median.

### At the top of the ranking

None of the top 25 rows has a measured period, and only 21 of the 407 NEOs in
the top thousand do (two of them faster than the lift period). The first 16
rows are bodies of H 31.5 to 32.5, about a metre across. Prospector's way to
fill a missing period is the share of bodies within 0.75 mag of the same H,
and at H 32 that window holds one or two measured bodies, too few to use. The
nearest measured population is the H >= 28 bin, where all 51 spin faster than
the lift period, at a median of 1.4 minutes. Bodies four magnitudes fainter
are smaller still, and small bodies spin faster. The metre-scale bodies that
lead the ranking are almost certainly bare, fast-spinning rock with no loose
layer to dig.

What that costs the model, this probe cannot say. It may mean a different
excavation method, or no mission at all, and either is a modelling decision
in economicspace. The cohesion they would need is small (a metre-scale body
spinning once a minute needs about 1 Pa), so they need not fly apart. What they
lack is the loose surface layer the excavation model assumes.

### Cautions

Lightcurve surveys find short periods more easily than long ones on small,
faint bodies, so the measured fast share at faint H is likely an overestimate.
At H >= 25 it is above 94%, so that would need to be very large to change the
reading. Every body is treated as a sphere of its catalog density estimate.
An elongated body sheds material from the tips of its long axis at a somewhat
slower spin than a sphere does.

No module, config or CSV was touched, and no number a run produces moved.

---

## F14. Orbit uncertainty barely moves delta-v; it moves where the body is

Measured on 2026-10-29 with [`probe_covariance.py`](probe_covariance.py),
against catalog release `data-2026-10-21b` and NEODyS as served that day.

SPACE-MAP.md found poorly determined orbits enriched 3.1x at the top of the
ranking and left open how to respond. F9 proposed pushing NEODyS's orbit
covariance through the delta-v estimator to get an error bar. The probe
fetches the Keplerian elements and 6x6 covariance of every NEO in the top 100
rows (69 bodies, all found), draws 4,000 times from the Gaussian on (a, e, i),
and prices each draw with the shipped estimator.

Before measuring, the shipped copy reproduces the five docstring figures to
0.0081 km/s. For every body fetched, the square roots of the parsed
covariance's diagonal match the RMS line NEODyS prints, to 4.5e-04 at worst.
That check failed when the 21 upper-triangle values were planted in
column-major order.

```
DELTA-V UNCERTAINTY  (shipped dv_out, 5th-95th percentile, 4,000 draws each)
  90% interval width: median 0.001 km/s   p90 0.088   max 4.443
  wider than 0.1 km/s:   6 of 69
  wider than 0.5 km/s:   2 of 69
  wider than 1.0 km/s:   2 of 69
  draws rejected as unbound or unphysical: median 0.00%   max 79.3%
  spread of nominal dv across these bodies: IQR 1.306 km/s (6.108 to 7.414)

  by orbit condition code:
    U 0-2  n= 19  width90 median    0.000 km/s   max    0.000
    U 3-4  n=  4  width90 median    0.000 km/s   max    0.000
    U 5-6  n= 16  width90 median    0.000 km/s   max    0.006
    U 7-9  n= 30  width90 median    0.035 km/s   max    4.443
```

For 67 of the 69, the 90% interval on `dv_out` is under 0.5 km/s, and for
most it is a few metres per second, well inside the 1.3 km/s spread between
bodies. Even at U 7 to 9 the median width is 35 m/s. The two exceptions are
U 9 with one-day arcs: 2015 KW157, 5.30 to 9.74 km/s, and 2000 EA116, where
79% of the draws are not bound orbits at all.

The reason is in the covariance. 2022 FB21, rank 1 at U 8, has its semi-major
axis to 0.002 AU, its eccentricity to 0.0012 and its inclination to 0.019
degrees, but its mean anomaly only to 2.2 degrees at an epoch in 2026. A short
arc taken at a close approach fixes the orbit's size, shape and tilt well. What
it leaves loose is where along the orbit the body is, and the closed-form
estimator never looks at that.

So the answer to SPACE-MAP.md's open question is narrower than F9 hoped. An
error bar on delta-v is cheap, and for almost every body at the top it is
negligible. What a poor orbit actually threatens is finding the body again: a
metre-scale object whose position along its orbit is uncertain by degrees today
can be lost within a few apparitions. That is a recovery and targeting cost,
which nothing in the model represents, and an arc floor or a U cutoff is the
crude stand-in for it.

### Cautions

A Gaussian on Keplerian elements is a linearisation, and for the one-day arcs
it is a poor one; their intervals mean "badly known", not a calibrated 90%.
Only (a, e, i) are drawn, because only they enter the estimator. NEODyS
elements are its own fit at its own epoch, so the nominal figures differ a
little from the catalog's.

No module, config or CSV was touched, and no number a run produces moved.

---

## F15. Most of the top of the ranking cannot be found again before 2036

Measured on 2026-10-29 with [`probe_recovery.py`](probe_recovery.py), against
catalog release `data-2026-10-21b` and NEODyS as served that day (replayable
from that day's snapshot folder; see README).

F14 found that a poorly determined orbit moves where a body is, not its
delta-v. This asks what that costs. A mission needs its target found again
before launch, and a metre-scale body is visible only during a close approach.
For each NEO in the top 100 rows the probe takes NEODyS's elements and full
covariance, propagates the nominal orbit two-body day by day against a
Keplerian Earth from 2026-10-29 to 2036, and looks for days when the body is
at V 24 or brighter (a Rubin-class limit, optimistic for today's surveys) and
at least 60 degrees from the Sun. Where there is such a window, 200 draws from
the covariance give the patch of sky a telescope would have to search: the
95th percentile of their angular distance from the nominal position.

Checks before measuring. The Keplerian Earth sits at heliocentric longitude
100.380 degrees at J2000, against 100.376 from the Sun's true longitude (its
mean longitude 280.460 plus the equation of centre, -0.084). For all 69 bodies
the linear growth of the mean anomaly's uncertainty over ten years matches
4,000 draws to 5%. Both checks failed on planted defects: Earth's mean anomaly
moved by 0.1 degree, and the a-M covariance term dropped, which moves the
ten-year spread of 2000 SG344 by 20%. The MPC's condition code U is defined from
that same ten-year drift, and the probe's drift ranks with the catalog's U at
Spearman 0.979.

```
top 100 rows: 69 NEOs; searching 2026-10-29 to 2036-01-01
LINEAR vs MONTE CARLO  sigma_M after ten years: 69 of 69 agree to 5% (OK)
  10-year mean-anomaly drift against the catalog's U, Spearman: 0.979 (n = 69)

RECOVERY BEFORE 2036-01-01  (V <= 24.0, elongation >= 60 deg)
  recoverable    21   median H  19.8   median p_observable 1.00
  search          1   median H  21.8   median p_observable 1.00
  lost            5   median H  22.5   median p_observable 0.96
  no window      42   median H  28.1   median p_observable 0.00

  by condition code:
    U 0-4  n= 23  recoverable 19  search 0  lost 0  no window 4
    U 5-6  n= 16  recoverable 2  search 1  lost 0  no window 13
    U 7-9  n= 30  recoverable 0  search 0  lost 5  no window 25
```

Recoverable means a window exists and the search radius is under a degree;
lost, that it is over ten. 42 of the 69 never come within V 24 of Earth's
telescopes before 2036 at all, the 16 metre-scale bodies that lead the ranking
among them (their brightest is V 29 to 33). Every one of the 30 bodies with
U 7 to 9 is either never observable or lost. The 21 recoverable bodies are
nearly all U 0 to 4, at a median H of 19.8.

So the bodies that lead the ranking are not just uncertain: within a decade no
telescope could find them to target. The model carries no recovery or
targeting cost and no notion of a body being lost, which is the gap SPACE-MAP.md
found and F14 sharpened. An arc floor or a U cutoff is the crude form of this
probe.

Two-body propagation ignores the planets, and every one of these bodies makes
close Earth approaches, which stretch an orbit's uncertainty further. The
search radii are therefore lower bounds, so "lost" is safe and "recoverable"
may be optimistic. NEODyS's own covariance is a linearisation, poor for the
one- and two-day arcs. The window search follows the nominal orbit, and
`p_observable` shows how often the draws agree with it.

No module, config or CSV was touched, and no number a run produces moved.

---

## F16. The guessed NEOs have no measured albedo anywhere the catalog missed

Measured on 2026-10-29 with [`probe_albedo_coverage.py`](probe_albedo_coverage.py),
against catalog release `data-2026-10-21b`, NEOWISE V2.0 and SsODNet as served
that day.

F11 found a class guessed from a measured albedo right about 70% of the time
and one guessed from an assumed albedo about half. The builder assumes an
albedo for 38,616 of the 42,775 NEOs. If any of them has a measured albedo
somewhere, using it is the cheapest improvement to F2. The probe looked in two
independent places.

NEOWISE V2.0, re-read directly from IRSA (all 183,412 rows, every orbit type,
joined on number and provisional designation against both the catalog's
designation and its `provisional_designation`, because V2.0 predates many
numberings). Before trusting an absence, the join had to find again at least
95% of the 1,180 NEOs the catalog itself matched to NEOWISE. It found 1,127
(95.5%). The first versions of the join found 50% and 89% and were refused:
they keyed on the wrong flag, then missed bodies numbered after V2.0.

SsODNet's datacloud, which holds every published diameter and albedo SsODNet
knows for a body, for all 364 guessed NEOs in the top 1,000 rows and 300 more
at random.

```
NEOWISE V2.0  257,011 designations listed, 255,132 with an albedo, all orbit types
  join check: NEOs the catalog matched to NEOWISE, found again: 1,127 / 1,180 (95.5%) (OK)
  albedo_assumed NEOs with a NEOWISE V2.0 albedo: 0

SsODNet DATACLOUD  364 albedo_assumed NEOs in the top 1,000 rows, plus 300 at random
  answered 653; failed 11; answers that also carried some other body, dropped: 27
  random     n= 297  with any albedo    0 (0.0%)
  top 1,000  n= 356  with any albedo    0 (0.0%)
```

None. F18 adds the lightcurve database, which holds no unflagged albedo for any
of them either. The NEOs typed from an assumed albedo have no measured albedo
in any of the three, so no reclassification is available from existing data.

The first datacloud reading found ten albedos, and all ten were wrong.
SsODNet resolves names loosely and can answer with several bodies: asked for
2019 XG1, it also returned (3104) Durer, which carries the alias 1955 XG1 and
eight albedos from 2011 to 2022. The probe now keeps only the body whose name,
number or alias is the designation asked for, and counts the others it drops
(27 answers). Anything else that reads the datacloud by name should do the
same.

No module, config or CSV was touched, and no number a run produces moved.

---

## F17. A debiased population model agrees: the guessed Cs are mostly not dark

Measured on 2026-10-29 with [`probe_source_regions.py`](probe_source_regions.py),
against catalog release `data-2026-10-21b` and Granvik et al.'s low-resolution
NEO model (downloaded that day, 33,342,324 bytes, sha256 recorded in the probe).

F11's prior is biased toward bright S-types, because bodies with a measured
class are the ones spectroscopy reached first. Granvik et al. (2018) model the
whole NEO population corrected for what surveys can see. For each cell of
(a, e, i, H), their model file gives the probability that a NEO there came from
each of seven source regions, and the debiased probability that its albedo is
at most 0.1, which is exactly the builder's line between C and S.

The model's grid is laid out like Toliou et al.'s, so the probe reuses F12's
lookup and its check against prospector's full scan (exact on 1,500 catalog
NEOs). The file has 773,280 cells, as its header says, and every cell's source
probabilities sum to 1 within 3e-5. The check failed on a planted defect: a
probability nudged by 0.01, or a row missing.

```
CALIBRATION  NEOs with a measured albedo and H inside the grid: 2,047
  dark (albedo <= 0.1): observed 0.384, model mean P 0.180
    model P 0.0-0.1  n=   80  mean P 0.076  observed dark 0.237
    model P 0.1-0.2  n=1,416  mean P 0.133  observed dark 0.319
    model P 0.2-0.3  n=  319  mean P 0.244  observed dark 0.495
    model P 0.3-0.5  n=  195  mean P 0.376  observed dark 0.682
    model P 0.5-1.0  n=   37  mean P 0.621  observed dark 0.676

THE GUESS  NEOs typed from an assumed albedo 38,614 (H inside the grid 23,303)
  assigned C-complex  n= 3,460  model P(dark) median 0.251  mean 0.325;  model disagrees on 2,959
  assigned S-complex  n=35,154  model P(dark) median 0.112  mean 0.117;  model disagrees on 0
  the builder calls 3,460 of them dark (C); the model expects 5,238 (13.6%)
```

The model ranks bodies well: the darker it predicts, the darker they are. Its
level runs low against the NEOs with a measured albedo, 0.18 against 0.38
observed. That is the direction the sample's own bias predicts, since most of
those albedos come from NEOWISE, and an infrared survey finds dark bodies more
easily than an optical one at the same H. The model is the debiased one; the
measured sample is not.

On the guessed NEOs it says what F11 said, with an instrument that has no
spectroscopic bias. Of the 3,460 NEOs the builder types C from an assumed
albedo, the model gives 2,959 better-than-even odds of not being dark at all
(median P(dark) 0.25, against F11's 0.29). And across all 38,614 it expects
5,238 dark bodies where the builder assigns 3,460. The builder's no-albedo rule
puts dark bodies in the wrong places: too many beyond 2.5 AU, too few inside.

At the top of the ranking, the two Jupiter-family-comet bodies near rows 21 and
23 (2015 KW157, 2026 EX2) get P(dark) of 0.99 and 1.00, as dark as anything in
the file. The metre-scale bodies that lead it lie far outside the grid's H
(15 to 25), so their values are the faint edge's, an extrapolation.

No module, config or CSV was touched, and no number a run produces moved.

---

## F18. The lightcurve database adds no spins, and trims F13's headline

Measured on 2026-10-29 with [`probe_lightcurves.py`](probe_lightcurves.py),
against catalog release `data-2026-10-21b` and the Asteroid Lightcurve Database
(Warner, Harris and Pravec), release of 2023-10-01, downloaded that day
(39,966,552 bytes, sha256 recorded in the probe).

The hope was more measured periods for F13, especially at faint H. The probe
reads the LCDB summary table from its own header line, so its fixed-width
layout is not copied from anyone, and before measuring it reproduces four rows
read off the table by eye (Ceres, Itokawa, Bennu, 2008 TC3). That check failed
when the header was shifted by one column.

```
COVERAGE  NEOs 42,775; with a catalog period 2,440; with a reliable LCDB period (U >= 2, no limit flag) 1,955
  LCDB periods for NEOs the catalog has none for: 0
  the catalog's 2,440 NEO periods in the LCDB: 2,176; by quality: U 3 927, U 2 1,031, U 1 110, a limit (> or <) 76, no U 66

AGREEMENT  where both have a period (n = 1,955): within 1% 99.9%, within 5% 99.9%; a factor 2 apart 0.1%

SPIN  (LCDB periods; faster than the lift period for the body's catalog density)
  H [0, 18)    n=  526  median period    5.8370 h   faster than lift   1.7%
  H [18, 20)   n=  532  median period    7.2315 h   faster than lift   3.2%
  H [20, 22)   n=  367  median period    6.7400 h   faster than lift  13.9%
  H [22, 25)   n=  303  median period    0.4670 h   faster than lift  68.6%
  H [25, 28)   n=  193  median period    0.0760 h   faster than lift  94.3%
  H [28, 30)   n=   31  median period    0.0239 h   faster than lift 100.0%
  H [30, 40)   n=    3  median period    0.0121 h   faster than lift 100.0%
  all     n=1,955  faster than lift 25.6%
  faintest: H 30.72; at H >= 30 the slowest period is 1.6 min
```

There are no new periods: the catalog's NEO periods already are the LCDB's,
agreeing to 1% on all but 0.1% of the overlap. What the LCDB adds is a quality
code. 264 of the catalog's 2,440 NEO periods are not in its 2023 release, and
of those that are, 110 are rated U 1 (tentative), 76 are only limits and 66
carry no rating. F13 used every one. On the 1,955 reliable periods alone, the
share faster than the lift period is 25.6%, not 33.2%: the fast-spinning third
was partly made of tentative periods. The faint bins do not move. Past H 25 it
is still 94% or more, and all 34 measured NEOs past H 28 still spin in under
two minutes, so F13's conclusion about the top of the ranking stands.
Public spin data run out at H 30.72, with three bodies past H 30.

The LCDB also carries albedos. Of its 2,265 that it does not itself flag as
assumed, none belongs to a NEO the catalog typed from an assumed albedo (F16).

No module, config or CSV was touched, and no number a run produces moved.

---

## F19. F4's oracle is not the problem: its two-body error is 0.0035 km/s

Measured on 2026-10-29 with [`probe_horizons.py`](probe_horizons.py), against
catalog release `data-2026-10-21b` and JPL Horizons as served that day.

F4 measured the shipped closed form about 1.3 km/s below a Lambert porkchop
oracle, and SPACE-MAP.md section 4 left open how much of that the oracle's own
two-body propagation contributes. SECOND-PASS.md proposed bounding it with
brahe and the `sb441` kernels, which reach only the 373 most massive asteroids,
and THIRD-PASS.md proposed ASSIST, which needs a Linux build. JPL Horizons
serves state vectors from its full n-body solution for any body, so neither is
needed. The probe runs F4's porkchop twice per body on the same grid: once with
the target on Horizons' trajectory, once two-body from Horizons' osculating
elements at the window's start. Earth is Horizons' Earth both times.

The sample is the first 40 bodies of probe_lambert's seed-1 sample on today's
catalog (mostly main belt, like F4's), and 30 NEOs, seed 1. Each window opens on
2027-01-01 and runs one synodic period of departures plus four years of
flight, with Horizons vectors every two days, interpolated by cubic Hermite
splines.

Checks before measuring. probe_lambert's docstring and solver checks pass.
Two-body propagation from Horizons' elements lands on Horizons' own start
position to 0.0000 km. Interpolation error is measured on nodes a sparser track
never saw, at four and eight days. The ratio between the two is a median 16.0,
so the fourth-power scaling holds, and carried down to two days the worst
error is 1,187 km, at most 1.5e-4 km/s in a 0.25-year Lambert solve. A first
version gated on the raw four-day error and failed at 16,128 km, on eccentric
orbits swinging through perihelion. That was the right failure for the wrong
spacing, and the check was rebuilt around the spacing the oracle uses. A linear
interpolant planted in its place fails it.

```
THE ORACLE'S TWO-BODY ERROR  (oracle with Horizons' trajectory minus two-body, km/s)
  F4 sample  n= 40  median +0.0008  mean +0.0014  median |delta| 0.0043  p90 |delta| 0.0165  max |delta| 0.036
                    two-body position drift over the window: median 0.0143 AU, max 0.103 AU
  NEOs       n= 30  median -0.0011  mean -0.0056  median |delta| 0.0024  p90 |delta| 0.0078  max |delta| 0.113
                    two-body position drift over the window: median 0.0108 AU, max 0.126 AU
  all        n= 70  median -0.0002  mean -0.0016  median |delta| 0.0035  p90 |delta| 0.0152  max |delta| 0.113
```

Two-body orbits drift a median 0.013 AU from the real ones over a window of
five years or more, and up to 0.126 AU, yet the porkchop's minimum barely
notices. Its median error is 0.0035 km/s, centred on zero, against F4's gap of
1.3 km/s, and the worst case, 2016 KB4 (e 0.67), is 0.113 km/s. F4's gap is the
closed form's, not an artefact of how the oracle was propagated. SPACE-MAP.md
section 4's open item is closed, without kernels and without a Linux build.

This bounds the oracle's two-body error only. F4's other caveat, that the
oracle is an unconstrained two-impulse optimum and so itself optimistic,
stands as written.

No module, config or CSV was touched, and no number a run produces moved.
