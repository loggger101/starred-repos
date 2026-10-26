# Citations, sources and attribution

Where the data, borrowed code and algorithms in this repository came from,
and what each one obliges you to do. Extracted on 2026-10-20 from
[economicspace's `CITATIONS.md`](https://github.com/loggger101/economicspace/blob/main/CITATIONS.md), which keeps the
citations for the pipeline itself; the section numbers here follow it.

---

## 1. Sources read only by the probes

None of these is an economicspace pipeline input. No stage fetches them and
none of them reaches a CSV; they are the external oracles
[`SECOND-PASS.md`](SECOND-PASS.md) measures the model against.

- NASA JPL NHATS, `https://ssd-api.jpl.nasa.gov/nhats.api`. Optimised
  round-trip delta-v and mission duration for the near-Earth asteroids
  accessible to human spaceflight. Read by
  `probe_nhats.py`.
- NASA JPL Horizons, `https://ssd.jpl.nasa.gov/api/horizons.api`, and the
  NAIF kernel archive. Recommended in F7 for bounding the two-body error in
  F4's Lambert oracle; not yet used.
- NEODyS (University of Pisa and SpaceDyS),
  `https://newton.spacedys.com/neodys/`. Per-object orbital elements with the
  full 6x6 covariance. Recommended in F9; not yet used.
- Asterank, `http://www.asterank.com/api/asterank`, MIT. Used for its `dv`
  column only, which is the Shoemaker-Helin closed-form delta-v
  approximation. Its economic columns return degenerate values through the
  public endpoint and are not used.

  > Shoemaker, E. M. and Helin, E. F. (1978). *Earth-approaching asteroids as
  > targets for exploration.* NASA CP-2053, pp. 245-256.

- Bus-DeMeo asteroid taxonomy, PDS Small Bodies Node bundle
  `urn:nasa:pds:ast.bus-demeo.taxonomy`. Fetched and inspected in the second
  pass; not committed and not used by any stage.

  > DeMeo, F. E., Binzel, R. P., Slivan, S. M., and Bus, S. J. (2009). *An
  > extension of the Bus asteroid taxonomy into the near-infrared.* Icarus,
  > 202(1), 160-180.

---

## 2. Data committed to this repository

### SDSS-based Asteroid Taxonomy V1.1

`sdsstax_ast.tab`, 7.5 MB, is committed because it is static and small. It
gives 63,468 numbered asteroids an independent taxonomic class from u'g'r'i'z'
photometry.

Dataset, as the archive's own `CITATION_DESC` gives it:

> Hasselmann, P. H., Carvano, J. M., and Lazzaro, D., *SDSS-based Asteroid
> Taxonomy V1.1.* EAR-A-I0035-5-SDSSTAX-V1.1. NASA Planetary Data System,
> 2012.

The method paper, which is what a scientific claim should cite:

> Carvano, J. M., Hasselmann, P. H., Lazzaro, D., and Mothe-Diniz, T. (2010).
> *SDSS-based taxonomic classification and orbital distribution of main belt
> asteroids.* Astronomy & Astrophysics 510, A43.
> https://doi.org/10.1051/0004-6361/200913322

The underlying survey:

> Ivezic, Z. et al. (2001). *Solar System Objects Observed in the Sloan
> Digital Sky Survey Commissioning Data.* The Astronomical Journal 122, 2749.

- Archive: `https://sbnarchive.psi.edu/pds3/non_mission/EAR_A_I0035_5_SDSSTAX_V1_1/`
- sha256 of the committed file, unchanged since it was first committed on
  2026-09-03:
  `aba21eff3f973ca81d5b460fb158600690fe0f4223dd6c8b6e754216a9e2f9b5`. CI
  checks it on every push.

A published result that rests on this table should cite the papers above, and
not this repository. NASA PDS data is in the public domain, but the scientific
credit belongs to its authors.

The column layout comes from the archive's own `.lbl` files, by way of
`juliensimon/space-datasets` (MIT), which had already transcribed them. See
section 5.

---

## 3. Code adapted from other projects

### skyfield, MIT

`orbital.py`, the functions `kepler_E`, `true_anomaly`
and `elements_to_state`, are adapted from `skyfield/keplerlib.py`.

Upstream's notice, verbatim from its `LICENSE`:

```
Copyright © 2013–2018 Brandon Rhodes and available under the MIT license:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

What changed: `kepler_E` keeps skyfield's starter and quartic Newton step but
drops the scalar early return so it broadcasts over arrays, and clips `e` away
from 0 and 1 rather than raising. `elements_to_state` is `ele_to_vec`
restricted to the elliptical case and re-expressed in terms of the semi-major
axis rather than the semi-latus rectum, because that is the column this catalog
stores.

The obligation is that the copyright notice and the MIT permission notice
travel with any distribution. The notice above, in this file, meets it.

- `https://github.com/skyfielders/python-skyfield`

---

## 4. Algorithms implemented from published papers

An algorithm described in a paper is not a copyrightable expression, so
implementing one from its equations carries no licence obligation. It carries a
citation obligation, which is what this section is.

### Lambert's problem

`lambert_izzo` and its helpers in `orbital.py` are implemented from:

> Izzo, D. (2015). *Revisiting Lambert's problem.* Celestial Mechanics and
> Dynamical Astronomy 121(1), 1-15.

It was deliberately not copied from an implementation. The two obvious ones
to take it from were pykep (MPL-2.0) and nyx (AGPL-3.0), and nyx would have
made economicspace AGPL, which is why it came from the paper.

Validated against a published worked example:

> Vallado, D. A. *Fundamentals of Astrodynamics and Applications*, Example 7-5.

Agreement is 7.4e-07 km/s, which is the precision the reference is quoted to.
See [`FINDINGS.md`](FINDINGS.md) F4.

### Kepler's equation

The starter and quartic Newton iteration in `kepler_E` follow:

> arXiv:2108.03215

reached by way of skyfield's implementation, which is credited in section 3.

---

## 5. Repositories audited, and what was taken from each

This repository records a survey of 17 repositories. Most were read and
nothing was taken from them, and this table gives the licence position for
all of them.

They were read twice. The second pass, in [`SECOND-PASS.md`](SECOND-PASS.md),
moved rows here, because a repository can add exactly the thing it was
rejected for lacking. The licence column does not go stale.

| repository | licence | what was taken |
|---|---|---|
| skyfielders/python-skyfield | MIT | code, adapted; see section 3 |
| juliensimon/space-datasets | MIT | the PDS3 column layout for the SDSS table (section 2), and the JPL NHATS endpoint used by `probe_nhats.py` |
| astropy/astroquery | BSD-3 | nothing. It showed that IRSA's async TAP was the right approach; the implementation here is written against `requests`. Its NEODyS client is where the covariance service was found |
| julie-dujardin/space-map | AGPL-3.0 | nothing. It showed that SBDB publishes a field list, which is how `condition_code` was found missing |
| Z3Prover/z3 | MIT | nothing yet |
| Pyomo/pyomo | BSD-3 | nothing |
| pymc-devs/pymc | Apache-2.0 | nothing yet |
| mesa/mesa | Apache-2.0 | nothing |
| OpenSCvx/OpenSCvx | Apache-2.0 | nothing |
| esa/pygmo2 | MPL-2.0 | nothing |
| pola-rs/polars | MIT | nothing |
| typpo/spacekit | MIT | nothing. Its per-body Kepler solve is the technique behind economicspace's orbit diagram (`ui.py`), but that diagram is written against `orbital.py` and Altair, and no line of spacekit was read into it |
| duncaneddy/brahe | MIT | nothing yet. Its native SPICE reader and Horizons SPK client are the recommended way to bound F4's oracle |
| Small-Bodies-Node/pds4_tools | BSD-3 (Univ. of Maryland) | nothing |
| nyx-space/nyx | AGPL-3.0 | nothing |
| cuspaceflight/CamPyRoS | GPL-3.0 | nothing |
| CelestiaProject/Celestia | GPL-2.0 | nothing |

Nothing is derived from any of the four copyleft projects. They are listed
because they were read, and the listing is the record that they were read and
not copied.

Read each licence file rather than GitHub's label. `pds4_tools` shows no
licence on GitHub because its BSD-3 text sits in a file named `LICENSES`,
which the detector does not recognise. "No licence" means all rights reserved,
so the label errs in the dangerous direction.
