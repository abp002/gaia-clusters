# Protocol: cluster ages from isochrones

Written and committed before any test fit. It is fixed in two steps: this text, and then, after the
calibration below, the frozen settings appended at the end in a separate commit, still before the test
run. Anything done after the test run that departs from it is reported as post hoc, next to the result
it would replace, never instead of it.

## Question

Do ages fitted directly to our groups' colour-magnitude diagrams agree with Hunt & Reffert's (2023)?
Their ages come from a neural network applied to each cluster's diagram, trained on synthetic clusters.
Ours come from fitting theoretical isochrones to the same kind of data, star by star. Agreement says
the ages the viewer shows can be reproduced by a second, plainer method; disagreement says where they
cannot.

This is not a test of HDBSCAN, and our ages are not better than theirs by construction: both rest on
stellar models.

## Clusters

- Every cluster recovered in Carina and in Auriga (the score's view: each from the run that recovers it
  best), with at least **30 stars with a BP − RP colour** in our group. That is 115 in Carina and 58 in
  Auriga, counted from the exported diagrams before any fit.
- The stars are our group's, not Hunt & Reffert's members: field stars our group takes in are part of
  what the fit has to cope with.
- **Calibration clusters**, used to develop and tune the fit and left out of the test: NGC 3532 and
  Trumpler 19 (Carina) and NGC 2099, M37 (Auriga). They are reported, apart.

## Isochrones

PARSEC v1.2S with COLIBRI, from the CMD 3.9 service (stev.oapd.inaf.it), Gaia EDR3 passbands
(YBC, Vega), Kroupa IMF, at **solar metallicity, Z = 0.0152, fixed**: log10(age / yr) from 6.60 to
10.05 in steps of 0.05. Open clusters near the Sun are close to solar; metallicity is not fitted, so its
spread goes into the scatter.

Extinction is applied with the service's coefficients for a G2V star and the Cardelli et al. (1989) +
O'Donnell (1994) curve, R_V = 3.1: A_G = 0.83627 A_V, A_BP = 1.08337 A_V, A_RP = 0.63439 A_V. They do
not change with the star's colour, a simplification for hot and very red stars.

## Fit

For each cluster and each age on the grid, the extinction A_V (0 to 6 mag) and the distance modulus are
fitted by maximum likelihood; the age is the grid point with the highest likelihood.
- Each star's likelihood is a mixture: with weight 1 − f, a Gaussian in its distance to the isochrone in
  the colour-magnitude plane (widths σ_colour and σ_G); with weight f, a uniform density over the
  diagram's box, for field stars and other outliers.
- The distance modulus has a Gaussian prior from the group's median parallax, with 0.017 mas added
  for Gaia's zero point and an uncertainty that adds the median's error and 0.01 mas of systematics.
- The interval is the 16th to 84th percentile of the best ages over 50 bootstrap resamples of the
  group's stars.
- Settings fixed at calibration: σ_colour, σ_G, f, the stellar phases of the isochrone kept, and any
  magnitude limit. Nothing else is tuned.

**No age** is given, and the cluster is counted as such, when the best age is at either end of the grid
or the interval is wider than 0.5 dex.

## Comparison

Against Hunt & Reffert's median age, logAge50, for every test cluster with an age from us:
Δ = our age − logAge50, in dex.

**Primary result:** the share of test clusters with |Δ| ≤ 0.3 dex.

- **The ages agree if at least two thirds of them are within 0.3 dex**, and the Spearman rank
  correlation is at least 0.8.
- Below that, they do not, and the README says so next to the age colouring in the viewer.

0.3 dex is a factor of two in age: a loose bar, chosen before seeing any fit, because the viewer only
shows ages in bands a factor of ten wide. A second method that misses by more than that for a third of
the clusters would not support even those bands.

## Secondary

Reported in full whatever they show:
- the share of the clusters given an age, and why the rest were not;
- the median Δ (a bias), by age band of the viewer (under 10 Myr, 10–100 Myr, 0.1–1 Gyr, over 1 Gyr);
- Δ against the number of stars, distance and Hunt & Reffert's extinction;
- our extinction against their AV50, and our distance modulus against their MOD50;
- each region on its own;
- the calibration clusters, apart.

## Expectations

Written down so the result can contradict them. They are not part of the test above.
- The youngest clusters agree worst: with no turnoff, the age rests on a few bright stars and the
  pre-main sequence, which our photometry barely reaches.
- Where extinction is high, age and extinction trade off, and the intervals widen.
- Most clusters with fewer than 50 stars are given an age, but a loose one.

## Frozen settings (after calibration, before the test run)

Committed with the code that fits (`src/gaia_clusters/ages.py`, `Settings`), after fitting only the
three calibration clusters and before any test cluster.

- The likelihood became a model of where stars should be, not only how far each is from the curve: the
  isochrone's points are weighted by the share of stars the Kroupa IMF puts at each (PARSEC's int_IMF)
  and normalised over the magnitudes kept. With plain distance to the curve, NGC 3532 came out at 8.10,
  on an isochrone whose turnoff sits a magnitude above its brightest star and whose giants miss the
  three it has: a younger isochrone fits the main sequence just as well and pays nothing for the stars it
  predicts and the data lack.
- Widths: σ_colour = 0.05 mag, σ_G = 0.15 mag. Uniform component: f = 0.2. Each star sums the 32
  nearest isochrone points, densified to 0.2 widths apart.
- Phases: PARSEC labels 0 to 7, from the pre-main sequence to the early AGB; thermal pulses left out.
- Stars and model cut to 5 ≤ G ≤ 18: Gaia saturates at the bright end, and at the faint end our
  parallax cut (parallax_over_error > 5) loses stars in a way that depends on distance, which the IMF
  weighting would read as a shortage of low-mass stars.
- No binary sequence: unresolved binaries go to the uniform component.

What calibration showed, for the record: NGC 3532 8.65, Trumpler 19 9.60, M37 8.85, against Hunt &
Reffert's 8.38, 9.34 and 8.69, so 0.16 to 0.27 dex older in all three. In NGC 3532 and M37 Hunt &
Reffert's isochrone puts the turnoff above the brightest stars in our groups; in M37 ours runs below
its giants. The settings were not tuned to close that gap.
