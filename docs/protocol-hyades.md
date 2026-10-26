# Protocol: the Hyades' tidal tails

Written and committed before the sample was downloaded. Anything done after the first run that departs
from it is reported as post hoc, next to the result it would replace, never instead of it.

## Question

Can the Hyades' tidal tails be traced with the tools of this project, a velocity cut and HDBSCAN, and how
much of the published tails comes back? This is a reproduction: Röser et al. (2019), Meingast & Alves
(2019) and Jerabkova et al. (2021) found the tails, the last out to about 800 pc end to end.

Why the regions' method cannot work here: at 47 pc the Hyades and their tails spread over much of the
sky, so stars that share one space velocity have proper motions that change with their position. In
(position, parallax, proper motion) they are not a clump.

## Sample

Gaia DR3 from the GAVO Data Center (`gaia.dr3lite`), all sky: parallax > 2 mas (within 500 pc),
parallax_over_error > 10, RUWE < 1.4. All 862 tail members of Jerabkova et al. (2021, Gaia EDR3, the
same source IDs as DR3) pass these cuts; 99 % of them are within 470 pc. Their parallaxes were looked up
to choose the volume; nothing else was.

To keep the download small, the query keeps only stars whose tangential velocity is within 10 km/s of
the Hyades' (below). The selection itself is stricter.

## The Hyades' velocity

From Hunt & Reffert (2023), Melotte 25: (α, δ) = (66.7086°, 16.0808°), parallax 21.232 mas, proper
motion (104.136, −28.732) mas/yr, radial velocity 39.26 km/s. In ICRS Cartesian coordinates,
**V = (−5.735, 45.478, 4.712) km/s**. Its centre is at their position and distance, 47.19 pc; their
tidal radius is 8.62 pc.

## Selection

1. **Convergent point.** For each star, the Hyades' velocity projected on the sky at the star's position
   gives the tangential velocity (along α and δ, km/s) the star would have if it moved with the cluster.
   The star's own is 4.74047 μ / ϖ. A star is kept when the two differ by Δ ≤ √(3² + (2σ)²) km/s,
   where σ is the error that its parallax puts on its tangential velocity, |v_tan| σ_ϖ / ϖ. With no
   radial velocities for most stars, this is all the velocity there is.
2. **HDBSCAN** (`hdbscan` 0.8.44) on the kept stars' heliocentric positions in parsecs (x, y, z from
   (α, δ) and 1/ϖ), unscaled: min_cluster_size = 20, min_samples = 10, **excess of mass**. The Hyades'
   structure is the group holding the most Gaia Catalogue of Nearby Stars Hyades members. Excess of
   mass, not leaf: the aim is the whole structure, core and tails as one group, which leaf would cut
   into pieces.
3. **Our tail stars** are the stars of that group more than two tidal radii (17.2 pc) from the centre.
   Röser et al. (2019) drew their core at 18 pc.

## References

- **Jerabkova et al. (2021)**, Table B.2 (`J/A+A/647/A137/tableb2`): 862 tail members in Gaia EDR3.
- **Röser et al. (2019)** (`J/A+A/621/L2`): 292 stars in the leading tail and 237 in the trailing one,
  Gaia DR2. DR2 source IDs are matched to DR3 by equal ID; the share that matches is reported.
- **Gaia Catalogue of Nearby Stars** (Gaia Collaboration, Smart et al. 2021, `J/A+A/649/A6/hyacomb`): 920
  Hyades members within 100 pc, used to name our group and as a check on the core.

A reference star counts only if it is in our sample (the cuts above, before the velocity cut), and as a
tail star only if it is more than 17.2 pc from the centre.

## Results

**Primary: the recovery of Jerabkova et al.'s tail**: the share of their tail stars in our sample that
are among our tail stars.

- **The tails are reproduced if that recovery is at least 50 %, and at least half of our tail stars
  within their volume (500 pc) are among their 862.**
- Below either, they are not, and the README says which.

Secondary, reported in full whatever they show:
- the same recovery for Röser et al.'s leading and trailing tails, each;
- the recovery of the GCNS core members (within 17.2 pc);
- how far the tails reach: the 95th percentile of our tail stars' distance from the centre, against
  Jerabkova et al.'s, for each tail: the stars ahead of the centre and behind it along the direction of
  Galactic rotation (l = 90°, b = 0°), which the leading tail points to;
- the number of stars the velocity cut keeps, and the group's size.

## Expectations

Written down so the result can contradict them. They are not part of the test above.
- The core comes back whole, and the near parts of the tails (within 200 pc of the Sun) mostly do.
- The far ends do not: there the tails are a few stars in a large volume, and excess of mass will not
  attach them to the cluster.
- Purity falls with distance, as field stars that happen to share the velocity outnumber tail stars.
