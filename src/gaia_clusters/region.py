"""A box in galactic coordinates: the unit a region is downloaded, clustered and scored in."""
from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord


@dataclass(frozen=True)
class GalacticBox:
    l_min: float
    l_max: float
    b_min: float
    b_max: float

    @property
    def centre(self) -> tuple[float, float]:
        return (self.l_min + self.l_max) / 2, (self.b_min + self.b_max) / 2

    def contains(self, l, b) -> np.ndarray:
        l = np.asarray(l) % 360
        return (l >= self.l_min) & (l < self.l_max) & (np.asarray(b) >= self.b_min) & (np.asarray(b) < self.b_max)

    def icrs_corners(self) -> list[tuple[float, float]]:
        """Corners as (ra, dec). Joined by great circles they enclose the box: the edges of constant
        l are great circles, and those of constant b bow out by under 0.02° for a 10° side at |b| = 5°.
        Filter with `contains` after downloading."""
        c = SkyCoord(l=[self.l_min, self.l_max, self.l_max, self.l_min] * u.deg,
                     b=[self.b_min, self.b_min, self.b_max, self.b_max] * u.deg, frame="galactic").icrs
        return list(zip(c.ra.deg, c.dec.deg))


def galactic(ra, dec) -> tuple[np.ndarray, np.ndarray]:
    g = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg).galactic
    return g.l.deg, g.b.deg


AURIGA = GalacticBox(170, 180, -5, 5)
