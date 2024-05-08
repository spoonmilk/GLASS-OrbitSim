from astropy.time import Time
from astropy import units as u

from matplotlib import pyplot as plt

from poliastro.bodies import Earth
from poliastro.ephem import Ephem
from poliastro.frames import Planes
from poliastro.plotting import StaticOrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range

EPOCH = Time("2018-02-18 12:00:00", scale="tdb")

roadster = Ephem.from_horizons(
    "SpaceX Roadster",
    epochs=time_range(EPOCH, end=EPOCH + 360 * u.day),
    attractor=Earth,
    plane=Planes.EARTH_ECLIPTIC,
)
roadster

frame = StaticOrbitPlotter(plane = Planes.EARTH_ECLIPTIC)

frame.plot_body_orbit(Earth, EPOCH)
frame.plot_ephem(roadster, EPOCH, label="SpaceX Roadster", color="black")