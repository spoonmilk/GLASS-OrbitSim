from astropy import units as u
from astropy.time import Time

import numpy as np

from poliastro import iod
from poliastro.bodies import Earth, Moon
from poliastro.frames import Planes
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter
from poliastro.util import time_range
from poliastro.spacecraft import Spacecraft
# More info: https://plotly.com/python/renderers/
import plotly.io as pio
pio.renderers.default = "plotly_mimetype+notebook_connected"
from astropy.coordinates import solar_system_ephemeris

epoch = Time("2024-05-01 13:05:50", scale = "tdb")

plotter = OrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
plotter.plot_body_orbit(Earth, epoch, label = "Earth")
plotter.plot_body_orbit(Moon, epoch, label = "Moon")

