import numpy as np
import math

from astropy import units as u
from astropy.time import Time

from poliastro.bodies import Earth
from poliastro.earth import EarthSatellite
from poliastro.maneuver import Maneuver
from poliastro.earth.enums import EarthGravity
from poliastro.spacecraft import Spacecraft
from poliastro.twobody.orbit import Orbit
from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter2D

# ------- HELPER METHODS --------

# Units: m, m/s, kg
# Evil math
def compute_time_int(al, vi, mas) -> float:
    return 8.50938 * (10**-9) * (
        14609987571408820 * (al**2) * mas
        + math.sqrt(
            213451736836720196366365882318848 * (al**4) * (mas**2)
            + 7452623680028430827520 * (al**3) * (mas**3) * (vi**3)
        )
    ) ** (1 / 3) - (0.166214 * al * mas * vi) / (
        14609987571408820 * (al**2) * mas
        + math.sqrt(
            213451736836720196366365882318848 * (al**4) * (mas**2)
            + 7452623680028430827520 * (al**3) * (mas**3) * (vi**3)
        )
    ) ** (
        1 / 3
    )

# Units: s, kg, m/s
def compute_vf(dt, mas, vi) -> float:
    return math.sqrt(((4 * (10**8)) * dt) / mas + vi ^ 2)


def compute_dv(vi, al, mas) -> float:
    dt = compute_time_int(al, vi, mas)
    vf = compute_vf(dt, mas, vi)
    return vf - vi


# ------- DEFINE PARAMETERS --------

# Attracting body
attr = Earth
# Radius
sat_radius = [-6045, -3490, 2500] << u.km
# Velocity
sat_velocity = [-3.457, 6.618, 2.533] << u.km / u.s
# Coefficient of drag unapplicable, made constant
sat_C_D = 1
# Mass
sat_mass = 100 * u.kg
# Area
sat_area = (3.4 * (u.m**2)).to(u.km**2)
# Epoch
sat_t = Time.now
# Plane
sat_plane = Planes.EARTH_EQUATOR

# ------- DEFINE AN ORBIT & FRAME --------

sat_orbit = Orbit.circular(Earth, alt=400 << u.km)

# ------- DEFINE A CRAFT --------

glass_spacecraft = Spacecraft(sat_area, sat_C_D, sat_mass)
# Place craft in orbit
earth_satellite = EarthSatellite(sat_orbit, glass_spacecraft)
plotter = OrbitPlotter2D()

# ------- MANEUVERS --------

# Impulse template

# (time * u.unit, [x comp, y comp, z comp] * u.unit) -> Impulses applied as dv at instant t(i)

# Maneuver template

# man_x = Maneuver(impulse, impulse, impulse, impulse)

# man = Maneuver(
#     (0 * u.s, [0, 0.1, 0] * u.km / u.s),
#     (20 * u.s, [0, 0.1, 0] * u.km / u.s),
#     (40 * u.s, [0, 0.1, 0] * u.km / u.s),
#     (60 * u.s, [0, 0.1, 0] * u.km / u.s),
#     (80 * u.s, [0, 0.1, 0] * u.km / u.s),
# )
# man = Maneuver(
#     (100 * u.s, [0.3, 0, 0] * u.km / u.s),
# )
# plotter.plot(sat_orbit, label="Init orbit")
# sat_orbit2 = sat_orbit.apply_maneuver(man)
# plotter.plot(sat_orbit2, label="After maneuver")

print(compute_time_int(400000, 7000, 500))
print(sat_orbit.classical.__str__)
