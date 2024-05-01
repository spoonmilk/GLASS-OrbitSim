# GLASS-OrbitSim
Final Orbital Simulations for the GLASS Mission, ENGN1760 Spring 2024

# Mission Spec
A simulation for the physics of an intrasolar mission propelled by a ground based laser

# Acknowledgments
[poliastro](https://docs.poliastro.space/en/stable/index.html), the open-source MIT astrodynamics toolkit, was our main library used for simulation. (From the website) "an open source (MIT) pure Python library for interactive Astrodynamics and Orbital Mechanics, with a focus on ease of use, speed, and quick visualization. It provides a simple and intuitive API, and handles physical quantities with units." Thanks for the contributors for developing such a useful tool.

# The GLASS Hard Problem
As proposed in GLASS' mission-spec, the payload-carrying sail is propulsed through the use of several ground-based laser stations utilizing powerful 500-1000 Kw lasers. The astrodynamics problem of this approach comes in the method of propulsion from a fixed point on the ground.

As we make a burn (or as the laser is shot at the probe) the probe's apogee increases, but its perigee decreases. The issue here is that, without enoguh stations to boost the probe, the perigee's altitude decrease will literally slam the satellite back to Earth. Our central research question is optimizing the number of ground stations such that the orbit can be continually boosted while avoiding the failure scenario described. 

# Compromises
<!-- TODO: Fill in info about prograde compromise burn -->


