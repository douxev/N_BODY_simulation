# N_BODY_simulation
On its way to code an NBODY simulation

The project is to simulate an NBODY system between a star, a planet and some planetesimals.

The main simulation is NBODY_RK_planet_att.py
Then you can run Repopulate_planetesimals_orbits.py to visualize the orbits of each planetesimal better (adds several planetesimals on existing orbits).

Solar_system_sim.py is a big approx of how the objects in the stellar system (before the asteroid belt) would interact.
It uses an non precise centroid to calculate global attraction forces to each object. Still needs a quick fix.

2body_RK.py is what helped me create the NBody simulation.
NBody_RK.py is quite the same as NBODY_RK_planet_att.py, but the planet does not interact with any of the planetesimals.

cartesian_to_kepler_elements.py is just a quick python module that extracts Orbital elements from speed and relative position to a massive object. It is thus used in the Repopulate_planetesimals_orbits.py script.
