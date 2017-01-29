import posidonius
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    initial_time = 4.5e6*365.25 # time [days] where simulation starts
    time_step = 0.08 # days
    #time_limit   = 4*time_step # days
    time_limit   = 365.25 * 1.0e8 # days
    historic_snapshot_period = 100.*365.25 # days
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    consider_tides = True
    consider_rotational_flattening = True
    consider_general_relativy = True
    consider_all_body_interactions = False
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativy, consider_all_body_interactions)

    star_mass = 0.08 # Solar masses
    star_radius_factor = 0.845649342247916
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_dissipation_factor = 2.006*3.845764e4 # -60+64
    star_dissipation_factor_scale = 1.0
    star_radius_of_gyration_2 = 0.2376400515700609
    star_love_number = 0.307
    fluid_love_number = star_love_number
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 70.0 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)
    #star_evolution_type = posidonius.BrownDwarf(star_mass)
    #star_evolution_type = posidonius.MDwarf()
    #star_evolution_type = posidonius.Jupyter()
    #star_evolution_type = posidonius.SolarLikeEvolvingDissipation(star_mass)
    #star_evolution_type = posidonius.SolarLikeConstantDissipation()
    star_evolution_type = posidonius.NonEvolving()
    universe.add_particle(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, star_love_number, fluid_love_number, star_position, star_velocity, star_spin, star_evolution_type)

    ############################################################################
    inner_planet_mass_factor = 1.0
    inner_planet_mass = inner_planet_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)
    # inner_planetary radius in AU (rearth in AU) Rocky inner_planet
    inner_planet_radius_factor = 1.
    inner_planet_radius = inner_planet_radius_factor * posidonius.constants.R_EARTH
    # Terrestrial:
    k2pdelta = 2.465278e-3 # Terrestrial inner_planets (no gas)
    inner_planet_dissipation_factor = 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(inner_planet_radius, 5))
    inner_planet_dissipation_factor_scale = 1.0
    inner_planet_radius_of_gyration_2 = 0.3308
    inner_planet_love_number = 0.305
    inner_planet_fluid_love_number = inner_planet_love_number

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.018;                             # semi-major axis (in AU)
    e = 0.1;                               # eccentricity
    i = 0. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0.;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n) * posidonius.constants.DEG2RAD;                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    gm = posidonius.constants.G*(inner_planet_mass+star_mass);
    x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
    inner_planet_position = posidonius.Axes(x, y, z)
    inner_planet_velocity = posidonius.Axes(vx, vy, vz)

    #////// Initialization of inner_planetary spin
    inner_planet_obliquity = 11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    inner_planet_rotation_period = 24. # hours
    inner_planet_angular_frequency = posidonius.constants.TWO_PI/(inner_planet_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*inner_planet_mass, inner_planet_position, inner_planet_velocity)
    #semi_major_axis = inner_planet_keplerian_orbital_elements[0]
    #eccentricity = inner_planet_keplerian_orbital_elements[2]
    #semi_major_axis = a
    #eccentricity = e
    #pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(inner_planet_semi_major_axis, inner_planet_eccentricity, star_mass, inner_planet_mass)
    #angular_frequency = posidonius.constants.TWO_PI/(inner_planet_pseudo_synchronization_period/24.) # days^-1
    inner_planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*inner_planet_mass, inner_planet_position, inner_planet_velocity)
    inner_planet_inclination = inner_planet_keplerian_orbital_elements[3]
    inner_planet_spin = posidonius.calculate_spin(inner_planet_angular_frequency, inner_planet_inclination, inner_planet_obliquity, inner_planet_position, inner_planet_velocity)

    inner_planet_evolution_type = posidonius.NonEvolving()
    universe.add_particle(inner_planet_mass, inner_planet_radius, inner_planet_dissipation_factor, inner_planet_dissipation_factor_scale, inner_planet_radius_of_gyration_2, inner_planet_love_number, inner_planet_fluid_love_number, inner_planet_position, inner_planet_velocity, inner_planet_spin, inner_planet_evolution_type)

    ############################################################################
    outer_planet_mass_factor = 1.0
    outer_planet_mass = outer_planet_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)
    # outer_planetary radius in AU (rearth in AU) Rocky outer_planet
    outer_planet_radius_factor = 1.
    outer_planet_radius = outer_planet_radius_factor * posidonius.constants.R_EARTH
    # Terrestrial:
    k2pdelta = 2.465278e-3 # Terrestrial outer_planets (no gas)
    outer_planet_dissipation_factor = 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(outer_planet_radius, 5))
    outer_planet_dissipation_factor_scale = 1.0
    outer_planet_radius_of_gyration_2 = 0.3308
    outer_planet_love_number = 0.305
    outer_planet_fluid_love_number = outer_planet_love_number

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.025;                             # semi-major axis (in AU)
    e = 0.01;                               # eccentricity
    i = 1. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0.;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n) * posidonius.constants.DEG2RAD;                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    gm = posidonius.constants.G*(outer_planet_mass+star_mass);
    x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
    outer_planet_position = posidonius.Axes(x, y, z)
    outer_planet_velocity = posidonius.Axes(vx, vy, vz)

    #////// Initialization of outer_planetary spin
    outer_planet_obliquity = 23. * posidonius.constants.DEG2RAD # 0.2 rad
    outer_planet_rotation_period = 24. # hours
    outer_planet_angular_frequency = posidonius.constants.TWO_PI/(outer_planet_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*outer_planet_mass, outer_planet_position, outer_planet_velocity)
    #semi_major_axis = outer_planet_keplerian_orbital_elements[0]
    #eccentricity = outer_planet_keplerian_orbital_elements[2]
    #semi_major_axis = a
    #eccentricity = e
    #pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(outer_planet_semi_major_axis, outer_planet_eccentricity, star_mass, outer_planet_mass)
    #angular_frequency = posidonius.constants.TWO_PI/(outer_planet_pseudo_synchronization_period/24.) # days^-1
    outer_planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*outer_planet_mass, outer_planet_position, outer_planet_velocity)
    outer_planet_inclination = outer_planet_keplerian_orbital_elements[3]
    outer_planet_spin = posidonius.calculate_spin(outer_planet_angular_frequency, outer_planet_inclination, outer_planet_obliquity, outer_planet_position, outer_planet_velocity)

    outer_planet_evolution_type = posidonius.NonEvolving()
    universe.add_particle(outer_planet_mass, outer_planet_radius, outer_planet_dissipation_factor, outer_planet_dissipation_factor_scale, outer_planet_radius_of_gyration_2, outer_planet_love_number, outer_planet_fluid_love_number, outer_planet_position, outer_planet_velocity, outer_planet_spin, outer_planet_evolution_type)

    universe.write(filename)


