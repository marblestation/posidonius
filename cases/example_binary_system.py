import posidonius
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case3.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 1.0e6*365.25 # time [days] where simulation starts
    time_step = 0.08 # days
    #time_limit   = 4*time_step # days
    time_limit   = 365.25 * 1.0e8 # days
    historic_snapshot_period = 100.*365.25 # days
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    consider_tides = True
    consider_rotational_flattening = True
    consider_general_relativity = "Newhall1983"
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativity)

    star1_mass = 1.0 # Solar masses
    star1_radius_factor = 1.0
    star1_radius = star1_radius_factor * posidonius.constants.R_SUN

    star1_dissipation_factor = 4.992*3.845764e-2 # -66+64
    star1_dissipation_factor_scale = 1.0
    star1_radius_of_gyration_2 = 5.9e-2 # Sun
    star1_love_number = 0.03
    star1_fluid_love_number = star1_love_number
    star1_position = posidonius.Axes(0., 0., 0.)
    star1_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star1_rotation_period = 24.0 # hours
    star1_angular_frequency = posidonius.constants.TWO_PI/(star1_rotation_period/24.) # days^-1
    star1_spin = posidonius.Axes(0., 0., star1_angular_frequency)
    #star1_evolution_type = posidonius.GalletBolmont2017(star1_mass) # mass = 0.30 .. 1.40
    #star1_evolution_type = posidonius.BolmontMathis2016(star1_mass) # mass = 0.40 .. 1.40
    star1_evolution_type = posidonius.Baraffe2015(star1_mass) # mass = 0.01 .. 1.40
    #star1_evolution_type = posidonius.Leconte2011(star1_mass) # mass = 0.01 .. 0.08
    #star1_evolution_type = posidonius.Baraffe1998(star1_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star1_evolution_type = posidonius.LeconteChabrier2013() # Jupiter
    #star1_evolution_type = posidonius.NonEvolving()
    universe.add_particle(star1_mass, star1_radius, star1_dissipation_factor, star1_dissipation_factor_scale, star1_radius_of_gyration_2, star1_love_number, star1_fluid_love_number, star1_position, star1_velocity, star1_spin, star1_evolution_type)


    ############################################################################
    star2_mass = 0.08 # Solar masses
    star2_radius_factor = 0.845649342247916
    star2_radius = star2_radius_factor * posidonius.constants.R_SUN

    star2_dissipation_factor = 2.006*3.845764e4 # -60+64
    star2_dissipation_factor_scale = 1.0
    star2_radius_of_gyration_2 = 1.94e-1 # Brown dwarf
    star2_love_number = 0.307
    star2_fluid_love_number = star2_love_number

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 5.0;                             # semi-major axis (in AU)
    e = 0.1;                               # eccentricity
    i = 0.5 * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0.;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n) * posidonius.constants.DEG2RAD;                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    gm = posidonius.constants.G*(star2_mass+star1_mass);
    x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
    star2_position = posidonius.Axes(x, y, z)
    star2_velocity = posidonius.Axes(vx, vy, vz)

    #////// Initialization of star2ary spin
    star2_obliquity = 5.0 * posidonius.constants.DEG2RAD # 0.2 rad
    star2_rotation_period = 70. # hours
    star2_angular_frequency = posidonius.constants.TWO_PI/(star2_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #star2_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star1_mass*star2_mass, star2_position, star2_velocity)
    #star2_semi_major_axis = star2_keplerian_orbital_elements[0]
    #star2_eccentricity = star2_keplerian_orbital_elements[2]
    #star2_semi_major_axis = a
    #star2_eccentricity = e
    #star2_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(star2_semi_major_axis, star2_eccentricity, star1_mass, star2_mass) # days
    #star2_angular_frequency = posidonius.constants.TWO_PI/(star2_pseudo_synchronization_period) # days^-1
    star2_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star1_mass*star2_mass, star2_position, star2_velocity)
    star2_inclination = star2_keplerian_orbital_elements[3]
    star2_spin = posidonius.calculate_spin(star2_angular_frequency, star2_inclination, star2_obliquity)

    #star2_evolution_type = posidonius.GalletBolmont2017(star2_mass) # mass = 0.30 .. 1.40
    #star2_evolution_type = posidonius.BolmontMathis2016(star2_mass) # mass = 0.40 .. 1.40
    star2_evolution_type = posidonius.Baraffe2015(star2_mass) # mass = 0.01 .. 1.40
    #star2_evolution_type = posidonius.Leconte2011(star2_mass) # mass = 0.01 .. 0.08
    #star2_evolution_type = posidonius.Baraffe1998(star2_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star2_evolution_type = posidonius.LeconteChabrier2013() # Jupiter
    #star2_evolution_type = posidonius.NonEvolving()
    universe.add_particle(star2_mass, star2_radius, star2_dissipation_factor, star2_dissipation_factor_scale, star2_radius_of_gyration_2, star2_love_number, star2_fluid_love_number, star2_position, star2_velocity, star2_spin, star2_evolution_type)

    universe.write(filename, whfast_alternative_coordinates="DemocraticHeliocentric")


