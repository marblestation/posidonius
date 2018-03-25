import posidonius
import numpy as np
import argparse

def calculate_spin(angular_frequency, inclination, obliquity, position, velocity):
    """
    Old version of calculate spin used in Bolmont's thesis 2013
    """
    if inclination == 0.:
        # No inclination, spin can already be calculated:
        x = angular_frequency * np.sin(obliquity) # zero if there is no obliquity
        y = 0.
        z = angular_frequency * np.cos(obliquity)
    else:
        # Calculation of orbital angular momentum (without mass and in AU^2/day)
        horb_x = position.y() * velocity.z() - position.z() * velocity.y()
        horb_y = position.z() * velocity.x() - position.x() * velocity.z()
        horb_z = position.x() * velocity.y() - position.y() * velocity.x()
        horbn = np.sqrt(np.power(horb_x, 2) + np.power(horb_y, 2) + np.power(horb_z, 2))
        # Spin taking into consideration the inclination:
        x = angular_frequency * (horb_x / (horbn * np.sin(inclination))) * np.sin(obliquity+inclination)
        y = angular_frequency * (horb_y / (horbn * np.sin(inclination))) * np.sin(obliquity+inclination)
        z = angular_frequency * np.cos(obliquity+inclination)
    return posidonius.Axes(x, y, z)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 1.0e6*365.25 # time [days] where simulation starts
    time_step = 0.08 # days
    #time_limit   = 4*time_step # days
    time_limit   = 365.25 * 1.0e8 # days
    historic_snapshot_period = 100.*365.25 # days
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    consider_tides = True
    consider_rotational_flattening = False
    #consider_general_relativity = False
    consider_general_relativity = "Kidder1995" # Assumes one central massive body
    #consider_general_relativity = "Anderson1975" # Assumes one central massive body
    #consider_general_relativity = "Newhall1983" # Considers all bodies
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativity)

    star_mass = 0.08 # Solar masses
    star_radius_factor = 0.845649342247916
    # [start correction] -------------------------------------------------------
    # To reproduce Bolmont's thesis 2013 (although in this case it is not strictly needed because the star evolves):
    #   Mercury-T was using R_SUN of 4.67920694e-3 AU which is not the IAU accepted value
    #   thus, to reproduce Mercury-T results the star radius factor should be slighly modified:
    star_radius_factor = star_radius_factor*(4.67920694e-3 / posidonius.constants.R_SUN)
    # [end correction] ---------------------------------------------------------
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_dissipation_factor = 2.006*3.845764e4 # -60+64
    star_dissipation_factor_scale = 10.
    star_radius_of_gyration_2 = 1.94e-1 # Brown dwarf
    star_love_number = 0.307
    fluid_love_number = star_love_number
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 70.0 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)
    #star_evolution_type = posidonius.GalletBolmont2017(star_mass) # mass = 0.30 .. 1.40
    #star_evolution_type = posidonius.BolmontMathis2016(star_mass) # mass = 0.40 .. 1.40
    #star_evolution_type = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
    star_evolution_type = posidonius.Leconte2011(star_mass) # mass = 0.01 .. 0.08
    #star_evolution_type = posidonius.Baraffe1998(star_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star_evolution_type = posidonius.LeconteChabrier2013() # Jupiter
    #star_evolution_type = posidonius.NonEvolving()

    #universe.add_particle(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, star_love_number, fluid_love_number, star_position, star_velocity, star_spin, star_evolution_type)
    universe.add_brown_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_evolution_type)

    ############################################################################
    planet1_mass_factor = 1.0
    # [start correction] -------------------------------------------------------
    # To reproduce Bolmont's thesis 2013:
    #   Mercury-T was using planet1_mass as 3.00e-6 M_SUN and that's not exactly 1 M_EARTH (as accepted by IAU)
    #   thus, to reproduce Mercury-T results the mass factor should be slighly modified:
    planet1_mass_factor = planet1_mass_factor * (3.00e-6 / posidonius.constants.M_EARTH) # 0.999000999000999
    # [end correction] ---------------------------------------------------------
    planet1_mass = planet1_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    # Earth-like => mass-radius relationship from Fortney 2007
    planet1_radius_factor = posidonius.tools.mass_radius_relation(planet1_mass_factor, planet_mass_type='factor', planet_percent_rock=0.70)
    # [start correction] -------------------------------------------------------
    # To reproduce Bolmont's thesis 2013:
    #   Mercury-T defined a different M2EARTH from the IAU accepted value
    #   and M2EARTH was used to compute planet1_radius_factor, thus to reproduce
    #   Mercury-T results the planet1_radius_factor has to be corrected:
    planet1_radius_factor = planet1_radius_factor * 0.999756053794 # 1.0097617465214679
    # [end correction] ---------------------------------------------------------
    planet1_radius = planet1_radius_factor * posidonius.constants.R_EARTH

    # Terrestrial:
    k2pdelta = 2.465278e-3 # Terrestrial planet1s (no gas)
    planet1_dissipation_factor = 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet1_radius, 5))
    planet1_dissipation_factor_scale = 1.0
    planet1_radius_of_gyration_2 = 0.3308
    planet1_love_number = 0.305
    planet1_fluid_love_number = planet1_love_number

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.018;                             # semi-major axis (in AU)
    e = 0.00001;                              # eccentricity
    i = 5. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    gm = posidonius.constants.G*(planet1_mass+star_mass);
    x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
    planet1_position = posidonius.Axes(x, y, z)
    planet1_velocity = posidonius.Axes(vx, vy, vz)

    #////// Initialization of planet1ary spin
    planet1_obliquity = 0.01#11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    planet1_rotation_period = 24. # hours
    planet1_angular_frequency = posidonius.constants.TWO_PI/(planet1_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #planet1_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*planet1_mass, planet1_position, planet1_velocity)
    #planet1_semi_major_axis = planet1_keplerian_orbital_elements[0]
    #planet1_eccentricity = planet1_keplerian_orbital_elements[2]
    #planet1_semi_major_axis = a
    #planet1_eccentricity = e
    #planet1_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(planet1_semi_major_axis, planet1_eccentricity, star_mass, planet1_mass) # days
    #planet1_angular_frequency = posidonius.constants.TWO_PI/(planet1_pseudo_synchronization_period) # days^-1
    planet1_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*planet1_mass, planet1_position, planet1_velocity)
    planet1_inclination = planet1_keplerian_orbital_elements[3]
    planet1_spin = calculate_spin(planet1_angular_frequency, planet1_inclination, planet1_obliquity, planet1_position, planet1_velocity)

    planet1_evolution_type = posidonius.NonEvolving()
    #universe.add_particle(planet1_mass, planet1_radius, planet1_dissipation_factor, planet1_dissipation_factor_scale, planet1_radius_of_gyration_2, planet1_love_number, planet1_fluid_love_number, planet1_position, planet1_velocity, planet1_spin, planet1_evolution_type)
    universe.add_earth_like(planet1_mass, planet1_dissipation_factor_scale, planet1_position, planet1_velocity, planet1_spin, planet1_evolution_type)

    ############################################################################
    planet2_mass_factor = 10.0
    # [start correction] -------------------------------------------------------
    # To reproduce Bolmont's thesis 2013:
    #   Mercury-T was using planet2_mass as 3.00e-6 M_SUN and that's not exactly 1 M_EARTH (as accepted by IAU)
    #   thus, to reproduce Mercury-T results the mass factor should be slighly modified:
    planet2_mass_factor = planet2_mass_factor * (3.00e-6 / posidonius.constants.M_EARTH) # 0.999000999000999
    # [end correction] ---------------------------------------------------------
    planet2_mass = planet2_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    # Earth-like => mass-radius relationship from Fortney 2007
    planet2_radius_factor = posidonius.tools.mass_radius_relation(planet2_mass_factor, planet_mass_type='factor', planet_percent_rock=0.70)
    # [start correction] -------------------------------------------------------
    # To reproduce Bolmont's thesis 2013:
    #   Mercury-T defined a different M2EARTH from the IAU accepted value
    #   and M2EARTH was used to compute planet2_radius_factor, thus to reproduce
    #   Mercury-T results the planet2_radius_factor has to be corrected:
    planet2_radius_factor = planet2_radius_factor * 1.00046285582 # 1.8070338480688148
    # [end correction] ---------------------------------------------------------
    planet2_radius = planet2_radius_factor * posidonius.constants.R_EARTH

    # Terrestrial:
    k2pdelta = 2.465278e-3 # Terrestrial planet2s (no gas)
    planet2_dissipation_factor = 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet2_radius, 5))
    planet2_dissipation_factor_scale = 10.0
    planet2_radius_of_gyration_2 = 0.3308
    planet2_love_number = 0.305
    planet2_fluid_love_number = planet2_love_number

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.1;                             # semi-major axis (in AU)
    e = 0.5;                               # eccentricity
    i = 1. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    gm = posidonius.constants.G*(planet2_mass+star_mass);
    x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
    planet2_position = posidonius.Axes(x, y, z)
    planet2_velocity = posidonius.Axes(vx, vy, vz)

    #////// Initialization of planet2ary spin
    planet2_obliquity = 0.2 #23. * posidonius.constants.DEG2RAD # 0.2 rad
    planet2_rotation_period = 24. # hours
    planet2_angular_frequency = posidonius.constants.TWO_PI/(planet2_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #planet2_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*planet2_mass, planet2_position, planet2_velocity)
    #planet2_semi_major_axis = planet2_keplerian_orbital_elements[0]
    #planet2_eccentricity = planet2_keplerian_orbital_elements[2]
    #planet2_semi_major_axis = a
    #planet2_eccentricity = e
    #planet2_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(planet2_semi_major_axis, planet2_eccentricity, star_mass, planet2_mass) # days
    #planet2_angular_frequency = posidonius.constants.TWO_PI/(planet2_pseudo_synchronization_period) # days^-1
    planet2_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*planet2_mass, planet2_position, planet2_velocity)
    planet2_inclination = planet2_keplerian_orbital_elements[3]
    planet2_spin = calculate_spin(planet2_angular_frequency, planet2_inclination, planet2_obliquity, planet2_position, planet2_velocity)

    planet2_evolution_type = posidonius.NonEvolving()
    #universe.add_particle(planet2_mass, planet2_radius, planet2_dissipation_factor, planet2_dissipation_factor_scale, planet2_radius_of_gyration_2, planet2_love_number, planet2_fluid_love_number, planet2_position, planet2_velocity, planet2_spin, planet2_evolution_type)
    universe.add_earth_like(planet2_mass, planet2_dissipation_factor_scale, planet2_position, planet2_velocity, planet2_spin, planet2_evolution_type)

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    #universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


