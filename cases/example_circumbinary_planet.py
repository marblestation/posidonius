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
    time_step = 20.0 # days
    #time_limit   = 4*time_step # days
    time_limit   = 365.25 * 1.0e8 # days
    historic_snapshot_period = 100.*365.25 # days
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": True,
        "rotational_flattening": True,
        "general_relativity": True,
        "disk": False,
        "wind": False,
        "evolution": True,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star1_mass = 10.0 # Solar masses
    star1_radius_factor = 1.0
    star1_radius = star1_radius_factor * posidonius.constants.R_SUN
    star1_radius_of_gyration = 2.43e-01 # Sun

    star1_position = posidonius.Axes(0., 0., 0.)
    star1_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star1_rotation_period = 24.0 # hours
    star1_angular_frequency = posidonius.constants.TWO_PI/(star1_rotation_period/24.) # days^-1
    star1_spin = posidonius.Axes(0., 0., star1_angular_frequency)

    star1_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 4.992*3.845764e-2,
        "love_number": 0.03,
    }
    star1_tides_model = posidonius.effects.tides.ConstantTimeLag(star1_tides_parameters)
    star1_tides = posidonius.effects.tides.CentralBody(star1_tides_model)
    #star1_tides = posidonius.effects.tides.OrbitingBody(star1_tides_model)
    #star1_tides = posidonius.effects.tides.Disabled()
    #
    star1_rotational_flattening_parameters = {"love_number": star1_tides_parameters["love_number"] }
    star1_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star1_rotational_flattening_parameters)
    star1_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star1_rotational_flattening_model)
    #star1_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(star1_rotational_flattening_model)
    #star1_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    star1_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #star1_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #star1_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    #star1_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    #star1_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    #star1_wind = posidonius.effects.wind.Interaction({
        ## Solar wind parametrisation (Bouvier 1997)
        #"k_factor": 4.0e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        #"rotation_saturation": 1.7592918860102842, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    #})
    star1_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
        #'inner_edge_distance': 0.01,  # AU
        #'outer_edge_distance': 100.0, # AU
        #'lifetime': 1.0e5 * 365.25e0, # days
        #'alpha': 1.0e-2,
        #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        #'mean_molecular_weight': 2.4,
    #}
    #star1_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #star1_disk = posidonius.effects.disk.OrbitingBody()
    star1_disk = posidonius.effects.disk.Disabled()
    #
    #star1_evolution = posidonius.GalletBolmont2017(star1_mass) # mass = 0.30 .. 1.40
    #star1_evolution = posidonius.BolmontMathis2016(star1_mass) # mass = 0.40 .. 1.40
    #star1_evolution = posidonius.Baraffe2015(star1_mass) # mass = 0.01 .. 1.40
    #star1_evolution = posidonius.Leconte2011(star1_mass) # mass = 0.01 .. 0.08
    #star1_evolution = posidonius.Baraffe1998(star1_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star1_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #star1_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
    star1_evolution = posidonius.NonEvolving()
    #
    star1 = posidonius.Particle(star1_mass, star1_radius, star1_radius_of_gyration, star1_position, star1_velocity, star1_spin)
    star1.set_tides(star1_tides)
    star1.set_rotational_flattening(star1_rotational_flattening)
    star1.set_general_relativity(star1_general_relativity)
    star1.set_wind(star1_wind)
    star1.set_disk(star1_disk)
    star1.set_evolution(star1_evolution)
    universe.add_particle(star1)



    ############################################################################
    star2_mass = 1.0 # Solar masses
    star2_radius_factor = 0.845649342247916
    star2_radius = star2_radius_factor * posidonius.constants.R_SUN
    star2_radius_of_gyration = 4.41e-01 # Brown dwarf

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 5.0;                             # semi-major axis (in AU)
    e = 0.1;                               # eccentricity
    i = 0.5 * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    star2_position, star2_velocity = posidonius.calculate_cartesian_coordinates(star2_mass, q, e, i, p, n, l, masses=[star1_mass], positions=[star1_position], velocities=[star1_velocity])

    #////// Initialization of star2ary spin
    star2_obliquity = 5.0 * posidonius.constants.DEG2RAD # 0.2 rad
    star2_rotation_period = 70. # hours
    star2_angular_frequency = posidonius.constants.TWO_PI/(star2_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #star2_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(star2_mass, star2_position, star2_velocity, masses=[star1_mass], positions=[star1_position], velocities=[star1_velocity])
    #star2_semi_major_axis = star2_keplerian_orbital_elements[0]
    #star2_eccentricity = star2_keplerian_orbital_elements[2]
    #star2_semi_major_axis = a
    #star2_eccentricity = e
    #star2_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(star2_semi_major_axis, star2_eccentricity, star1_mass, star2_mass) # days
    #star2_angular_frequency = posidonius.constants.TWO_PI/(star2_pseudo_synchronization_period) # days^-1
    star2_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(star2_mass, star2_position, star2_velocity, masses=[star1_mass], positions=[star1_position], velocities=[star1_velocity])
    star2_inclination = star2_keplerian_orbital_elements[3]
    star2_spin = posidonius.calculate_spin(star2_angular_frequency, star2_inclination, star2_obliquity)

    star2_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 2.006*3.845764e4,
        "love_number": 0.307,
    }
    star2_tides_model = posidonius.effects.tides.ConstantTimeLag(star2_tides_parameters)
    #star2_tides = posidonius.effects.tides.CentralBody(star2_tides_model)
    star2_tides = posidonius.effects.tides.OrbitingBody(star2_tides_model)
    #star2_tides = posidonius.effects.tides.Disabled()
    #
    star2_rotational_flattening_parameters = {"love_number": star2_tides_parameters["love_number"] }
    star2_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star2_rotational_flattening_parameters)
    #star2_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star2_rotational_flattening_model)
    star2_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(star2_rotational_flattening_model)
    #star2_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    #star2_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #star2_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #star2_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    star2_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    #star2_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    #star2_wind = posidonius.effects.wind.Interaction({
        ## Solar wind parametrisation (Bouvier 1997)
        #"k_factor": 4.0e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        #"rotation_saturation": 1.7592918860102842, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    #})
    star2_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
        #'inner_edge_distance': 0.01,  # AU
        #'outer_edge_distance': 100.0, # AU
        #'lifetime': 1.0e5 * 365.25e0, # days
        #'alpha': 1.0e-2,
        #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        #'mean_molecular_weight': 2.4,
    #}
    #star2_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #star2_disk = posidonius.effects.disk.OrbitingBody()
    star2_disk = posidonius.effects.disk.Disabled()
    #
    #star2_evolution = posidonius.GalletBolmont2017(star2_mass) # mass = 0.30 .. 1.40
    #star2_evolution = posidonius.BolmontMathis2016(star2_mass) # mass = 0.40 .. 1.40
    star2_evolution = posidonius.Baraffe2015(star2_mass) # mass = 0.01 .. 1.40
    #star2_evolution = posidonius.Leconte2011(star2_mass) # mass = 0.01 .. 0.08
    #star2_evolution = posidonius.Baraffe1998(star2_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star2_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #star2_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
    #star2_evolution = posidonius.NonEvolving()
    #
    star2 = posidonius.Particle(star2_mass, star2_radius, star2_radius_of_gyration, star2_position, star2_velocity, star2_spin)
    star2.set_tides(star2_tides)
    star2.set_rotational_flattening(star2_rotational_flattening)
    star2.set_general_relativity(star2_general_relativity)
    star2.set_wind(star2_wind)
    star2.set_disk(star2_disk)
    star2.set_evolution(star2_evolution)
    universe.add_particle(star2)

    ############################################################################
    planet_mass_factor = 1.0
    planet_mass = planet_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    # Earth-like => mass-radius relationship from Fortney 2007
    planet_radius_factor = posidonius.tools.mass_radius_relation(planet_mass_factor, planet_mass_type='factor', planet_percent_rock=0.70)
    planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
    planet_radius_of_gyration = 5.75e-01 # Earth type planet

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 40
    e = 0.;                               # eccentricity
    i = 0. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star1_mass, star2_mass], positions=[star1_position, star2_position], velocities=[star1_velocity, star2_velocity])


    #////// Initialization of planetary spin
    planet_obliquity = 11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    planet_rotation_period = 24. # hours
    planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period/24.) # days^-1
    # Pseudo-synchronization period
    #planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star1_mass, star2_mass], positions=[star1_position, star2_position], velocities=[star2_velocity, star2_velocity])
    #planet_semi_major_axis = planet_keplerian_orbital_elements[0]
    #planet_eccentricity = planet_keplerian_orbital_elements[2]
    #planet_semi_major_axis = a
    #planet_eccentricity = e
    #planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(planet_semi_major_axis, planet_eccentricity, star1_mass+star2_mass, planet_mass)
    #planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period/24.) # days^-1
    planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star1_mass, star2_mass], positions=[star1_position, star2_position], velocities=[star2_velocity, star2_velocity])
    planet_inclination = planet_keplerian_orbital_elements[3]
    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

    k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
    planet_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        "love_number": 0.305,
    }
    planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)
    #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_model)
    #planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
    planet_tides = posidonius.effects.tides.Disabled()
    #
    planet_rotational_flattening_parameters = {"love_number": planet_tides_parameters["love_number"]}
    planet_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(planet_rotational_flattening_parameters)
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(planet_rotational_flattening_model)
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_model)
    planet_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    #planet_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    planet_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    #planet_wind = posidonius.effects.wind.Interaction({
        ## Solar wind parametrisation (Bouvier 1997)
        #"k_factor": 4.0e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        #"rotation_saturation": 1.7592918860102842, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    #})
    planet_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
        #'inner_edge_distance': 0.01,  # AU
        #'outer_edge_distance': 100.0, # AU
        #'lifetime': 1.0e5 * 365.25e0, # days
        #'alpha': 1.0e-2,
        #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        #'mean_molecular_weight': 2.4,
    #}
    #planet_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #planet_disk = posidonius.effects.disk.OrbitingBody()
    planet_disk = posidonius.effects.disk.Disabled()
    #
    #planet_evolution = posidonius.GalletBolmont2017(planet_mass) # mass = 0.30 .. 1.40
    #planet_evolution = posidonius.BolmontMathis2016(planet_mass) # mass = 0.40 .. 1.40
    #planet_evolution = posidonius.Baraffe2015(planet_mass) # mass = 0.01 .. 1.40
    #planet_evolution = posidonius.Leconte2011(planet_mass) # mass = 0.01 .. 0.08
    #planet_evolution = posidonius.Baraffe1998(planet_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #planet_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #planet_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
    planet_evolution = posidonius.NonEvolving()
    #
    planet = posidonius.Particle(planet_mass, planet_radius, planet_radius_of_gyration, planet_position, planet_velocity, planet_spin)
    planet.set_tides(planet_tides)
    planet.set_rotational_flattening(planet_rotational_flattening)
    planet.set_general_relativity(planet_general_relativity)
    planet.set_wind(planet_wind)
    planet.set_disk(planet_disk)
    planet.set_evolution(planet_evolution)
    universe.add_particle(planet)


    #whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    #universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


