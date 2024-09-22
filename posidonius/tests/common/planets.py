import posidonius
import numpy as np

def _posvel(star_mass, planet_mass, a, e):
    zeros = posidonius.Axes(0., 0., 0.)
    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    #a = 0.018;                             # semi-major axis (in AU)
    #e = 0.1;                               # eccentricity
    i = 5. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    position, velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[zeros], velocities=[zeros])
    return position, velocity

def _spin(obliquity, rotation_period, star_mass, planet_mass, planet_position, planet_velocity):
    zeros = posidonius.Axes(0., 0., 0.)
    #////// Initialization of planetary spin
    angular_frequency = posidonius.constants.TWO_PI/(rotation_period/24.) # days^-1
    keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[zeros], velocities=[zeros])
    inclination = keplerian_orbital_elements[3]
    spin = posidonius.calculate_spin(angular_frequency, inclination, obliquity)
    return spin

def basic_configuration(universe):
    if len(universe._data['particles']) == 0:
        raise Exception("Star needs to be already present!")
    star_mass = universe._data['particles'][0]['mass']
    ############################################################################
    planet1_mass_factor = 1.0
    planet1_mass = planet1_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    planet1_dissipation_factor_scale = 1.0
    planet1_semimajoraxis = 0.5
    planet1_eccentricity = 0.1
    planet1_position, planet1_velocity = _posvel(star_mass, planet1_mass, planet1_semimajoraxis, planet1_eccentricity)

    planet1_obliquity = 11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    planet1_rotation_period = 24. # hours
    planet1_spin = _spin(planet1_obliquity, planet1_rotation_period, star_mass, planet1_mass, planet1_position, planet1_velocity)

    planet1_evolution = posidonius.NonEvolving()
    planet1 = earth_like(planet1_mass, planet1_dissipation_factor_scale, planet1_position, planet1_velocity, planet1_spin, planet1_evolution)
    universe.add_particle(planet1)

    ############################################################################
    planet2_mass_factor = 0.00095
    planet2_mass = planet2_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    planet2_dissipation_factor_scale = 1.0
    planet2_semimajoraxis = 5.0
    planet2_eccentricity = 0.1
    planet2_position, planet2_velocity = _posvel(star_mass, planet2_mass, planet2_semimajoraxis, planet2_eccentricity)

    planet2_obliquity = 45.0 * posidonius.constants.DEG2RAD
    planet2_rotation_period = 24. # hours
    planet2_spin = _spin(planet2_obliquity, planet2_rotation_period, star_mass, planet2_mass, planet2_position, planet2_velocity)

    planet2_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    planet2 = jupiter_like(planet2_mass, planet2_dissipation_factor_scale, planet2_position, planet2_velocity, planet2_spin, planet2_evolution)
    universe.add_particle(planet2)

    ############################################################################
    planet3_mass_factor = 0.00095
    planet3_mass = planet3_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    planet3_dissipation_factor_scale = 1.0
    planet3_semimajoraxis = 10.0
    planet3_eccentricity = 0.1 # this helps testing Kaula special case for small obliquities & eccentricity different than zero
    planet3_position, planet3_velocity = _posvel(star_mass, planet3_mass, planet3_semimajoraxis, planet3_eccentricity)

    planet3_obliquity = 0. * posidonius.constants.DEG2RAD # 0.0 rad => this helps testing Kaula special case for small obliquities
    planet3_rotation_period = 24. # hours
    planet3_spin = _spin(planet3_obliquity, planet3_rotation_period, star_mass, planet3_mass, planet3_position, planet3_velocity)

    planet3_evolution = posidonius.NonEvolving()
    planet3 = jupiter_like(planet3_mass, planet3_dissipation_factor_scale, planet3_position, planet3_velocity, planet3_spin, planet3_evolution)
    universe.add_particle(planet3)

    ############################################################################
    planet4_mass_factor = 0.00095
    planet4_mass = planet4_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    planet4_dissipation_factor_scale = 1.0
    planet4_semimajoraxis = 15.0
    planet4_eccentricity = 0.0 # this helps testing Kaula special case for small obliquities & eccentricity exactly zero
    planet4_position, planet4_velocity = _posvel(star_mass, planet4_mass, planet4_semimajoraxis, planet4_eccentricity)

    planet4_obliquity = 0. * posidonius.constants.DEG2RAD # 0.0 rad => this helps testing Kaula special case for small obliquities
    planet4_rotation_period = 24. # hours
    planet4_spin = _spin(planet4_obliquity, planet4_rotation_period, star_mass, planet4_mass, planet4_position, planet4_velocity)

    planet4_evolution = posidonius.NonEvolving()
    planet4 = jupiter_like(planet4_mass, planet4_dissipation_factor_scale, planet4_position, planet4_velocity, planet4_spin, planet4_evolution)
    universe.add_particle(planet4)


def earth_like(mass, dissipation_factor_scale, position, velocity, spin, evolution):
    if type(evolution) not in (posidonius.NonEvolving,):
        raise Exception("Evolution type should be NonEvolving!")


    # Earth-like => mass-radius relationship from Fortney 2007
    radius_factor = posidonius.tools.mass_radius_relation(mass, planet_mass_type='AU', planet_percent_rock=0.70)
    radius = radius_factor * posidonius.constants.R_EARTH
    radius_of_gyration = 5.75e-01 # Earth type planet

    # Typical rotation period: 24 hours
    love_number = 0.299 # Earth
    fluid_love_number = 0.9532 # Earth
    k2pdelta = 2.465278e-3 # Terrestrial planets
    dissipation_factor = 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(radius, 5))
    tides = posidonius.effects.tides.OrbitingBody(
                posidonius.effects.tides.ConstantTimeLag({
                    "dissipation_factor_scale": dissipation_factor_scale,
                    "dissipation_factor": dissipation_factor,
                    "love_number": love_number,
                })
            )
    rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(
                                posidonius.effects.rotational_flattening.OblateSpheroid({
                                    "love_number": fluid_love_number
                                })
                            )
    general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    wind = posidonius.effects.wind.Disabled()
    disk = posidonius.effects.disk.OrbitingBody()
    particle = posidonius.Particle(mass, radius, radius_of_gyration, position, velocity, spin)
    particle.set_tides(tides)
    particle.set_rotational_flattening(rotational_flattening)
    particle.set_general_relativity(general_relativity)
    particle.set_wind(wind)
    particle.set_disk(disk)
    particle.set_evolution(evolution)
    return particle


def jupiter_like(mass, dissipation_factor_scale, position, velocity, spin, evolution):
    if type(evolution) not in (posidonius.LeconteChabrier2013, posidonius.NonEvolving):
        raise Exception("Evolution type should be LeconteChabrier2013 or NonEvolving!")

    radius_factor = 10.9 # Jupiter in R_EARTH
    radius = radius_factor * posidonius.constants.R_EARTH
    radius_of_gyration = 5.04e-01 # Gas giant

    # Typical rotation period: 9.8 hours
    love_number = 0.380 # Gas giant
    # TODO: What k2pdelta/dissipation_factor is the recommended?
    #k2pdelta = 8.101852e-9 # Gas giant
    k2pdelta = 2.893519e-7 # Gas giant for Jupiter: 2-3d-2 s, here in day (Leconte)
    dissipation_factor = 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(radius, 5))
    #dissipation_factor = 2.006*3.845764e4 // Gas giant
    tides = posidonius.effects.tides.OrbitingBody(
                posidonius.effects.tides.ConstantTimeLag({
                    "dissipation_factor_scale": dissipation_factor_scale,
                    "dissipation_factor": dissipation_factor,
                    "love_number": love_number,
                })
            )
    rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(
                                posidonius.effects.rotational_flattening.OblateSpheroid({
                                    "love_number": love_number
                                })
                            )
    general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    wind = posidonius.effects.wind.Disabled()
    disk = posidonius.effects.disk.OrbitingBody()
    particle = posidonius.Particle(mass, radius, radius_of_gyration, position, velocity, spin)
    particle.set_tides(tides)
    particle.set_rotational_flattening(rotational_flattening)
    particle.set_general_relativity(general_relativity)
    particle.set_wind(wind)
    particle.set_disk(disk)
    particle.set_evolution(evolution)
    return particle
