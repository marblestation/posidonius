import posidonius

def _posvel(star_mass, planet_mass, a):
    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    #a = 0.018;                             # semi-major axis (in AU)
    e = 0.1;                               # eccentricity
    i = 5. * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;                      # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    gm = posidonius.constants.G*(planet_mass+star_mass);
    x, y, z, vx, vy, vz = posidonius.calculate_cartesian_coordinates(gm, q, e, i, p, n, l);
    position = posidonius.Axes(x, y, z)
    velocity = posidonius.Axes(vx, vy, vz)
    return position, velocity

def _spin(obliquity, rotation_period, star_mass, planet_mass, position, velocity):
    #////// Initialization of planetary spin
    angular_frequency = posidonius.constants.TWO_PI/(rotation_period/24.) # days^-1
    keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(posidonius.constants.G*star_mass*planet_mass, position, velocity)
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
    planet1_semimajoraxis = 0.018
    planet1_position, planet1_velocity = _posvel(star_mass, planet1_mass, planet1_semimajoraxis)

    planet1_obliquity = 11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    planet1_rotation_period = 24. # hours
    planet1_spin = _spin(planet1_obliquity, planet1_rotation_period, star_mass, planet1_mass, planet1_position, planet1_velocity)

    planet1_evolution_type = posidonius.NonEvolving()
    universe.add_earth_like(planet1_mass, planet1_dissipation_factor_scale, planet1_position, planet1_velocity, planet1_spin, planet1_evolution_type)

    ############################################################################
    planet2_mass_factor = 0.00095
    planet2_mass = planet2_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    planet2_dissipation_factor_scale = 1.0
    planet2_semimajoraxis = 0.18
    planet2_position, planet2_velocity = _posvel(star_mass, planet2_mass, planet2_semimajoraxis)

    planet2_obliquity = 11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    planet2_rotation_period = 24. # hours
    planet2_spin = _spin(planet2_obliquity, planet2_rotation_period, star_mass, planet2_mass, planet2_position, planet2_velocity)

    planet2_evolution_type = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    universe.add_jupiter_like(planet2_mass, planet2_dissipation_factor_scale, planet2_position, planet2_velocity, planet2_spin, planet2_evolution_type)

    ############################################################################
    planet3_mass_factor = 0.00095
    planet3_mass = planet3_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    planet3_dissipation_factor_scale = 1.0
    planet3_semimajoraxis = 1.8
    planet3_position, planet3_velocity = _posvel(star_mass, planet3_mass, planet3_semimajoraxis)

    planet3_obliquity = 11.459156 * posidonius.constants.DEG2RAD # 0.2 rad
    planet3_rotation_period = 24. # hours
    planet3_spin = _spin(planet3_obliquity, planet3_rotation_period, star_mass, planet3_mass, planet3_position, planet3_velocity)

    planet3_evolution_type = posidonius.NonEvolving()
    universe.add_jupiter_like(planet3_mass, planet3_dissipation_factor_scale, planet3_position, planet3_velocity, planet3_spin, planet3_evolution_type)
