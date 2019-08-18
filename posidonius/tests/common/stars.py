import posidonius

def solar_like(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842):
    """
    Wind parametrisation (Bouvier 1997):

    wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    """
    disk = posidonius.effects.disk.Disabled()
    return _solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution, disk, wind_k_factor, wind_rotation_saturation)

def solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842):
    """
    Wind parametrisation (Bouvier 1997):

    wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    """
    disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    disk_properties = {
        'inner_edge_distance': 0.01,  # AU
        'outer_edge_distance': 100.0, # AU
        'lifetime': 1.0e5 * 365.25e0, # days
        'alpha': 1.0e-2,
        'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        'mean_molecular_weight': 2.4,
    }
    disk = posidonius.effects.disk.CentralBody(disk_properties)
    return _solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution, disk, wind_k_factor, wind_rotation_saturation)

def _solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution, disk, wind_k_factor, wind_rotation_saturation):
    if type(evolution) not in (posidonius.BolmontMathis2016, posidonius.Leconte2011, posidonius.Baraffe2015, posidonius.GalletBolmont2017, posidonius.NonEvolving) and not (type(evolution) is posidonius.Baraffe1998 and mass == 1.0):
        raise Exception("Evolution type should be BolmontMathis2016 Leconte2011 Baraffe1998 (mass = 0.10) Baraffe2015 GalletBolmont2017 or NonEvolving!")

    # Typical rotation period: 24 hours
    angular_frequency = posidonius.constants.TWO_PI/(rotation_period/24.) # days^-1
    inclination = 0.
    obliquity = 0.
    spin = posidonius.tools.calculate_spin(angular_frequency, inclination, obliquity)

    radius_factor = 1.
    radius = radius_factor * posidonius.constants.R_SUN
    radius_of_gyration = 2.43e-01 # Sun

    love_number = 0.03
    # Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    dissipation_factor = 4.992*3.845764e-2 # -66+64
    tides = posidonius.effects.tides.CentralBody({
        "dissipation_factor_scale": dissipation_factor_scale,
        "dissipation_factor": dissipation_factor,
        "love_number": love_number,
    })
    rotational_flattening = posidonius.effects.rotational_flattening.CentralBody({"love_number": love_number})
    general_relativity = posidonius.effects.general_relativity.CentralBody(general_relativity_implementation)
    if wind_k_factor == 0:
        wind = posidonius.effects.wind.Disabled()
    else:
        wind = posidonius.effects.wind.Interaction({
            "k_factor": wind_k_factor,
            "rotation_saturation": wind_rotation_saturation,
        })
    particle = posidonius.Particle(mass, radius, radius_of_gyration, position, velocity, spin)
    particle.set_tides(tides)
    particle.set_rotational_flattening(rotational_flattening)
    particle.set_general_relativity(general_relativity)
    particle.set_wind(wind)
    particle.set_disk(disk)
    particle.set_evolution(evolution)
    return particle


def m_dwarf(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution, wind_k_factor=0., wind_rotation_saturation=0.):
    if type(evolution) not in (posidonius.Baraffe2015, posidonius.NonEvolving) and not (type(evolution) is posidonius.Baraffe1998 and mass == 0.10):
        raise Exception("Evolution type should be Baraffe2015 Baraffe1998 (mass = 0.10) or NonEvolving!")

    # Typical rotation period: 70 hours
    angular_frequency = posidonius.constants.TWO_PI/(rotation_period/24.) # days^-1
    inclination = 0.
    obliquity = 0.
    spin = posidonius.tools.calculate_spin(angular_frequency, inclination, obliquity)

    radius_factor = 0.845649342247916
    radius = radius_factor * posidonius.constants.R_SUN
    radius_of_gyration = 4.47e-01 # M-dwarf

    love_number = 0.307 # M Dwarf
    # BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    dissipation_factor = 2.006*3.845764e4 # -60+64
    tides = posidonius.effects.tides.CentralBody({
        "dissipation_factor_scale": dissipation_factor_scale,
        "dissipation_factor": dissipation_factor,
        "love_number": love_number,
    })
    rotational_flattening = posidonius.effects.rotational_flattening.CentralBody({"love_number": love_number})
    general_relativity = posidonius.effects.general_relativity.CentralBody(general_relativity_implementation)
    if wind_k_factor == 0:
        wind = posidonius.effects.wind.Disabled()
    else:
        wind = posidonius.effects.wind.Interaction({
            "k_factor": wind_k_factor,
            "rotation_saturation": wind_rotation_saturation,
        })
    disk = posidonius.effects.disk.Disabled()
    particle = posidonius.Particle(mass, radius, radius_of_gyration, position, velocity, spin)
    particle.set_tides(tides)
    particle.set_rotational_flattening(rotational_flattening)
    particle.set_general_relativity(general_relativity)
    particle.set_wind(wind)
    particle.set_disk(disk)
    particle.set_evolution(evolution)
    return particle

def brown_dwarf(mass, dissipation_factor_scale, position, velocity, general_relativity_implementation, evolution, wind_k_factor=0., wind_rotation_saturation=0.):
    rotation_period = None
    love_number = None
    if type(evolution) == posidonius.NonEvolving:
        rotation_period = 70.0 # hours
        love_number = 0.307 # BrownDwarf
    elif type(evolution) in (posidonius.Leconte2011, posidonius.Baraffe2015):
        mass = evolution._data[evolution.__class__.__name__]

        if mass <= 0.0101 and mass >= 0.0099:
            rotation_period = 8.0
            love_number = 0.3790
        elif mass <= 0.0121 and mass >= 0.0119:
            rotation_period = 13.0
            love_number = 0.3780
        elif mass <= 0.0151 and mass >= 0.0149:
            rotation_period = 19.0
            love_number = 0.3760
        elif mass <= 0.0201 and mass >= 0.0199:
            rotation_period = 24.0
            love_number = 0.3690
        elif mass <= 0.0301 and mass >= 0.0299:
            rotation_period = 30.0
            love_number = 0.3550
        elif mass <= 0.0401 and mass >= 0.0399:
            rotation_period = 36.0
            love_number = 0.3420
        elif mass <= 0.0501 and mass >= 0.0499:
            rotation_period = 41.0
            love_number = 0.3330
        elif mass <= 0.0601 and mass >= 0.0599:
            rotation_period = 47.0
            love_number = 0.3250
        elif mass <= 0.0701 and mass >= 0.0699:
            rotation_period = 53.0
            love_number = 0.3110
        elif mass <= 0.0721 and mass >= 0.0719:
            rotation_period = 58.0
            love_number = 0.3080
        elif mass <= 0.0751 and mass >= 0.0749:
            rotation_period = 64.0
            love_number = 0.3070
        elif mass <= 0.0801 and mass >= 0.0799:
            rotation_period = 70.0
            love_number = 0.3070
        else:
            raise Exception("The evolution type Leconte2011 does not support a mass of {} Msun and Baraffe2015 with higher masses is possible but then it should not be considered as a Brown Dwarf!".format(mass))
    else:
        raise Exception("Evolution type should be Leconte2011, Baraffe2015 or NonEvolving!")

    angular_frequency = posidonius.constants.TWO_PI/(rotation_period/24.) # days^-1
    inclination = 0.
    obliquity = 0.
    spin = posidonius.tools.calculate_spin(angular_frequency, inclination, obliquity)

    # BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    dissipation_factor = 2.006*3.845764e4 # -60+64

    radius_factor = 0.845649342247916
    radius = radius_factor * posidonius.constants.R_SUN
    radius_of_gyration = 4.41e-01 # Brown dwarf

    tides = posidonius.effects.tides.CentralBody({
        "dissipation_factor_scale": dissipation_factor_scale,
        "dissipation_factor": dissipation_factor,
        "love_number": love_number,
    })
    rotational_flattening = posidonius.effects.rotational_flattening.CentralBody({"love_number": love_number})
    general_relativity = posidonius.effects.general_relativity.CentralBody(general_relativity_implementation)
    if wind_k_factor == 0:
        wind = posidonius.effects.wind.Disabled()
    else:
        wind = posidonius.effects.wind.Interaction({
            "k_factor": wind_k_factor,
            "rotation_saturation": wind_rotation_saturation,
        })
    disk = posidonius.effects.disk.Disabled()
    particle = posidonius.Particle(mass, radius, radius_of_gyration, position, velocity, spin)
    particle.set_tides(tides)
    particle.set_rotational_flattening(rotational_flattening)
    particle.set_general_relativity(general_relativity)
    particle.set_wind(wind)
    particle.set_disk(disk)
    particle.set_evolution(evolution)
    return particle
