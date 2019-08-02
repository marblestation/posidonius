extern crate posidonius;

pub fn solar_like(star_mass: f64, star_evolution_type: posidonius::EvolutionType) -> posidonius::Particle {
    match star_evolution_type {
        posidonius::EvolutionType::BolmontMathis2016(_) => { },
        posidonius::EvolutionType::GalletBolmont2017(_) => { },
        posidonius::EvolutionType::Baraffe1998(_) => { },
        posidonius::EvolutionType::Baraffe2015(_) => { },
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"); }
    }
    let radius_factor: f64 = 1.0;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    let star_love_number: f64 = 0.03;  // Sun
    let star_fluid_love_number: f64 = star_love_number;
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 4.992*3.845764e-2; // -66+64
    // Radius of gyration
    let star_radius_of_gyration_2: f64 = 5.9e-2; // Sun
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes{x:0., y:0., z:0.};
    let star_velocity = posidonius::Axes{x:0., y:0., z:0.};
    let star_acceleration = posidonius::Axes{x:0., y:0., z:0.};
    // Initialization of stellar spin
    let star_rotation_period: f64 = 8.; // hours
    let star_angular_frequency = posidonius::constants::TWO_PI/(star_rotation_period/24.); // days^-1
    let star_spin = posidonius::Axes{x:0., y:0., z:star_angular_frequency };
    // Disk
    let star_disk = posidonius::Disk::None;

    // Wind parametrisation (Bouvier 1997):
    //      wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    //      wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    let star_wind_k_factor = 4.0e-18;
    let star_wind_rotation_saturation = 1.7592918860102842;
    let star = posidonius::Particle::new(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, 
                                            star_love_number, star_fluid_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin,
                                            star_evolution_type, star_disk, star_wind_k_factor, star_wind_rotation_saturation);
    star
}

pub fn solar_like_with_disk(star_mass: f64, star_evolution_type: posidonius::EvolutionType) -> posidonius::Particle {
    match star_evolution_type {
        posidonius::EvolutionType::BolmontMathis2016(_) => { },
        posidonius::EvolutionType::GalletBolmont2017(_) => { },
        posidonius::EvolutionType::Baraffe1998(_) => { },
        posidonius::EvolutionType::Baraffe2015(_) => { },
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"); }
    }
    let radius_factor: f64 = 1.0;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    let star_love_number: f64 = 0.03;  // Sun
    let star_fluid_love_number: f64 = star_love_number;
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 4.992*3.845764e-2; // -66+64
    // Radius of gyration
    let star_radius_of_gyration_2: f64 = 5.9e-2; // Sun
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes{x:0., y:0., z:0.};
    let star_velocity = posidonius::Axes{x:0., y:0., z:0.};
    let star_acceleration = posidonius::Axes{x:0., y:0., z:0.};
    // Initialization of stellar spin
    let star_rotation_period: f64 = 8.; // hours
    let star_angular_frequency = posidonius::constants::TWO_PI/(star_rotation_period/24.); // days^-1
    let star_spin = posidonius::Axes{x:0., y:0., z:star_angular_frequency };
    // Disk
    let surface_density_normalization_gcm = 1000.; // g.cm^-2
    let surface_density_normalization_si = surface_density_normalization_gcm * 1.0e-3 * 1.0e4; // kg.m^-2
    let star_disk = posidonius::Disk::Host(posidonius::DiskProperties {
        inner_edge_distance: 0.01,  // AU
        outer_edge_distance: 100.0, // AU
        lifetime: 1.0e5 * 365.25e0, // days
        alpha: 1.0e-2,
        surface_density_normalization: surface_density_normalization_si * (1.0/posidonius::constants::M_SUN) * posidonius::constants::AU.powi(2), // Msun.AU^-2
        mean_molecular_weight: 2.4,
    });

    // Wind parametrisation (Bouvier 1997):
    //      wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    //      wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    let star_wind_k_factor = 4.0e-18;
    let star_wind_rotation_saturation = 1.7592918860102842;
    let star = posidonius::Particle::new(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, 
                                            star_love_number, star_fluid_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin,
                                            star_evolution_type, star_disk, star_wind_k_factor, star_wind_rotation_saturation);
    star
}


pub fn brown_dwarf(star_mass: f64, star_evolution_type: posidonius::EvolutionType) -> posidonius::Particle {
    match star_evolution_type {
        posidonius::EvolutionType::Leconte2011(_) => { },
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be Leconte2011 or non evolving to create a Brown Dwarf body!"); }
    }

    let (star_rotation_period, star_love_number) = match star_evolution_type {
        posidonius::EvolutionType::NonEvolving => { 
            let rotation_period: f64 = 70.0; // hours
            let love_number: f64 = 0.307; // BrownDwarf
            (rotation_period, love_number)
        },
        posidonius::EvolutionType::Leconte2011(mass) => { 
            let rotation_period; // hours
            let love_number;
            if mass <= 0.0101 && mass >= 0.0099 {
                rotation_period = 8.0;
                love_number = 0.3790;
            } else if mass <= 0.0121 && mass >= 0.0119 {
                rotation_period = 13.0;
                love_number = 0.3780;
            } else if mass <= 0.0151 && mass >= 0.0149 {
                rotation_period = 19.0;
                love_number = 0.3760;
            } else if mass <= 0.0201 && mass >= 0.0199 {
                rotation_period = 24.0;
                love_number = 0.3690;
            } else if mass <= 0.0301 && mass >= 0.0299 {
                rotation_period = 30.0;
                love_number = 0.3550;
            } else if mass <= 0.0401 && mass >= 0.0399 {
                rotation_period = 36.0;
                love_number = 0.3420;
            } else if mass <= 0.0501 && mass >= 0.0499 {
                rotation_period = 41.0;
                love_number = 0.3330;
            } else if mass <= 0.0601 && mass >= 0.0599 {
                rotation_period = 47.0;
                love_number = 0.3250;
            } else if mass <= 0.0701 && mass >= 0.0699 {
                rotation_period = 53.0;
                love_number = 0.3110;
            } else if mass <= 0.0721 && mass >= 0.0719 {
                rotation_period = 58.0;
                love_number = 0.3080;
            } else if mass <= 0.0751 && mass >= 0.0749 {
                rotation_period = 64.0;
                love_number = 0.3070;
            } else if mass <= 0.0801 && mass >= 0.0799 {
                rotation_period = 70.0;
                love_number = 0.3070;
            } else {
                panic!("The evolution type Leconte2011 does not support a mass of {} Msun!", mass);
            }
            (rotation_period, love_number)
        },
        _ => { panic!("Evolution type should be solar like or non evolving to create a solar like body!"); }
    };
    let star_fluid_love_number: f64 = star_love_number;

    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 2.006*3.845764e4; // -60+64
    // Radius of gyration
    let star_radius_of_gyration_2: f64 = 1.94e-1; // Brown dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes{x:0., y:0., z:0.};
    let star_velocity = posidonius::Axes{x:0., y:0., z:0.};
    let star_acceleration = posidonius::Axes{x:0., y:0., z:0.};
    // Initialization of stellar spin
    let star_angular_frequency = posidonius::constants::TWO_PI/(star_rotation_period/24.); // days^-1
    let star_spin = posidonius::Axes{x:0., y:0., z:star_angular_frequency };
    // Disk
    let star_disk = posidonius::Disk::None;

    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;

    let star = posidonius::Particle::new(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, 
                                            star_love_number, star_fluid_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin,
                                            star_evolution_type, star_disk, star_wind_k_factor, star_wind_rotation_saturation);
    star
}

pub fn m_dwarf(star_mass: f64, star_evolution_type: posidonius::EvolutionType) -> posidonius::Particle {
    match star_evolution_type {
        posidonius::EvolutionType::Baraffe1998(_) => { },
        posidonius::EvolutionType::Baraffe2015(_) => { },
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be Baraffe1998, Baraffe2015 or non evolving to create a M-Dwarf body!"); }
    }

    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    let star_love_number: f64 = 0.307; // Brown Dwarf / M Dwarf
    let star_fluid_love_number: f64 = star_love_number;
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 2.006*3.845764e4; // -60+64
    // Radius of gyration
    let star_radius_of_gyration_2: f64 = 2.0e-1; // M-dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes{x:0., y:0., z:0.};
    let star_velocity = posidonius::Axes{x:0., y:0., z:0.};
    let star_acceleration = posidonius::Axes{x:0., y:0., z:0.};
    // Initialization of stellar spin
    let star_rotation_period: f64 = 70.0; // hours
    let star_angular_frequency = posidonius::constants::TWO_PI/(star_rotation_period/24.); // days^-1
    let star_spin = posidonius::Axes{x:0., y:0., z:star_angular_frequency };
    // Disk
    let star_disk = posidonius::Disk::None;

    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;
    let star = posidonius::Particle::new(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, 
                                            star_love_number, star_fluid_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin,
                                            star_evolution_type, star_disk, star_wind_k_factor, star_wind_rotation_saturation);
    star
}

