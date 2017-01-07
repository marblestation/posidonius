use super::super::particles::{Axes, Particle, Universe};
use super::super::particles::{EvolutionType};
use super::super::constants::{R_SUN, TWO_PI, R_EARTH, K2, G, DEG2RAD, INTEGRATOR};
use super::super::tools::{calculate_cartesian_coordinates, calculate_keplerian_orbital_elements};
//use super::super::tools::{calculate_pseudo_synchronization_period};
use super::super::tools::{calculate_spin};



pub fn case3() -> Universe {
    ////////////////////////////////////////////////////////////////////////////
    // Initial conditions from CASE 3 in Bolmont et al. 2015
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    //---- Star
    let star_mass: f64 = 0.08; // Solar masses
    let star_dissipation_factor_scale: f64 = 1.;
    let star_evolution_type = EvolutionType::NonEvolving;
    let star_position = Axes{x:0., y:0., z:0.};
    let star_velocity = Axes{x:0., y:0., z:0.};
    let star_acceleration = Axes{x:0., y:0., z:0.};
    let star = Particle::new_brown_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_acceleration, star_evolution_type);

    ////////////////////////////////////////////////////////////////////////////
    //---- Planet
    let planet_radius_factor = 1.;
    let planet_mass: f64 = 3.0e-6; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet_dissipation_factor_scale: f64 = 1.;
    let planet_love_number: f64 = 0.305; // Earth
    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = 0.018;                             // semi-major axis (in AU)
    let e: f64 = 0.1;                               // eccentricity
    let i: f64 = 5. * DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0.;                            // argument of pericentre (degrees)
    let n: f64 = 0. * DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * DEG2RAD;                      // mean anomaly (degrees)
    p = (p + n) * DEG2RAD;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = G*(planet_mass+star_mass);
    let (x, y, z, vx, vy, vz) = calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    let planet_position = Axes{x:x, y:y, z:z};
    let planet_velocity = Axes{x:vx, y:vy, z:vz};
    let planet_acceleration = Axes{x:0., y:0., z:0.};
    
    ////// Initialization of planetary spin
    let planet_obliquity: f64 = 11.459156 * DEG2RAD; // 0.2 rad
    let planet_rotation_period: f64 = 24.; // hours
    //
    let planet_angular_frequency = TWO_PI/(planet_rotation_period/24.); // days^-1
    let planet_keplerian_orbital_elements = calculate_keplerian_orbital_elements(G*star_mass*planet_mass, planet_position, planet_velocity);
    let planet_inclination = planet_keplerian_orbital_elements.3;
    let planet_spin = calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_position, planet_velocity);

    let mut planet = Particle::new_terrestrial(planet_mass, planet_radius_factor, planet_dissipation_factor_scale, planet_position, planet_velocity, planet_acceleration);
    // Replace default values:
    planet.spin = planet_spin;
    planet.love_number = planet_love_number;
    ////////////////////////////////////////////////////////////////////////////

    let universe = Universe::new(vec![star, planet], INTEGRATOR);
    //println!("{:?}", universe);
    universe
}

pub fn case3_without_helpers() -> Universe {
    ////////////////////////////////////////////////////////////////////////////
    // Initial conditions from CASE 3 in Bolmont et al. 2015
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    //---- Star (central body)
    let star_mass: f64 = 0.08; // Solar masses
    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * R_SUN;
    let star_love_number: f64 = 0.307; // M Dwarf
    //let star_love_number: f64 = 0.03;  // Sun
    let star_fluid_love_number: f64 = star_love_number;
    ////// Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    //// Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    //let star_dissipation_factor: f64 = 4.992*3.845764e-2; // -66+64
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 2.006*3.845764e4; // -60+64
    ////// Radius of gyration
    //let star_radius_of_gyration_2: f64 = 5.9e-2; // Sun
    //let star_radius_of_gyration_2: f64 = 2.0e-1; // M-dwarf
    let star_radius_of_gyration_2: f64 = 1.94e-1; // Brown dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = Axes{x:0., y:0., z:0.};
    let star_velocity = Axes{x:0., y:0., z:0.};
    let star_acceleration = Axes{x:0., y:0., z:0.};
    ////// Initialization of stellar spin
    let star_rotation_period: f64 = 70.0; // hours
    let star_angular_frequency = TWO_PI/(star_rotation_period/24.); // days^-1
    let star_spin = Axes{x:0., y:0., z:star_angular_frequency };

    //let stellar_evolution_type = EvolutionType::BrownDwarf(star_mass);
    //let stellar_evolution_type = EvolutionType::MDwarf;
    //let stellar_evolution_type = EvolutionType::SolarLike(SolarEvolutionType::ConstantDissipation);
    //let stellar_evolution_type = EvolutionType::SolarLike(SolarEvolutionType::EvolvingDissipation(star_mass));
    //let star_mass = 1.0;
    //let stellar_evolution_type = EvolutionType::Jupiter;
    let stellar_evolution_type = EvolutionType::NonEvolving;
    let star = Particle::new(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, 
                                            star_love_number, star_fluid_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin, stellar_evolution_type);
    ////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////
    //---- Planet
    let planet_mass: f64 = 3.0e-6; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    let planet_radius_factor: f64 = 1.;
    let planet_radius: f64 = planet_radius_factor * R_EARTH;
    let planet_love_number: f64 = 0.305; // Earth
    //let planet_love_number: f64 = 0.38; // Gas
    let planet_fluid_love_number: f64 = planet_love_number;
    ////// Disipation factor (sigma)
    let planet_dissipation_factor_scale: f64 = 1.;
    //// Hot Gas Giant:
    //let planet_dissipation_factor: f64 = 2.006*3.845764d4;
    //// Terrestrial:
    let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets (no gas)
    let planet_dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * planet_radius.powi(5));
    ////// Radius of gyration
    let planet_radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet
    //let planet_radius_of_gyration_2: f64 = 2.54e-1; // Gas planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = 0.018;                             // semi-major axis (in AU)
    let e: f64 = 0.1;                               // eccentricity
    let i: f64 = 5. * DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0.;                            // argument of pericentre (degrees)
    let n: f64 = 0. * DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * DEG2RAD;                      // mean anomaly (degrees)
    p = (p + n) * DEG2RAD;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = G*(planet_mass+star_mass);
    let (x, y, z, vx, vy, vz) = calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    //////// Cartesian coordinates
    //let gm: f64 = posidonius::constants::G*(planet_mass+star_mass);
    //let x: f64 = 0.1;
    //let y: f64 = 0.;
    //let z: f64 = 0.;
    //let vx: f64 = 0.;
    //let vy: f64 = core::num::Float::sqrt(gm/x); // Stable orbit
    //let vz: f64 = 0.;
    //println!("** {}", core::num::Float::sqrt(gm/x));

    let planet_position = Axes{x:x, y:y, z:z};
    let planet_velocity = Axes{x:vx, y:vy, z:vz};
    let planet_acceleration = Axes{x:0., y:0., z:0.};
    
    //// t: Orbital period
    //// https://en.wikipedia.org/wiki/Orbital_period#Small_body_orbiting_a_central_body
    //let sun_mass      =  1.98892e30;               // kg
    //let gravitational_constant         =  6.6742367e-11;            // m^3.kg^-1.s^-2
    //let astronomical_units        =  1.49598e11;               // m
    //let t = TWO_PI * ((a*astronomical_units).powi(3)/(star_mass*sun_mass*gravitational_constant)).sqrt(); // seconds
    //let t = t/(60.*60.*24.); // days
    //println!("Orbital period in days: {:e} {}", t, t);
    //// Fig. 2
    //// https://arxiv.org/pdf/1506.01084v1.pdf
    //let recommended_timestep = t/1000.;
    //println!("Recommended time step in days for WHFastHelio: {:e} {}", recommended_timestep, recommended_timestep);

    ////// Initialization of planetary spin
    // Planets obliquities in rad
    let planet_obliquity: f64 = 11.459156 * DEG2RAD; // 0.2 rad
    //// Custom period
    let planet_rotation_period: f64 = 24.; // hours
    let planet_angular_frequency = TWO_PI/(planet_rotation_period/24.); // days^-1
    let planet_keplerian_orbital_elements = calculate_keplerian_orbital_elements(G*star_mass*planet_mass, planet_position, planet_velocity);
    let planet_inclination = planet_keplerian_orbital_elements.3;
    //// Pseudo-synchronization period
    ////let planet_keplerian_orbital_elements = calculate_keplerian_orbital_elements(G*star_mass*planet_mass, planet_position, planet_velocity);
    ////let planet_semi_major_axis = planet_keplerian_orbital_elements.0;
    ////let planet_eccentricity = planet_keplerian_orbital_elements.2;
    //let planet_semi_major_axis = a;
    //let planet_eccentricity = e;
    //let planet_pseudo_synchronization_period = calculate_pseudo_synchronization_period(planet_semi_major_axis, planet_eccentricity, star_mass, planet_mass);
    //let planet_angular_frequency = TWO_PI/(planet_pseudo_synchronization_period/24.); // days^-1
    //
    let mut planet_spin = Axes{x:0., y:0., z:0. };
    if planet_inclination == 0. {
        // No inclination, spin can already be calculated:
        planet_spin.x = planet_angular_frequency * planet_obliquity.sin();
        planet_spin.y = 0.;
        planet_spin.z = planet_angular_frequency * planet_obliquity.cos();
    } else {
        // Calculation of orbital angular momentum (without mass and in AU^2/day)
        let horb_x = planet_position.y * planet_velocity.z - planet_position.z * planet_velocity.y;
        let horb_y = planet_position.z * planet_velocity.x - planet_position.x * planet_velocity.z;
        let horb_z = planet_position.x * planet_velocity.y - planet_position.y * planet_velocity.x;
        let horbn = (horb_x.powi(2) + horb_y.powi(2) + horb_z.powi(2)).sqrt();
        // Spin taking into consideration the inclination:
        planet_spin.x = planet_angular_frequency * (horb_x / (horbn * planet_inclination.sin())) * (planet_obliquity+planet_inclination).sin();
        planet_spin.y = planet_angular_frequency * (horb_y / (horbn * planet_inclination.sin())) * (planet_obliquity+planet_inclination).sin();
        planet_spin.z = planet_angular_frequency * (planet_obliquity+planet_inclination).cos();
    }

    //let planetary_evolution_type = EvolutionType::Jupiter;
    let planetary_evolution_type = EvolutionType::NonEvolving;
    let planet = Particle::new(planet_mass, planet_radius, planet_dissipation_factor, planet_dissipation_factor_scale, 
                                            planet_radius_of_gyration_2, planet_love_number, planet_fluid_love_number,
                                            planet_position, planet_velocity, planet_acceleration, planet_spin, planetary_evolution_type);
    ////////////////////////////////////////////////////////////////////////////

    let universe = Universe::new(vec![star, planet], INTEGRATOR);

    //println!("{:?}", universe);
    universe
}

