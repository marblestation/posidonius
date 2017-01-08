use super::super::super::particles::{Axes, Particle, Universe};
use super::super::super::particles::{EvolutionType};
use super::super::super::constants::{M_EARTH, R_SUN, TWO_PI, R_EARTH, K2, G, DEG2RAD, INTEGRATOR};
use super::super::super::tools::{calculate_cartesian_coordinates, calculate_keplerian_orbital_elements};
use super::super::super::tools::{calculate_spin};

//pub const TIME_STEP: f64 = 0.08; // in days
//pub const TIME_LIMIT: f64 = 365.25 * 1.0e8;
//pub const PRINT_EVERY_N_DAYS: f64 = 100.*365.25;
//pub const INITIAL_TIME: f64 = 1.0e6*365.25; // time [days] where simulation starts
//
//pub const CONSIDER_EVERY_BODY_COMBINATIONS: bool = false;
//pub const TIDES: bool = true;
//pub const ROTATIONAL_FLATTENING: bool = false;
//pub const GENERAL_RELATIVITY: bool = false;

pub fn case3() -> Universe {
    ////////////////////////////////////////////////////////////////////////////
    // Initial conditions from CASE 3 in Bolmont et al. 2015
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    //---- Star (central body)
    let star_mass: f64 = 0.08; // Solar masses
    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * R_SUN;
    let star_love_number: f64 = 0.307; // Brown Dwarf / M Dwarf
    let star_fluid_love_number: f64 = star_love_number;
    ////// Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 2.006*3.845764e4; // -60+64
    ////// Radius of gyration
    let star_radius_of_gyration_2: f64 = 1.94e-1; // Brown dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = Axes{x:0., y:0., z:0.};
    let star_velocity = Axes{x:0., y:0., z:0.};
    let star_acceleration = Axes{x:0., y:0., z:0.};
    ////// Initialization of stellar spin
    let star_rotation_period: f64 = 70.0; // hours
    let star_angular_frequency = TWO_PI/(star_rotation_period/24.); // days^-1
    let star_spin = Axes{x:0., y:0., z:star_angular_frequency };

    let stellar_evolution_type = EvolutionType::NonEvolving;
    let star = Particle::new(star_mass, star_radius, star_dissipation_factor, star_dissipation_factor_scale, star_radius_of_gyration_2, 
                                            star_love_number, star_fluid_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin, stellar_evolution_type);
    ////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////
    //---- Planet
    let planet_mass_factor: f64 = 1.0;
    let planet_mass: f64 = planet_mass_factor * M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    let planet_radius_factor: f64 = 1.;
    let planet_radius: f64 = planet_radius_factor * R_EARTH;
    let planet_love_number: f64 = 0.305; // Earth
    let planet_fluid_love_number: f64 = planet_love_number;
    ////// Disipation factor (sigma)
    let planet_dissipation_factor_scale: f64 = 1.;
    let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets (no gas)
    let planet_dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * planet_radius.powi(5));
    ////// Radius of gyration
    let planet_radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet

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
    
    let planetary_evolution_type = EvolutionType::NonEvolving;
    let planet = Particle::new(planet_mass, planet_radius, planet_dissipation_factor, planet_dissipation_factor_scale, 
                                            planet_radius_of_gyration_2, planet_love_number, planet_fluid_love_number,
                                            planet_position, planet_velocity, planet_acceleration, planet_spin, planetary_evolution_type);
    ////////////////////////////////////////////////////////////////////////////

    let universe = Universe::new(vec![star, planet], INTEGRATOR);

    //println!("{:?}", universe);
    universe
}

