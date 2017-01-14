use super::super::super::particles::{Axes, Particle, Universe};
use super::super::super::{WHFastHelio};
use super::super::super::particles::{EvolutionType};
use super::super::super::constants::{M_EARTH, R_SUN, TWO_PI, R_EARTH, K2, G, DEG2RAD};
use super::super::super::tools::{calculate_cartesian_coordinates, calculate_keplerian_orbital_elements};
use super::super::super::tools::{calculate_spin};

pub fn case7() -> WHFastHelio {
    ////////////////////////////////////////////////////////////////////////////
    // Initial conditions from CASE 7 in Bolmont et al. 2015
    ////////////////////////////////////////////////////////////////////////////
    let time_step: f64 = 0.08; // in days
    let time_limit: f64 = 365.25 * 1.0e8; // days
    //let time_limit: f64 = time_step*4.;
    let initial_time: f64 = 1.0e6*365.25; // time [days] where simulation starts
    let historic_snapshot_period: f64 = 100.*365.25; // days
    let recovery_snapshot_period: f64 = 100.*historic_snapshot_period; // days
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_general_relativy = false;
    let consider_all_body_interactions = false;

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
                                            star_position, star_velocity, star_acceleration, star_spin, 
                                            stellar_evolution_type, initial_time, time_limit);
    ////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////
    //---- Inner Planet
    let inner_planet_mass_factor: f64 = 1.0;
    let inner_planet_mass: f64 = inner_planet_mass_factor * M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    let inner_planet_radius_factor: f64 = 1.;
    let inner_planet_radius: f64 = inner_planet_radius_factor * R_EARTH;
    let inner_planet_love_number: f64 = 0.305; // Earth
    let inner_planet_fluid_love_number: f64 = inner_planet_love_number;
    ////// Disipation factor (sigma)
    let inner_planet_dissipation_factor_scale: f64 = 1.;
    let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets (no gas)
    let inner_planet_dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * inner_planet_radius.powi(5));
    ////// Radius of gyration
    let inner_planet_radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = 0.018;                             // semi-major axis (in AU)
    let e: f64 = 0.01;                               // eccentricity
    let i: f64 = 0. * DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0.;                            // argument of pericentre (degrees)
    let n: f64 = 0. * DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * DEG2RAD;                      // mean anomaly (degrees)
    p = (p + n) * DEG2RAD;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = G*(inner_planet_mass+star_mass);
    let (x, y, z, vx, vy, vz) = calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    let inner_planet_position = Axes{x:x, y:y, z:z};
    let inner_planet_velocity = Axes{x:vx, y:vy, z:vz};
    let inner_planet_acceleration = Axes{x:0., y:0., z:0.};
    
    ////// Initialization of planetary spin
    let inner_planet_obliquity: f64 = 11.459156 * DEG2RAD; // 0.2 rad
    let inner_planet_rotation_period: f64 = 24.; // hours
    //
    let inner_planet_angular_frequency = TWO_PI/(inner_planet_rotation_period/24.); // days^-1
    let inner_planet_keplerian_orbital_elements = calculate_keplerian_orbital_elements(G*star_mass*inner_planet_mass, inner_planet_position, inner_planet_velocity);
    let inner_planet_inclination = inner_planet_keplerian_orbital_elements.3;
    let inner_planet_spin = calculate_spin(inner_planet_angular_frequency, inner_planet_inclination, inner_planet_obliquity, inner_planet_position, inner_planet_velocity);

    let planetary_evolution_type = EvolutionType::NonEvolving;
    let inner_planet = Particle::new(inner_planet_mass, inner_planet_radius, inner_planet_dissipation_factor, inner_planet_dissipation_factor_scale, 
                                            inner_planet_radius_of_gyration_2, inner_planet_love_number, inner_planet_fluid_love_number,
                                            inner_planet_position, inner_planet_velocity, inner_planet_acceleration, inner_planet_spin, 
                                            planetary_evolution_type, initial_time, time_limit);
    ////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////
    //---- Outer Planet
    let outer_planet_mass_factor: f64 = 1.0;
    let outer_planet_mass: f64 = outer_planet_mass_factor * M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    let outer_planet_radius_factor: f64 = 1.;
    let outer_planet_radius: f64 = outer_planet_radius_factor * R_EARTH;
    let outer_planet_love_number: f64 = 0.305; // Earth
    let outer_planet_fluid_love_number: f64 = outer_planet_love_number;
    ////// Disipation factor (sigma)
    let outer_planet_dissipation_factor_scale: f64 = 1.;
    let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets (no gas)
    let outer_planet_dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * outer_planet_radius.powi(5));
    ////// Radius of gyration
    let outer_planet_radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = 0.025;                             // semi-major axis (in AU)
    let e: f64 = 0.01;                               // eccentricity
    let i: f64 = 1. * DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0.;                            // argument of pericentre (degrees)
    let n: f64 = 0. * DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * DEG2RAD;                      // mean anomaly (degrees)
    p = (p + n) * DEG2RAD;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = G*(outer_planet_mass+star_mass);
    let (x, y, z, vx, vy, vz) = calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    let outer_planet_position = Axes{x:x, y:y, z:z};
    let outer_planet_velocity = Axes{x:vx, y:vy, z:vz};
    let outer_planet_acceleration = Axes{x:0., y:0., z:0.};
    
    ////// Initialization of planetary spin
    let outer_planet_obliquity: f64 = 23. * DEG2RAD; // 0.2 rad
    let outer_planet_rotation_period: f64 = 24.; // hours
    //
    let outer_planet_angular_frequency = TWO_PI/(outer_planet_rotation_period/24.); // days^-1
    let outer_planet_keplerian_orbital_elements = calculate_keplerian_orbital_elements(G*star_mass*outer_planet_mass, outer_planet_position, outer_planet_velocity);
    let outer_planet_inclination = outer_planet_keplerian_orbital_elements.3;
    let outer_planet_spin = calculate_spin(outer_planet_angular_frequency, outer_planet_inclination, outer_planet_obliquity, outer_planet_position, outer_planet_velocity);

    let planetary_evolution_type = EvolutionType::NonEvolving;
    let outer_planet = Particle::new(outer_planet_mass, outer_planet_radius, outer_planet_dissipation_factor, outer_planet_dissipation_factor_scale, 
                                            outer_planet_radius_of_gyration_2, outer_planet_love_number, outer_planet_fluid_love_number,
                                            outer_planet_position, outer_planet_velocity, outer_planet_acceleration, outer_planet_spin, 
                                            planetary_evolution_type, initial_time, time_limit);
    ////////////////////////////////////////////////////////////////////////////

    let universe = Universe::new(vec![star, inner_planet, outer_planet], consider_tides, consider_rotational_flattening, consider_general_relativy, consider_all_body_interactions);
    let universe_integrator = WHFastHelio::new(time_step, time_limit, recovery_snapshot_period, historic_snapshot_period, universe);

    universe_integrator
}

