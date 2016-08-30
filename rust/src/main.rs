extern crate posidonius;
extern crate time;
//use std::num;
//use std::io::timer;
//use std::time::Duration;
use std::path::Path;
use std::fs::{OpenOptions};
use std::io::{BufWriter};

use posidonius::Integrator;


fn main() {
    let t1 = time::precise_time_s();
    //timer::sleep(Duration::milliseconds(1000));

    ////////////////////////////////////////////////////////////////////////////
    // Initial conditions from CASE 3 in Bolmont et al. 2015
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    //---- Star (central body)
    let star_mass: f64 = 0.08; // Solar masses
    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    let star_love_number: f64 = 0.307; // M Dwarf
    //let star_love_number: f64 = 0.03;  // Sun
    ////// Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    //// Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    //let star_dissipation_factor: f64 = star_dissipation_factor_scale * 4.992*3.845764e-2; // -66+64
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = star_dissipation_factor_scale * 2.006*3.845764e4; // -60+64
    ////// Radius of gyration
    //let star_radius_of_gyration_2: f64 = 5.9e-2; // Sun
    //let star_radius_of_gyration_2: f64 = 2.0e-1; // M-dwarf
    let star_radius_of_gyration_2: f64 = 1.94e-1; // Brown dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes{x:0., y:0., z:0.};
    let star_velocity = posidonius::Axes{x:0., y:0., z:0.};
    let star_acceleration = posidonius::Axes{x:0., y:0., z:0.};
    ////// Initialization of stellar spin
    let star_rotation_period: f64 = 70.0; // days
    let star_spin0 = posidonius::constants::TWO_PI/star_rotation_period; // days^-1
    let star_spin = posidonius::Axes{x:0., y:0., z:star_spin0 };

    let star = posidonius::Particle::new(star_mass, star_radius, star_dissipation_factor, star_radius_of_gyration_2, star_love_number,
                                            star_position, star_velocity, star_acceleration, star_spin);
    ////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////
    //---- Planet
    let planet_mass: f64 = 3.0e-6; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    //// Mass-radius relationship for earth type planets only:
    //let planet_radius: f64 = posidonius::constants::R_EARTH * ((0.0592*0.7+0.0975) 
                //* core::num::Float::powi((core::num::Float::log10(m(j)) + core::num::Float::log10(posidonius::constants::M2EARTH) - core::num::Float::log10(posidonius::constants::K2)),2)
                //+ (0.2337*0.7+0.4938) 
                //* (core::num::Float::log10(m(j)) + core::num::Float::log10(posidonius::constants::M2EARTH) - core::num::Float::log10(posidonius::constants::K2))+0.3102*0.7+0.7932);
    // For other planet types, specify a factor:
    let planet_radius_factor: f64 = 1.;
    let planet_radius: f64 = planet_radius_factor * posidonius::constants::R_EARTH;
    let planet_love_number: f64 = 0.305; // Earth
    //let planet_love_number: f64 = 0.38; // Gas
    ////// Disipation factor (sigma)
    let planet_dissipation_factor_scale: f64 = 1.;
    //// Hot Gas Giant:
    //let planet_dissipation_factor: f64 = planet_dissipation_factor_scale * 2.006*3.845764d4;
    //// Terrestrial:
    let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets (no gas)
    let planet_dissipation_factor: f64 = planet_dissipation_factor_scale * 2. * posidonius::constants::K2 * k2pdelta/(3. * planet_radius.powf(5.));
    ////// Radius of gyration
    let planet_radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet
    //let planet_radius_of_gyration_2: f64 = 2.54e-1; // Gas planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = 0.018;                               // semi-major axis (in AU)
    let e: f64 = 0.1;                               // eccentricity
    let i: f64 = 5. * posidonius::constants::DEG2RAD;         // inclination (degrees)
    let mut p: f64 = 0.;                           // argument of pericentre (degrees)
    let n: f64 = 0. * posidonius::constants::DEG2RAD;         // longitude of the ascending node (degrees)
    let l: f64 = 0. * posidonius::constants::DEG2RAD;         // mean anomaly (degrees)
    p = (p + n) * posidonius::constants::DEG2RAD;          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = posidonius::constants::G*(planet_mass+star_mass);
    let (x, y, z, vx, vy, vz) = posidonius::general::calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    //////// Cartesian coordinates
    //let gm: f64 = posidonius::constants::G*(planet_mass+star_mass);
    //let x: f64 = 0.1;
    //let y: f64 = 0.;
    //let z: f64 = 0.;
    //let vx: f64 = 0.;
    //let vy: f64 = core::num::Float::sqrt(gm/x); // Stable orbit
    //let vz: f64 = 0.;
    //println!("** {}", core::num::Float::sqrt(gm/x));

    let planet_position = posidonius::Axes{x:x, y:y, z:z};
    let planet_velocity = posidonius::Axes{x:vx, y:vy, z:vz};
    let planet_acceleration = posidonius::Axes{x:0., y:0., z:0.};

    ////// Initialization of planetary spin
    // Planets obliquities in rad
    let planet_obliquity: f64 = 11.459156 * posidonius::constants::DEG2RAD; // 0.2 rad
    //// Custom period
    let planet_rotation_period: f64 = 24.; // hours
    let planet_spin0 = posidonius::constants::TWO_PI/(planet_rotation_period/24.); // days^-1
    let keplerian_orbital_elements = posidonius::general::calculate_keplerian_orbital_elements(posidonius::constants::G*star_mass*planet_mass, planet_position, planet_velocity);
    let inclination = keplerian_orbital_elements.3;
    //// Pseudo-synchronization period
    //let pseudo_sync_period_scale = 1.;
    // TODO: Change ** for core::num::Float::powi(x, 2)
    //let pseudo_sync_period = (1. + 15./2.*e**2 + 45./8.*e**4 + 5./16.*e**6)
                                //* 1. / (1.+3.*e**2 + 3./8.*e**4) * 1./(1-e**2)**1.5*core::num::Float::sqrt(star_mass + planet_mass)
                                //* (q/(1.-e))**(-1.5);
    //let planet_spin0 = pseudo_sync_period_scale * pseudo_sync_period
    //
    let mut planet_spin = posidonius::Axes{x:0., y:0., z:0. };
    if inclination == 0. {
        // No inclination, spin can already be calculated:
        planet_spin.x = planet_spin0 * planet_obliquity.sin();
        planet_spin.y = 0.;
        planet_spin.z = planet_spin0 * planet_obliquity.cos();
    } else {
        // Calculation of orbital angular momentum (without mass and in AU^2/day)
        let horb_x = planet_position.y * planet_velocity.z - planet_position.z * planet_velocity.y;
        let horb_y = planet_position.z * planet_velocity.x - planet_position.x * planet_velocity.z;
        let horb_z = planet_position.x * planet_velocity.y - planet_position.y * planet_velocity.x;
        let horbn = (horb_x.powf(2.) + horb_y.powf(2.) + horb_z.powf(2.)).sqrt();
        // Spin taking into consideration the inclination:
        planet_spin.x = planet_spin0 * (horb_x / (horbn * inclination.sin())) * (planet_obliquity+inclination).sin();
        planet_spin.y = planet_spin0 * (horb_y / (horbn * inclination.sin())) * (planet_obliquity+inclination).sin();
        planet_spin.z = planet_spin0 * (planet_obliquity+inclination).cos();
    }

    let planet = posidonius::Particle::new(planet_mass, planet_radius, planet_dissipation_factor, planet_radius_of_gyration_2, planet_love_number,
                                            planet_position, planet_velocity, planet_acceleration, planet_spin);
    ////////////////////////////////////////////////////////////////////////////



    //let particles_tmp : [posidonius::Particle; posidonius::N_PARTICLES] = [star, planet];
    let particles = posidonius::Particles::new([star, planet]);


    ////////////////////////////////////////////////////////////////////////////
    let path = Path::new("target/output.txt");
    // We create file options to write
    let mut options = OpenOptions::new();
    options.create(true).truncate(true).write(true);

    let output_file = match options.open(&path) {
        Ok(f) => f,
        Err(e) => panic!("file error: {}", e),
    };
    let mut output_writer = BufWriter::new(&output_file);
    ////////////////////////////////////////////////////////////////////////////

    // TODO: Improve (dynamic dispatching?)
    let mut leapfrog = posidonius::LeapFrog::new(posidonius::constants::TIME_STEP, posidonius::constants::TIME_LIMIT, particles);
    let mut ias15 = posidonius::Ias15::new(posidonius::constants::TIME_STEP, posidonius::constants::TIME_LIMIT, particles);

    loop {
        match posidonius::constants::INTEGRATOR {
            posidonius::IntegratorType::LeapFrog => {
                match leapfrog.iterate(&mut output_writer) {
                    Ok(_) => {},
                    Err(e) => { println!("{}", e); break; }
                };
            },
            posidonius::IntegratorType::Ias15 => {
                match ias15.iterate(&mut output_writer) {
                    Ok(_) => {},
                    Err(e) => { println!("{}", e); break; }
                };
            }
        };
    }

    let t2 = time::precise_time_s();
    let d = t2 - t1;
    println!("Execution time: {} seconds", d);
}

