extern crate posidonius;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

mod common;
use std::path::Path;

////////////////////////////////////////////////////////////////////////////////
fn enabled_planet_tides_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: true,
        general_relativity: true,
        disk: true,
        wind: true,
        evolution: true,
    };
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.08; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::brown_dwarf(star_mass, star_evolution, general_relativity_implementation);

    let mut particles = vec![star.clone()];
    particles.extend(common::planets::basic_configuration(&star));
    //////////////////////////////////////////////////////////////////////////////////
    // Disable all first or `check_uniform_viscosity_coefficient` will fail
    let disabled_tides = posidonius::Tides::new(posidonius::TidesEffect::Disabled);
    let disabled_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    for particle in particles.iter_mut() {
        particle.set_tides(disabled_tides);
        particle.set_rotational_flattening(disabled_rotational_flattening);
    }
    //////////////////////////////////////////////////////////////////////////////////
    let star_tides = disabled_tides;
    // Load planet data
    let planet_file = "./input/love_numbers/TRAPPIST-1_Earth-like/Results_Trappist1_b_Fe_90_Si_02_670K_freq_Imk2_posidonius.txt";
    let planet_tidal_model_params = common::load_kaula_parameters(planet_file).unwrap();
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::Kaula(planet_tidal_model_params)));
    //
    if let Some((star, planets)) = particles.split_first_mut() {
        // Change from constant time lag to kaula:
        star.set_tides(star_tides);
        for planet in planets.iter_mut() {
            planet.set_tides(planet_tides);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn enabled_planet_tides_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "enabled_planet_tides");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = enabled_planet_tides_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//// TODO:
//#[test]
//fn enabled_tides_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "enabled_tides");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = enabled_planet_tides_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::one_step(&mut universe_integrator);
    //let parallel_universe_from_json = common::universe::one_step_from_json(&python_data_dirname);
    //common::universe::loosly_assert(&universe_integrator.universe, &parallel_universe_from_json); // Built in situ vs Restored from JSON: it requires larger precision and only one step
    //let parallel_universe_from_json = common::universe::iterate_universe_from_json(&python_data_dirname);
    //let universe_from_json = common::universe::iterate_universe_from_json(&rust_data_dirname);
    //common::universe::assert(&universe_from_json, &parallel_universe_from_json);
//}

////////////////////////////////////////////////////////////////////////////////


fn enabled_both_tides_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: true,
        general_relativity: true,
        disk: true,
        wind: true,
        evolution: true,
    };
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.08; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::brown_dwarf(star_mass, star_evolution, general_relativity_implementation);

    let mut particles = vec![star.clone()];
    particles.extend(common::planets::basic_configuration(&star));
    //////////////////////////////////////////////////////////////////////////////////
    // Disable all first or `check_uniform_viscosity_coefficient` will fail
    let disabled_tides = posidonius::Tides::new(posidonius::TidesEffect::Disabled);
    let disabled_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    for particle in particles.iter_mut() {
        particle.set_tides(disabled_tides);
        particle.set_rotational_flattening(disabled_rotational_flattening);
    }
    //////////////////////////////////////////////////////////////////////////////////
    // Load star data
    let star_file = "./input/love_numbers/Aurelie_Stellar/alpha0.51600_P1p2_Ek1p5em6.txt";
    let star_tidal_model_params = common::load_kaula_parameters(star_file).unwrap();
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(posidonius::TidalModel::Kaula(star_tidal_model_params)));
    // Load planet data
    let planet_file = "./input/love_numbers/TRAPPIST-1_Earth-like/Results_Trappist1_b_Fe_90_Si_02_670K_freq_Imk2_posidonius.txt";
    let planet_tidal_model_params = common::load_kaula_parameters(planet_file).unwrap();
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::Kaula(planet_tidal_model_params)));
    //
    if let Some((star, planets)) = particles.split_first_mut() {
        // Change from constant time lag to kaula:
        star.set_tides(star_tides);
        for planet in planets.iter_mut() {
            planet.set_tides(planet_tides);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn enabled_both_tides_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "enabled_both_tides");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = enabled_both_tides_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//// TODO:
//#[test]
//fn enabled_tides_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "enabled_tides");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = enabled_both_tides_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::one_step(&mut universe_integrator);
    //let parallel_universe_from_json = common::universe::one_step_from_json(&python_data_dirname);
    //common::universe::loosly_assert(&universe_integrator.universe, &parallel_universe_from_json); // Built in situ vs Restored from JSON: it requires larger precision and only one step
    //let parallel_universe_from_json = common::universe::iterate_universe_from_json(&python_data_dirname);
    //let universe_from_json = common::universe::iterate_universe_from_json(&rust_data_dirname);
    //common::universe::assert(&universe_from_json, &parallel_universe_from_json);
//}

////////////////////////////////////////////////////////////////////////////////




