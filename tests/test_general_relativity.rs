extern crate posidonius;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

mod common;
use std::path::Path;

fn kidder1995_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.08; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::brown_dwarf(star_mass, star_evolution_type, general_relativity_implementation);

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles,
                                             consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity);

    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn kidder1995_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "kidder1995");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = kidder1995_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn kidder1995_rust_vs_python() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "kidder1995");
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = kidder1995_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::iterate(&mut universe_integrator);
    let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    common::universe::assert(&universe_integrator.universe, &parallel_universe);
}

////////////////////////////////////////////////////////////////////////////////

fn anderson1975_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.08; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::brown_dwarf(star_mass, star_evolution_type, general_relativity_implementation);

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles,
                                             consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity);

    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn anderson1975_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "anderson1975");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = anderson1975_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn anderson1975_rust_vs_python() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "anderson1975");
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = anderson1975_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::iterate(&mut universe_integrator);
    let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    common::universe::assert(&universe_integrator.universe, &parallel_universe);
}

////////////////////////////////////////////////////////////////////////////////

fn newhall1983_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.08; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::brown_dwarf(star_mass, star_evolution_type, general_relativity_implementation);

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles,
                                             consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity);

    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn newhall1983_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "newhall1983");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = newhall1983_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn newhall1983_rust_vs_python() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "newhall1983");
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = newhall1983_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::iterate(&mut universe_integrator);
    let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    common::universe::assert(&universe_integrator.universe, &parallel_universe);
}

////////////////////////////////////////////////////////////////////////////////

fn none_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = false;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.08; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::brown_dwarf(star_mass, star_evolution_type, general_relativity_implementation);

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles,
                                             consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity);

    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn none_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "none");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = none_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn none_rust_vs_python() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "none");
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = none_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::iterate(&mut universe_integrator);
    let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    common::universe::assert(&universe_integrator.universe, &parallel_universe);
}

////////////////////////////////////////////////////////////////////////////////
