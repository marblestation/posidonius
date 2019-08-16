extern crate posidonius;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

mod common;
use std::path::Path;

fn solar_like_baraffe1998_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 1.0; // Solar masses
    let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like(star_mass, star_evolution_type, general_relativity_implementation);

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
fn solar_like_baraffe1998_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_baraffe1998");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = solar_like_baraffe1998_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn solar_like_baraffe1998_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_baraffe1998");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = solar_like_baraffe1998_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

fn solar_like_baraffe2015_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 1.0; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like(star_mass, star_evolution_type, general_relativity_implementation);

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
fn solar_like_baraffe2015_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_baraffe2015");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = solar_like_baraffe2015_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn solar_like_baraffe2015_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_baraffe2015");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = solar_like_baraffe2015_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////


fn solar_like_bolmontmathis2016_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 1.0; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like(star_mass, star_evolution_type, general_relativity_implementation);

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
fn solar_like_bolmontmathis2016_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_bolmontmathis2016");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = solar_like_bolmontmathis2016_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn solar_like_bolmontmathis2016_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_bolmontmathis2016");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = solar_like_bolmontmathis2016_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////


fn solar_like_galletbolmont2017_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 1.0; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like(star_mass, star_evolution_type, general_relativity_implementation);

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
fn solar_like_galletbolmont2017_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_galletbolmont2017");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = solar_like_galletbolmont2017_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn solar_like_galletbolmont2017_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_galletbolmont2017");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = solar_like_galletbolmont2017_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

fn solar_like_non_evolving_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 1.0; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like(star_mass, star_evolution_type, general_relativity_implementation);

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
fn solar_like_non_evolving_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_non_evolving");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = solar_like_non_evolving_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn solar_like_non_evolving_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "solar_like_non_evolving");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = solar_like_non_evolving_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

fn brown_dwarf_leconte2011_case() -> posidonius::WHFast {
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
fn brown_dwarf_leconte2011_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "brown_dwarf_leconte2011");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = brown_dwarf_leconte2011_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn brown_dwarf_leconte2011_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "brown_dwarf_leconte2011");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = brown_dwarf_leconte2011_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////
    
fn brown_dwarf_non_evolving_case() -> posidonius::WHFast {
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
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    let star_evolution_type = posidonius::EvolutionType::NonEvolving;
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
fn brown_dwarf_non_evolving_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "brown_dwarf_non_evolving");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = brown_dwarf_non_evolving_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn brown_dwarf_non_evolving_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "brown_dwarf_non_evolving");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = brown_dwarf_non_evolving_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

fn m_dwarf_baraffe1998_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.10; // Solar masses
    let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::m_dwarf(star_mass, star_evolution_type, general_relativity_implementation);
    
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
fn m_dwarf_baraffe1998_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "m_dwarf_baraffe1998");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = m_dwarf_baraffe1998_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn m_dwarf_baraffe1998_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "m_dwarf_baraffe1998");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = m_dwarf_baraffe1998_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

fn m_dwarf_baraffe2015_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_tides = true;
    let consider_rotational_flattening = true;
    let consider_disk_interaction = false;
    let consider_general_relativity = true;
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let star_mass: f64 = 0.11; // Solar masses
    //let star_evolution_type = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    let star_evolution_type = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::m_dwarf(star_mass, star_evolution_type, general_relativity_implementation);
    
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
fn m_dwarf_baraffe2015_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "m_dwarf_baraffe2015");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = m_dwarf_baraffe2015_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn m_dwarf_baraffe2015_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "m_dwarf_baraffe2015");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = m_dwarf_baraffe2015_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

fn m_dwarf_non_evolving_case() -> posidonius::WHFast {
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
    //let star_evolution_type = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution_type = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution_type = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    let star_evolution_type = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::m_dwarf(star_mass, star_evolution_type, general_relativity_implementation);
    
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
fn m_dwarf_non_evolving_rust() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "m_dwarf_non_evolving");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = m_dwarf_non_evolving_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

//#[test]
//fn m_dwarf_non_evolving_rust_vs_python() {
    //let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "m_dwarf_non_evolving");
    //let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    //let mut universe_integrator = m_dwarf_non_evolving_case();
    //common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    //common::universe::iterate(&mut universe_integrator);
    //let parallel_universe = common::universe::iterate_universe_from_python_generated_json(&python_data_dirname);
    //common::universe::assert(&universe_integrator.universe, &parallel_universe);
//}

////////////////////////////////////////////////////////////////////////////////

