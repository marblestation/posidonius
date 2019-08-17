extern crate posidonius;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

mod common;
use std::path::Path;

fn whfast_star_first(alternative_coordinates_type: posidonius::whfast::CoordinatesType) -> posidonius::WHFast {
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

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn whfast_star_second(alternative_coordinates_type: posidonius::whfast::CoordinatesType) -> posidonius::WHFast {
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

    let mut particles = common::planets::basic_configuration(&star);
    // WARNING: With WHFast Jacobi does not compute gravity for the star and a particle,
    // the selection of that particle depends on the provided list of particles, hence
    // different orders will lead to slightly different results that will make the test fail
    // But only switching the star from position 0 to 1 makes the algorithm select
    // always the same pair of particles, and the test is valid in this case.
    particles.insert(1, star);
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn leapfrog_star_first() -> posidonius::LeapFrog {
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

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn leapfrog_star_second() -> posidonius::LeapFrog {
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

    let mut particles = common::planets::basic_configuration(&star);
    // WARNING: With WHFast Jacobi does not compute gravity for the star and a particle,
    // the selection of that particle depends on the provided list of particles, hence
    // different orders will lead to slightly different results that will make the test fail
    // But only switching the star from position 0 to 1 makes the algorithm select
    // always the same pair of particles, and the test is valid in this case.
    particles.insert(1, star);
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn ias15_star_first() -> posidonius::Ias15 {
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

    let mut particles = vec![star];
    particles.extend(common::planets::basic_configuration(&star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn ias15_star_second() -> posidonius::Ias15 {
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

    let mut particles = common::planets::basic_configuration(&star);
    // WARNING: With WHFast Jacobi does not compute gravity for the star and a particle,
    // the selection of that particle depends on the provided list of particles, hence
    // different orders will lead to slightly different results that will make the test fail
    // But only switching the star from position 0 to 1 makes the algorithm select
    // always the same pair of particles, and the test is valid in this case.
    particles.insert(1, star);
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn ias15_binary_setup1() -> posidonius::Ias15 {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: true,
        general_relativity: true,
        disk: true,
        wind: true,
        evolution: true,
    };
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let primary_star_mass: f64 = 1.0; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(primary_star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(primary_star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(primary_star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    let primary_star_evolution = posidonius::EvolutionType::BolmontMathis2016(primary_star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(primary_star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let primary_star = common::stars::solar_like_primary(primary_star_mass, primary_star_evolution);
    let secondary_star_mass: f64 = 1.2; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(secondary_star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(secondary_star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(secondary_star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    let secondary_star_evolution = posidonius::EvolutionType::BolmontMathis2016(secondary_star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(secondary_star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let secondary_star = common::stars::solar_like_secondary(&primary_star, secondary_star_mass, secondary_star_evolution, general_relativity_implementation);


    let mut particles = vec![primary_star, secondary_star];
    particles.extend(common::planets::basic_configuration(&primary_star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

fn ias15_binary_setup2() -> posidonius::Ias15 {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) = common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: true,
        general_relativity: true,
        disk: true,
        wind: true,
        evolution: true,
    };
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full

    let primary_star_mass: f64 = 1.0; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(primary_star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(primary_star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(primary_star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    let primary_star_evolution = posidonius::EvolutionType::BolmontMathis2016(primary_star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(primary_star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let primary_star = common::stars::solar_like_primary(primary_star_mass, primary_star_evolution);
    let secondary_star_mass: f64 = 1.2; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(secondary_star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(secondary_star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(secondary_star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    let secondary_star_evolution = posidonius::EvolutionType::BolmontMathis2016(secondary_star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(secondary_star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let secondary_star = common::stars::solar_like_secondary(&primary_star, secondary_star_mass, secondary_star_evolution, general_relativity_implementation);

    let mut particles = vec![secondary_star, primary_star];
    particles.extend(common::planets::basic_configuration(&primary_star));
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    //let universe_integrator = posidonius::WHFast::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe, alternative_coordinates_type);
    universe_integrator
}

#[test]
fn order_whfast_jacobi() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "whfast_jacobi");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    let mut universe_integrator = whfast_star_first(alternative_coordinates_type);
    let mut parallel_universe_integrator = whfast_star_second(alternative_coordinates_type);
    common::universe::iterate(&mut universe_integrator);
    common::universe::iterate(&mut parallel_universe_integrator);
    common::universe::assert_when_star_is_swapped(&universe_integrator.universe, &parallel_universe_integrator.universe);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn order_whfast_democraticheliocentric() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "whfast_democraticheliocentric");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    let mut universe_integrator = whfast_star_first(alternative_coordinates_type);
    let mut parallel_universe_integrator = whfast_star_second(alternative_coordinates_type);
    common::universe::iterate(&mut universe_integrator);
    common::universe::iterate(&mut parallel_universe_integrator);
    common::universe::assert_when_star_is_swapped(&universe_integrator.universe, &parallel_universe_integrator.universe);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn order_whfast_whds() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "whfast_whds");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    let mut universe_integrator = whfast_star_first(alternative_coordinates_type);
    let mut parallel_universe_integrator = whfast_star_second(alternative_coordinates_type);
    common::universe::iterate(&mut universe_integrator);
    common::universe::iterate(&mut parallel_universe_integrator);
    common::universe::assert_when_star_is_swapped(&universe_integrator.universe, &parallel_universe_integrator.universe);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn order_leapfrog() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "leapfrog");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = leapfrog_star_first();
    let mut parallel_universe_integrator = leapfrog_star_second();
    common::universe::iterate(&mut universe_integrator);
    common::universe::iterate(&mut parallel_universe_integrator);
    common::universe::assert_when_star_is_swapped(&universe_integrator.universe, &parallel_universe_integrator.universe);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn order_ias15() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "ias15");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = ias15_star_first();
    let mut parallel_universe_integrator = ias15_star_second();
    common::universe::iterate(&mut universe_integrator);
    common::universe::iterate(&mut parallel_universe_integrator);
    common::universe::assert_when_star_is_swapped(&universe_integrator.universe, &parallel_universe_integrator.universe);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn order_ias15_binary() {
    let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "ias15_binary");
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = ias15_binary_setup1();
    let mut parallel_universe_integrator = ias15_binary_setup2();
    common::universe::iterate(&mut universe_integrator);
    common::universe::iterate(&mut parallel_universe_integrator);
    common::universe::assert_when_star_is_swapped(&universe_integrator.universe, &parallel_universe_integrator.universe);
    common::universe::store_positions_unless_files_exist(&universe_integrator.universe, &rust_data_dirname);
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}
