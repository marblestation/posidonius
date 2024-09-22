extern crate posidonius;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

mod common;

#[macro_use]
extern crate criterion;

use criterion::Criterion;

fn create_universe() -> posidonius::Universe {
    let (_time_step, time_limit, initial_time, _historic_snapshot_period, _recovery_snapshot_period) = common::simulation_properties();
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

    let star_mass: f64 = 1.0; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like_with_disk(star_mass, star_evolution, general_relativity_implementation);

    let mut particles = vec![star.clone()];
    particles.extend(common::planets::basic_configuration(&star));
    let mut universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);
    let current_time = 0.;
    // initialize_physical_values
    let evolution = true;
    universe.calculate_spin_and_evolving_quantities(current_time, evolution); // Make sure we start with the good initial values
    universe.calculate_roche_radiuses(); // Needed for collision detection
    //
    universe.inertial_to_heliocentric();

    universe
}

fn create_kaula_universe() -> posidonius::Universe {
    let (_time_step, time_limit, initial_time, _historic_snapshot_period, _recovery_snapshot_period) = common::simulation_properties();
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

    let star_mass: f64 = 1.0; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    //let star_evolution = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like_with_disk(star_mass, star_evolution, general_relativity_implementation);

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
    let mut universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);
    let current_time = 0.;
    // initialize_physical_values
    let evolution = true;
    universe.calculate_spin_and_evolving_quantities(current_time, evolution); // Make sure we start with the good initial values
    universe.calculate_roche_radiuses(); // Needed for collision detection
    //
    universe.inertial_to_heliocentric();

    universe
}


fn criterion_benchmark_universe(c: &mut Criterion) {
    let mut universe = create_universe();
    let mut kaula_universe = create_kaula_universe();
    let ignore_gravity_terms = posidonius::IgnoreGravityTerms::None;
    let ignored_gravity_terms = ignore_gravity_terms;
    let evolution = true;
    let dangular_momentum_dt_per_moment_of_inertia = true;
    let accelerations = true;
    let current_time = 0.0;

    c.bench_function("gravity_calculate_acceleration", |b| b.iter(|| universe.gravity_calculate_acceleration(ignore_gravity_terms)));

    let mut group = c.benchmark_group("calculate_additional_effects");
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: true,
        general_relativity: true,
        disk: true,
        wind: true,
        evolution: true,
    };
    group.bench_function("all", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    group.bench_function("constant_time_lag_tides", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    kaula_universe.consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    group.bench_function("kaula_tides", |b| b.iter(|| kaula_universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: false,
        rotational_flattening: true,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    group.bench_function("rotational_flattening", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: false,
        rotational_flattening: false,
        general_relativity: true,
        disk: false,
        wind: false,
        evolution: false,
    };
    group.bench_function("general_relativity", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: false,
        rotational_flattening: false,
        general_relativity: false,
        disk: true,
        wind: false,
        evolution: false,
    };
    group.bench_function("disk", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: false,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: true,
        evolution: false,
    };
    group.bench_function("wind", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: false,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: true,
    };
    group.bench_function("evolution", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));
    
    universe.consider_effects = posidonius::ConsiderEffects {
        tides: false,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    group.bench_function("none", |b| b.iter(|| universe.calculate_additional_effects(current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms)));

    group.finish();
}

criterion_group!(benches, criterion_benchmark_universe);
criterion_main!(benches);
