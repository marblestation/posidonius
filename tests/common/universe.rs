extern crate assert_approx_eq;
extern crate serde;
extern crate serde_derive;
extern crate serde_json;
use std::fs;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::io::{Read};
use std::path::Path;
use self::assert_approx_eq::assert_approx_eq;
use serde::ser::Serialize;
use posidonius::Integrator;

pub fn store_positions_unless_files_exist(universe: &posidonius::Universe, dirname: &String) {
    let _ = fs::create_dir(&dirname);
    for (i, particle) in universe.particles[..universe.n_particles].iter().enumerate() {
        let expected_position_filename = format!("{0}/position_{1}.json", dirname, i);
        if ! Path::new(&expected_position_filename).exists() {
            let json_encoded = serde_json::to_string(&particle.position).unwrap();
            let mut writer = BufWriter::new(File::create(&expected_position_filename).unwrap());
            writer.write(json_encoded.as_bytes()).unwrap();
        }
    }
}

pub fn store_unless_files_exist<I: Serialize>(universe_integrator: &I, dirname: &String) {
    let _ = fs::create_dir(&dirname);
    let snapshot_filename = format!("{0}/case.json", dirname);
    let snapshot_path = Path::new(&snapshot_filename);
    if ! Path::new(&snapshot_path).exists() {
        posidonius::output::write_recovery_snapshot(snapshot_path, &universe_integrator);
    }
}


pub fn assert_stored_positions(universe: &posidonius::Universe, dirname: &String) {
    let precision = 1.0e-13;
    for (i, particle) in universe.particles[..universe.n_particles].iter().enumerate() {
        let expected_position_filename = format!("{0}/position_{1}.json", dirname, i);
        let expected_position_path = Path::new(&expected_position_filename);
        // Open the path in read-only mode, returns `io::Result<File>`
        let mut expected_position_file = File::open(&expected_position_path).unwrap();
        //// Deserialize using `json::decode`
        let mut expected_position_json = String::new();
        let _ = expected_position_file.read_to_string(&mut expected_position_json).unwrap();
        let expected_position: posidonius::Axes = serde_json::from_str(&expected_position_json).unwrap();
        assert_approx_eq!(particle.position.x, expected_position.x, precision);
        assert_approx_eq!(particle.position.y, expected_position.y, precision);
        assert_approx_eq!(particle.position.z, expected_position.z, precision);
    }
}

pub fn iterate<T>(universe_integrator: &mut T) where T: posidonius::Integrator {
    let universe_history_filename = "/tmp/delete_me.dump";
    let universe_history_path = Path::new(universe_history_filename);
    let expected_n_bytes = 0;
    let silent_mode = true;
    let mut universe_history_writer = posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);
    universe_integrator.initialize_physical_values();
    loop {
        match universe_integrator.iterate(&mut universe_history_writer, silent_mode) {
            Ok(_) => { },
            Err(e) => { break; }
        };
    }
    let _ = fs::remove_file(universe_history_filename);
}

fn iterate_box(universe_integrator: &mut Box<dyn posidonius::Integrator>) {
    let universe_history_filename = "/tmp/delete_me.dump";
    let universe_history_path = Path::new(universe_history_filename);
    let expected_n_bytes = 0;
    let silent_mode = true;
    let mut universe_history_writer = posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);
    universe_integrator.initialize_physical_values();
    loop {
        match universe_integrator.iterate(&mut universe_history_writer, silent_mode) {
            Ok(_) => { break; },
            Err(e) => { break; }
        };
    }
    let _ = fs::remove_file(universe_history_filename);
}

pub fn iterate_universe_from_python_generated_json(dirname: &String) -> posidonius::Universe {
    let snapshot_from_python_filename = format!("{0}/case.json", dirname);
    let snapshot_from_python_path = Path::new(&snapshot_from_python_filename);
   
    let mut boxed_universe_integrator_from_python : Box<dyn posidonius::Integrator> = posidonius::output::restore_snapshot(&snapshot_from_python_path).unwrap();
    iterate_box(&mut boxed_universe_integrator_from_python);
    // Extract universe from the integrator
    let universe_from_python: posidonius::Universe = match boxed_universe_integrator_from_python.as_any().downcast_ref::<posidonius::WHFast>() {
            Some(universe_integrator_from_python) => { universe_integrator_from_python.universe.clone() },
            None => {
                match boxed_universe_integrator_from_python.as_any().downcast_ref::<posidonius::Ias15>() {
                    Some(universe_integrator_from_python) => { universe_integrator_from_python.universe.clone() },
                    None => {
                        match boxed_universe_integrator_from_python.as_any().downcast_ref::<posidonius::LeapFrog>() {
                            Some(universe_integrator_from_python) => { universe_integrator_from_python.universe.clone() },
                            None => {
                                panic!("Unknown integrator!");
                            }
                        }
                    }
                }
            }
    };
    universe_from_python
}

pub fn assert(universe: &posidonius::Universe, parallel_universe: &posidonius::Universe) {
    let precision = 1.0e-13;

    for (particle, parallel_particle) in universe.particles[..universe.n_particles].iter()
                                                                .zip(parallel_universe.particles[..parallel_universe.n_particles].iter()) {
        assert_approx_eq!(particle.position.x, parallel_particle.position.x, precision);
        assert_approx_eq!(particle.position.y, parallel_particle.position.y, precision);
        assert_approx_eq!(particle.position.z, parallel_particle.position.z, precision);
    }
}
