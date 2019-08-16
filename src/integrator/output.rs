use std::fs::File;
use std::fs::{OpenOptions};
use std::io::{Write, BufWriter};
use std::io::{Read, BufReader};
use super::super::Integrator;
use super::super::particles::Universe;
use super::super::tools::calculate_keplerian_orbital_elements;
use bincode;
use time;
use serde::ser::Serialize;
use serde_json;
use std::path::Path;
use std::fs;
use super::super::constants::{MIN_ORBITAL_PERIOD_TIME_STEP_RATIO};

pub use super::whfast::*;
pub use super::ias15::*;
pub use super::leapfrog::*;


////////////////////////////////////////////////////////////////////////////////
//- Dump and restore functions
////////////////////////////////////////////////////////////////////////////////

pub fn write_recovery_snapshot<I: Serialize>(snapshot_path: &Path, universe_integrator: &I) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.

    if snapshot_path.exists() {
        //// Backup (keep only 1 per hour thanks to filename collision)
        //let new_extension = format!("{0}.bin", time::now().strftime("%Y-%m-%dT%Hh").unwrap());
        // Backup (keep only 1 every 12 hours thanks to filename collision)
        let new_extension = format!("{0}.bin", time::now().strftime("%Y%m%dT%p").unwrap());
        fs::rename(snapshot_path, snapshot_path.with_extension(new_extension)).unwrap();
    }

    // 1.- Serialize the integrator to be able to resume if the simulation is interrupted
    let mut writer = BufWriter::new(File::create(&snapshot_path).unwrap());

    if snapshot_path.extension().unwrap() == "json" {
        let json_encoded = serde_json::to_string_pretty(&universe_integrator).unwrap();
        writer.write(json_encoded.as_bytes()).unwrap();
    } else {
        // Binary
        bincode::serialize_into(&mut writer, &universe_integrator, bincode::Infinite).unwrap(); // bin
    }

}


pub fn n_bytes_per_particle_in_historic_snapshot() -> u64 {
    let n_stored_fields : u64 = 32;
    let n_bytes_per_particle = 8+8+4+8*(n_stored_fields-3);
    n_bytes_per_particle
}

pub fn get_universe_history_writer(universe_history_path: &Path, expected_n_bytes: u64) -> BufWriter<File> {
    // If history snapshot exists:
    // - It validates that it contains all the expected data
    // - If it contains more, the extra data is erased

    // We create file options to write
    let mut options_bin = OpenOptions::new();
    options_bin.create(true).write(true).append(true);

    let universe_history_file = match options_bin.open(&universe_history_path) {
        Ok(f) => f,
        Err(e) => panic!("[PANIC {} UTC] File error: {}", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), e),
    };

    let metadata = universe_history_file.metadata().unwrap();
    let current_n_bytes = metadata.len();
    if current_n_bytes < expected_n_bytes {
        panic!("[PANIC {} UTC] Historic snapshots do not contain all the expected history ({} bytes) as indicated by the recovery snapshot ({} bytes)", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), current_n_bytes, expected_n_bytes);
    }

    // Keep only historic data that saved until the restored snapshot (if there it is the case)
    universe_history_file.set_len(expected_n_bytes).unwrap();
    ////////////////////////////////////////////////////////////////////////////
    let universe_history_writer = BufWriter::new(universe_history_file);
    ////////////////////////////////////////////////////////////////////////////
    universe_history_writer
}

pub fn write_historic_snapshot<T: Write>(universe_history_writer: &mut BufWriter<T>, universe: &Universe, current_time: f64, time_step: f64) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.

    // 2.- Write accumulative output data to conserve the history of the simulation
    let total_energy = universe.compute_total_energy();
    let total_angular_momentum = universe.compute_total_angular_momentum();
    for particle in universe.particles[..universe.n_particles].iter() {
        // Serialize in chunks of maximum 12 elements or it fails
        let output = (
                        current_time,                           // days
                        time_step,                              // days
                        (particle.id as i32),
                        particle.heliocentric_position.x,       // AU
                        particle.heliocentric_position.y,
                        particle.heliocentric_position.z,
                        particle.spin.x,                        // AU/days
                        particle.spin.y,
                        particle.spin.z,
                        particle.heliocentric_velocity.x,       // AU/days
                        particle.heliocentric_velocity.y,
                        particle.heliocentric_velocity.z,
                    );
        bincode::serialize_into(universe_history_writer, &output, bincode::Infinite).unwrap();

        if particle.id != universe.most_massive_particle_index {
            //// Only for planets
            let star_id = universe.most_massive_particle_index;
            let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period) = calculate_keplerian_orbital_elements(universe.particles[star_id].mass_g+particle.mass_g, particle.heliocentric_position, particle.heliocentric_velocity);
            // Control once in a while (when historic point is written) that the
            // time step is small enough to correctly integrate an orbit
            if orbital_period <= time_step*MIN_ORBITAL_PERIOD_TIME_STEP_RATIO {
                println!("\n");
                panic!("[PANIC {} UTC] Time step is too large! Particle {} has an orbital period around particle {} of {:0.3} days which is less than the recommended limit ({:0.3} days) based on the current time step ({:0.3} days).", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), particle.id, star_id, orbital_period, time_step*MIN_ORBITAL_PERIOD_TIME_STEP_RATIO, time_step);
            }
            // Calculation of orbital angular momentum (without mass and in AU^2/day)
            let horb_x = particle.heliocentric_position.y * particle.heliocentric_velocity.z - particle.heliocentric_position.z * particle.heliocentric_velocity.y;
            let horb_y = particle.heliocentric_position.z * particle.heliocentric_velocity.x - particle.heliocentric_position.x * particle.heliocentric_velocity.z;
            let horb_z = particle.heliocentric_position.x * particle.heliocentric_velocity.y - particle.heliocentric_position.y * particle.heliocentric_velocity.x;
            let horbn = (horb_x.powf(2.) + horb_y.powf(2.) + horb_z.powf(2.)).sqrt();
            let output = (
                            semimajor_axis,                     // AU
                            perihelion_distance,                // AU
                            eccentricity,
                            inclination,                        // Degrees
                            longitude_of_perihelion,
                            longitude_of_ascending_node,
                            mean_anomaly,
                            horb_x,
                            horb_y,
                            horb_z,
                            horbn,
                            particle.tides.parameters.internal.denergy_dt,                // Msun.AU^2.day^-3
                            particle.disk.parameters.internal.migration_timescale,
                        );
            bincode::serialize_into(universe_history_writer, &output, bincode::Infinite).unwrap();
        } else {
            let output = (
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                            0.,
                        );
            bincode::serialize_into(universe_history_writer, &output, bincode::Infinite).unwrap();
        }

        let output = (
                        total_energy,
                        total_angular_momentum,
                        particle.mass,                          // Msun
                        particle.radius,                        // Rsun
                        particle.radius_of_gyration_2,
                        particle.tides.parameters.input.scaled_dissipation_factor,
                        particle.tides.parameters.input.love_number,
                        particle.tides.parameters.internal.lag_angle,
                    );
        bincode::serialize_into(universe_history_writer, &output, bincode::Infinite).unwrap();

    }
}

////////////////////////////////////////////////////////////////////////////////
//- Restore functions
////////////////////////////////////////////////////////////////////////////////

pub fn restore_snapshot(universe_integrator_snapshot_path: &Path) -> Result<Box<Integrator>, String> {
    let mut universe_integrator: Box<Integrator>;
    if universe_integrator_snapshot_path.exists() {
        // Open the path in read-only mode only to verify it exists, returns `io::Result<File>`
        if let Err(why) = File::open(&universe_integrator_snapshot_path) {
            return Err(format!("Couldn't open {}: {}", universe_integrator_snapshot_path.display(), why))
        }

        if universe_integrator_snapshot_path.extension().unwrap() == "json" { 
            universe_integrator = deserialize_json_snapshot(&universe_integrator_snapshot_path).unwrap()
        } else {
            universe_integrator = deserialize_bin_snapshot(&universe_integrator_snapshot_path).unwrap();
        }
        if universe_integrator.get_current_time() == 0. {
            println!("[INFO {} UTC] Created new simulation based on '{}'.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display());
            universe_integrator.initialize_physical_values();
        } else {
            println!("[INFO {} UTC] Restored previous simulation from '{}'.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display());
            let current_time_years = universe_integrator.get_current_time()/365.25;
            println!("[INFO {} UTC] Continuing from year {:0.0} ({:0.1e}).", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), current_time_years, current_time_years);
        }
        return Ok(universe_integrator);
    } else {
        return Err(format!("File does not exist"));
    }
}

fn deserialize_json_snapshot(snapshot_path: &Path) -> Result<Box<Integrator>, String> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let mut snapshot_file = File::open(&snapshot_path).unwrap();
    //// Deserialize using `json::decode`
    let mut json_encoded = String::new();
    match snapshot_file.read_to_string(&mut json_encoded) {
        Err(why) => return Err(format!("Couldn't read json snapshot file: {}", why)),
        Ok(_) => {}
    }

    let wrapped_universe_integrator: Result<WHFast, serde_json::Error> = serde_json::from_str(&json_encoded);
    match wrapped_universe_integrator {
        Ok(universe_integrator) => {
            println!("[INFO {} UTC] WHFAST Integrator.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
            Ok(Box::new(universe_integrator))
        },
        Err(_) => {
            let wrapped_universe_integrator: Result<Ias15, serde_json::Error> = serde_json::from_str(&json_encoded);
            match wrapped_universe_integrator {
                Ok(universe_integrator) => {
                    println!("[INFO {} UTC] IAS15 Integrator.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
                    Ok(Box::new(universe_integrator))
                },
                Err(_) => {
                    let wrapped_universe_integrator: Result<LeapFrog, serde_json::Error> = serde_json::from_str(&json_encoded);
                    match wrapped_universe_integrator {
                        Ok(universe_integrator) => {
                            println!("[INFO {} UTC] LeapFrog Integrator.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
                            Ok(Box::new(universe_integrator))
                        },
                        Err(_) => Err(format!("Unknown integrator!")),
                    }
                }
            }
        }
    }
}

fn deserialize_bin_snapshot(snapshot_path: &Path) -> Result<Box<Integrator>, String> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let snapshot_file = File::open(&snapshot_path).unwrap();
    let mut reader = BufReader::new(&snapshot_file);
    let wrapped_universe_integrator: Result<WHFast, bincode::Error> = bincode::deserialize_from(&mut reader, bincode::Infinite);
    match wrapped_universe_integrator {
        Ok(universe_integrator) => {
            println!("[INFO {} UTC] WHFAST Integrator.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
            Ok(Box::new(universe_integrator))
        },
        Err(_) => {
            // Re-open file because the previous File/BufReader was already consumed
            let snapshot_file = File::open(&snapshot_path).unwrap();
            let mut reader = BufReader::new(&snapshot_file);
            let wrapped_universe_integrator: Result<Ias15, bincode::Error> = bincode::deserialize_from(&mut reader, bincode::Infinite);
            match wrapped_universe_integrator {
                Ok(universe_integrator) => {
                    println!("[INFO {} UTC] IAS15 Integrator.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
                    Ok(Box::new(universe_integrator))
                },
                Err(_) => {
                    // Re-open file because the previous File/BufReader was already consumed
                    let snapshot_file = File::open(&snapshot_path).unwrap();
                    let mut reader = BufReader::new(&snapshot_file);
                    let wrapped_universe_integrator: Result<LeapFrog, bincode::Error> = bincode::deserialize_from(&mut reader, bincode::Infinite);
                    match wrapped_universe_integrator {
                        Ok(universe_integrator) => {
                            println!("[INFO {} UTC] LeapFrog Integrator.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
                            Ok(Box::new(universe_integrator))
                        },
                        Err(_) => Err(format!("Unknown integrator!")),
                    }
                }
            }
        }
    }
}

