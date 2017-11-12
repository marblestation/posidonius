extern crate time;
extern crate rustc_serialize;
extern crate bincode;
use std::fs::File;
use std::fs::{OpenOptions};
use std::io::{Write, BufWriter};
use super::super::Integrator;
use super::super::particles::Universe;
use super::super::tools::calculate_keplerian_orbital_elements;
use bincode::rustc_serialize::{encode_into};
use bincode::SizeLimit;
use rustc_serialize::Encodable;
use rustc_serialize::json;
use std::path::Path;
use std::fs;
use super::super::constants::{MIN_ORBITAL_PERIOD_TIME_STEP_RATIO};

pub fn write_recovery_snapshot<I: Integrator+Encodable>(snapshot_path: &Path, universe_integrator: &I) {
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
        let json_encoded = json::encode(&universe_integrator).unwrap();
        writer.write(json_encoded.as_bytes()).unwrap();
    } else {
        // Binary
        encode_into(&universe_integrator, &mut writer, SizeLimit::Infinite).unwrap(); // bin
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
                        particle.position.x,                    // AU
                        particle.position.y,
                        particle.position.z,
                        particle.spin.x,                        // AU/days
                        particle.spin.y,
                        particle.spin.z,
                        particle.velocity.x,                    // AU/days
                        particle.velocity.y,
                        particle.velocity.z,
                    );
        bincode::rustc_serialize::encode_into(&output, universe_history_writer, bincode::SizeLimit::Infinite).unwrap();

        if particle.id > 0 {
            //// Only for planets
            let star_id = 0;
            let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period) = calculate_keplerian_orbital_elements(universe.particles[star_id].mass_g+particle.mass_g, particle.position, particle.velocity);
            // Control once in a while (when historic point is written) that the
            // time step is small enough to correctly integrate an orbit
            if orbital_period <= time_step*MIN_ORBITAL_PERIOD_TIME_STEP_RATIO {
                panic!("[PANIC {} UTC] Time step is too small! Particle {} has an orbital period around particle {} of {:0.3} days which is less than the recommended limit ({:0.3} days) based on the current time step ({:0.3} days).", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), particle.id, star_id, orbital_period, time_step*MIN_ORBITAL_PERIOD_TIME_STEP_RATIO, time_step);
            }
            // Calculation of orbital angular momentum (without mass and in AU^2/day)
            let horb_x = particle.position.y * particle.velocity.z - particle.position.z * particle.velocity.y;
            let horb_y = particle.position.z * particle.velocity.x - particle.position.x * particle.velocity.z;
            let horb_z = particle.position.x * particle.velocity.y - particle.position.y * particle.velocity.x;
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
                            particle.denergy_dt,                // Msun.AU^2.day^-3
                        );
            bincode::rustc_serialize::encode_into(&output, universe_history_writer, bincode::SizeLimit::Infinite).unwrap();
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
                        );
            bincode::rustc_serialize::encode_into(&output, universe_history_writer, bincode::SizeLimit::Infinite).unwrap();
        }

        let output = (
                        total_energy,
                        total_angular_momentum,
                        particle.mass,                          // Msun
                        particle.radius,                        // Rsun
                        particle.radius_of_gyration_2,
                        particle.scaled_dissipation_factor,
                        particle.love_number,
                        particle.lag_angle,
                    );
        bincode::rustc_serialize::encode_into(&output, universe_history_writer, bincode::SizeLimit::Infinite).unwrap();

    }
}

