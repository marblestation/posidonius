extern crate time;
extern crate rustc_serialize;
extern crate bincode;
use std::fs::File;
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
            let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly) = calculate_keplerian_orbital_elements(universe.particles[star_id].mass_g+particle.mass_g, particle.position, particle.velocity);
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
                        particle.love_number
                    );
        bincode::rustc_serialize::encode_into(&output, universe_history_writer, bincode::SizeLimit::Infinite).unwrap();

    }
}

