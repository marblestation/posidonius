extern crate bincode;
use std::io::{Write, BufWriter};
use super::super::particle::Particles;
use super::super::tools::calculate_keplerian_orbital_elements;


pub fn write_bin_snapshot<T: Write>(output_bin: &mut BufWriter<T>, universe: &Particles, current_time: f64, time_step: f64) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.

    let total_energy = universe.compute_total_energy();
    let total_angular_momentum = universe.compute_total_angular_momentum();
    for (i, particle) in universe.particles.iter().enumerate() {
        // Serialize in chunks of maximum 12 elements or it fails
        let output = (
                        current_time,                           // days
                        time_step,                              // days
                        (i as i32),
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
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();
        let output = (
                        particle.acceleration.x,                // AU^2/days
                        particle.acceleration.y,
                        particle.acceleration.z,
                        particle.dspin_dt.x,                    // AU^2/days
                        particle.dspin_dt.y,
                        particle.dspin_dt.z,
                        particle.torque.x,
                        particle.torque.y,
                        particle.torque.z,
                    );
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();
        let output = (
                        particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide,
                        particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide,
                        particle.radial_component_of_the_tidal_force,
                        particle.tidal_acceleration.x,          // AU^2/days
                        particle.tidal_acceleration.y,
                        particle.tidal_acceleration.z,
                        particle.radial_velocity,               // AU/days
                        particle.norm_velocity_vector,
                        particle.distance,                      // AU
                    );
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();

        if i > 0 {
            //// Only for planets
            let star = 0;
            let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly) = calculate_keplerian_orbital_elements(universe.particles[star].mass_g+particle.mass_g, particle.position, particle.velocity);
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
            bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();
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
            bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();
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
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();

    }
}
