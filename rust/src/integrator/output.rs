extern crate bincode;
extern crate rusqlite;
use std;
use std::io::{Write, BufWriter};
use super::super::particle::Particles;
use super::super::general::calculate_keplerian_orbital_elements;

pub fn write_db_snapshot(output_db: &rusqlite::Connection, particles: &Particles, current_time: f64, time_step: f64, header: bool) {
    if header {
        output_db.execute("DROP TABLE IF EXISTS particles;", &[]).unwrap();
        output_db.execute("CREATE TABLE particles (
                      id                                                                INTEGER PRIMARY KEY,
                      current_time                                                      REAL,
                      time_step                                                         REAL,
                      particle                                                          INTEGER,
                      position_x                                                        REAL,
                      position_y                                                        REAL,
                      position_z                                                        REAL,
                      spin_x                                                            REAL,
                      spin_y                                                            REAL,
                      spin_z                                                            REAL,
                      velocity_x                                                        REAL,
                      velocity_y                                                        REAL,
                      velocity_z                                                        REAL,
                      acceleration_x                                                    REAL,
                      acceleration_y                                                    REAL,
                      acceleration_z                                                    REAL,
                      dspin_dt_x                                                        REAL,
                      dspin_dt_y                                                        REAL,
                      dspin_dt_z                                                        REAL,
                      torque_x                                                          REAL,
                      torque_y                                                          REAL,
                      torque_z                                                          REAL,
                      orthogonal_component_of_the_tidal_force_due_to_stellar_tide       REAL,
                      orthogonal_component_of_the_tidal_force_due_to_planetary_tide     REAL,
                      radial_component_of_the_tidal_force                               REAL,
                      radial_component_of_the_tidal_force_conservative_part             REAL,
                      radial_component_of_the_tidal_force_dissipative_part              REAL,
                      tidal_acceleration_x                                              REAL,
                      tidal_acceleration_y                                              REAL,
                      tidal_acceleration_z                                              REAL,
                      radial_velocity                                                   REAL,
                      norm_velocity_vector                                              REAL,
                      distance                                                          REAL,
                      semimajor_axis                                                    REAL,
                      perihelion_distance                                               REAL,
                      eccentricity                                                      REAL,
                      inclination                                                       REAL,
                      longitude_of_perihelion                                           REAL,
                      longitude_of_ascending_node                                       REAL,
                      mean_anomaly                                                      REAL,
                      orbital_angular_momentum_x                                        REAL,
                      orbital_angular_momentum_y                                        REAL,
                      orbital_angular_momentum_z                                        REAL,
                      orbital_angular_momentum                                          REAL,
                      denergy_dt                                                        REAL,
                      mass                                                              REAL,
                      radius                                                            REAL,
                      radius_of_gyration_2                                              REAL,
                      dissipation_factor                                                REAL,
                      love_number                                                       REAL
                      )", &[]).unwrap();
    }
    let current_time_years = current_time/365.25;
    for (i, particle) in particles.particles.iter().enumerate() {
        // TODO: Improve the destructure of the calculate_keplerian_orbital_elements results
        let mut semimajor_axis: f64 = 0.;
        let mut perihelion_distance: f64 = 0.;
        let mut eccentricity: f64 = 0.;
        let mut inclination: f64 = 0.;
        let mut longitude_of_perihelion: f64 = 0.;
        let mut longitude_of_ascending_node: f64 = 0.;
        let mut mean_anomaly: f64 = 0.;
        let mut horb_x: f64 = 0.;
        let mut horb_y: f64 = 0.;
        let mut horb_z: f64 = 0.;
        let mut horbn: f64 = 0.;
        if i > 0 {
            //// Not for the central body (star)
            let star = 0;
            let (semimajor_axis0, perihelion_distance0, eccentricity0, inclination0, longitude_of_perihelion0, longitude_of_ascending_node0, mean_anomaly0) = calculate_keplerian_orbital_elements(particles.particles[star].mass_g+particle.mass_g, particle.position, particle.velocity);
            semimajor_axis = semimajor_axis0;
            perihelion_distance = perihelion_distance0;
            eccentricity = eccentricity0;
            inclination = inclination0;
            longitude_of_perihelion = longitude_of_perihelion0;
            longitude_of_ascending_node = longitude_of_ascending_node0;
            mean_anomaly = mean_anomaly0;
            // Calculation of orbital angular momentum (without mass and in AU^2/day)
            horb_x = particle.position.y * particle.velocity.z - particle.position.z * particle.velocity.y;
            horb_y = particle.position.z * particle.velocity.x - particle.position.x * particle.velocity.z;
            horb_z = particle.position.x * particle.velocity.y - particle.position.y * particle.velocity.x;
            horbn = (horb_x.powf(2.) + horb_y.powf(2.) + horb_z.powf(2.)).sqrt();
        }
        output_db.execute("INSERT INTO particles (
                              current_time,
                              time_step,
                              particle,
                              position_x,
                              position_y,
                              position_z,
                              spin_x,
                              spin_y,
                              spin_z,
                              velocity_x,
                              velocity_y,
                              velocity_z,
                              acceleration_x,
                              acceleration_y,
                              acceleration_z,
                              dspin_dt_x,
                              dspin_dt_y,
                              dspin_dt_z,
                              torque_x,
                              torque_y,
                              torque_z,
                              orthogonal_component_of_the_tidal_force_due_to_stellar_tide,
                              orthogonal_component_of_the_tidal_force_due_to_planetary_tide,
                              radial_component_of_the_tidal_force,
                              radial_component_of_the_tidal_force_conservative_part,
                              radial_component_of_the_tidal_force_dissipative_part,
                              tidal_acceleration_x,
                              tidal_acceleration_y,
                              tidal_acceleration_z,
                              radial_velocity,
                              norm_velocity_vector,
                              distance,
                              semimajor_axis,
                              perihelion_distance,
                              eccentricity,
                              inclination,
                              longitude_of_perihelion,
                              longitude_of_ascending_node,
                              mean_anomaly,
                              orbital_angular_momentum_x,
                              orbital_angular_momentum_y,
                              orbital_angular_momentum_z,
                              orbital_angular_momentum,
                              denergy_dt,
                              mass,
                              radius,
                              radius_of_gyration_2,
                              dissipation_factor,
                              love_number                
                      )
                  VALUES (
                            $1, $2, $3, $4, $5, $6, $7, $8, $9, $10,
                            $11, $12, $13, $14, $15, $16, $17, $18, $19, $20,
                            $21, $22, $23, $24, $25, $26, $27, $28, $29, $30,
                            $31, $32, $33, $34, $35, $36, $37, $38, $39, $40,
                            $41, $42, $43, $44, $45, $46, $47, $48, $49
                          )",
                 &[
                    &current_time_years, 
                    &time_step, 
                    &(i as i32),
                    &particle.position.x,
                    &particle.position.y,
                    &particle.position.z,
                    &particle.spin.x,
                    &particle.spin.y,
                    &particle.spin.z,
                    &particle.velocity.x,
                    &particle.velocity.y,
                    &particle.velocity.z,
                    &particle.acceleration.x,
                    &particle.acceleration.y,
                    &particle.acceleration.z,
                    &particle.dspin_dt.x,
                    &particle.dspin_dt.y,
                    &particle.dspin_dt.z,
                    &particle.torque.x,
                    &particle.torque.y,
                    &particle.torque.z,
                    &particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide,
                    &particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide,
                    &particle.radial_component_of_the_tidal_force,
                    &particle.radial_component_of_the_tidal_force_conservative_part,
                    &particle.radial_component_of_the_tidal_force_dissipative_part,
                    &particle.tidal_acceleration.x,
                    &particle.tidal_acceleration.y,
                    &particle.tidal_acceleration.z,
                    &particle.radial_velocity,
                    &particle.norm_velocity_vector,
                    &particle.distance,
                    &semimajor_axis,
                    &perihelion_distance,
                    &eccentricity,
                    &inclination,
                    &longitude_of_perihelion,
                    &longitude_of_ascending_node,
                    &mean_anomaly,
                    &(horb_x/horbn),
                    &(horb_y/horbn),
                    &(horb_z/horbn),
                    &horbn,
                    &particle.denergy_dt,
                    &particle.mass,
                    &particle.radius,
                    &particle.radius_of_gyration_2,
                    &particle.dissipation_factor,
                    &particle.love_number
                 ]).unwrap();
    }
}

pub fn write_txt_snapshot<T: Write>(output_txt: &mut BufWriter<T>, particles: &Particles, current_time: f64, time_step: f64, header: bool) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.
    let debug = true;

    if header {
        // Header
        let _ = output_txt.write(b"current_time\ttime_step\tparticle\t");
        let _ = output_txt.write(b"position_x\tposition_y\tposition_z\t");
        let _ = output_txt.write(b"spin_x\tspin_y\tspin_z\t");
        if debug {
            let _ = output_txt.write(b"velocity_x\tvelocity_y\tvelocity_z\t");
            let _ = output_txt.write(b"acceleration_x\tacceleration_y\tacceleration_z\t");
            let _ = output_txt.write(b"dspin_dt_x\tdspin_dt_y\tdspin_dt_z\t");
            let _ = output_txt.write(b"torque_x\ttorque_y\ttorque_z\t");
            let _ = output_txt.write(b"orthogonal_component_of_the_tidal_force_due_to_stellar_tide\torthogonal_component_of_the_tidal_force_due_to_planetary_tide\t");
            let _ = output_txt.write(b"radial_component_of_the_tidal_force\tradial_component_of_the_tidal_force_conservative_part\tradial_component_of_the_tidal_force_dissipative_part\t");
            let _ = output_txt.write(b"tidal_acceleration_x\ttidal_acceleration_y\ttidal_acceleration_z\t");
            let _ = output_txt.write(b"radial_velocity\tnorm_velocity_vector\tdistance\t");
        }
        let _ = output_txt.write(b"semi-major_axis\tperihelion_distance\teccentricity\tinclination\tlongitude_of_perihelion\tlongitude_of_ascending_node\tmean_anomaly\t");
        let _ = output_txt.write(b"orbital_angular_momentum_x\torbital_angular_momentum_y\torbital_angular_momentum_z\torbital_angular_momentum\t");
        let _ = output_txt.write(b"denergy_dt\t");
        let _ = output_txt.write(b"mass\tradius\t");
        if debug {
            let _ = output_txt.write(b"radius_of_gyration_2\t");
            let _ = output_txt.write(b"dissipation_factor\tlove_number");
        } else {
            let _ = output_txt.write(b"radius_of_gyration_2");
        }

        let _ = output_txt.write(b"\n");
    }

    let current_time_years = current_time/365.25;
    print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
    let _ = std::io::stdout().flush();

    //// SBC
    //for i in 1..N_PARTICLES  {
    //for particle in particles.particles[1..].iter_mut() {
    for (i, particle) in particles.particles.iter().enumerate() {
        let _ = output_txt.write(format!("{}\t{}\t{}\t", current_time_years, time_step, i).as_bytes());
        let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.position.x, particle.position.y, particle.position.z).as_bytes());
        let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.spin.x, particle.spin.y, particle.spin.z).as_bytes());
        if debug {
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.velocity.x, particle.velocity.y, particle.velocity.z).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.acceleration.x, particle.acceleration.y, particle.acceleration.z).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.dspin_dt.x, particle.dspin_dt.y, particle.dspin_dt.z).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.torque.x, particle.torque.y, particle.torque.z).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t", particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.radial_component_of_the_tidal_force, particle.radial_component_of_the_tidal_force_conservative_part, particle.radial_component_of_the_tidal_force_dissipative_part).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.tidal_acceleration.x, particle.tidal_acceleration.y, particle.tidal_acceleration.z).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t", particle.radial_velocity, particle.norm_velocity_vector, particle.distance).as_bytes());
        }
        if i > 0 {
            //// Not for the central body (star)
            let star = 0;
            let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly) = calculate_keplerian_orbital_elements(particles.particles[star].mass_g+particle.mass_g, particle.position, particle.velocity);
            // Calculation of orbital angular momentum (without mass and in AU^2/day)
            let horb_x = particle.position.y * particle.velocity.z - particle.position.z * particle.velocity.y;
            let horb_y = particle.position.z * particle.velocity.x - particle.position.x * particle.velocity.z;
            let horb_z = particle.position.x * particle.velocity.y - particle.position.y * particle.velocity.x;
            let horbn = (horb_x.powf(2.) + horb_y.powf(2.) + horb_z.powf(2.)).sqrt();
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t", semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}\t{:e}\t{:e}\t", horb_x/horbn, horb_y/horbn, horb_z/horbn, horbn).as_bytes());
            let _ = output_txt.write(format!("{:e}\t", particle.denergy_dt).as_bytes());
        } else {
            let _ = output_txt.write(b"\t\t\t\t\t\t\t");
            let _ = output_txt.write(b"\t\t\t\t");
            let _ = output_txt.write(b"\t");
        }
        let _ = output_txt.write(format!("{:e}\t{:e}\t", particle.mass, particle.radius).as_bytes());
        if debug {
            let _ = output_txt.write(format!("{:e}\t", particle.radius_of_gyration_2).as_bytes());
            let _ = output_txt.write(format!("{:e}\t{:e}", particle.dissipation_factor, particle.love_number).as_bytes());
        } else {
            let _ = output_txt.write(format!("{:e}", particle.radius_of_gyration_2).as_bytes());
        }
        let _ = output_txt.write(b"\n");
    }
} 

pub fn write_bin_snapshot<T: Write>(output_bin: &mut BufWriter<T>, particles: &Particles, current_time: f64, time_step: f64) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.

    let current_time_years = current_time/365.25;
    for (i, particle) in particles.particles.iter().enumerate() {
        // Serialize in chunks of maximum 12 elements or it fails
        let output = (
                        current_time_years, 
                        time_step, 
                        (i as i32),
                        particle.position.x,
                        particle.position.y,
                        particle.position.z,
                        particle.spin.x,
                        particle.spin.y,
                        particle.spin.z,
                        particle.velocity.x,
                        particle.velocity.y,
                        particle.velocity.z,
                    );
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();
        let output = (
                        particle.acceleration.x,
                        particle.acceleration.y,
                        particle.acceleration.z,
                        particle.dspin_dt.x,
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
                        particle.radial_component_of_the_tidal_force_conservative_part,
                        particle.radial_component_of_the_tidal_force_dissipative_part,
                        particle.tidal_acceleration.x,
                        particle.tidal_acceleration.y,
                        particle.tidal_acceleration.z,
                        particle.radial_velocity,
                        particle.norm_velocity_vector,
                        particle.distance,
                    );
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();

        let output = (
                        particle.denergy_dt,
                        particle.mass,
                        particle.radius,
                        particle.radius_of_gyration_2,
                        particle.dissipation_factor,
                        particle.love_number
                    );
        bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();

        if i > 0 {
            //// Only for planets
            let star = 0;
            let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly) = calculate_keplerian_orbital_elements(particles.particles[star].mass_g+particle.mass_g, particle.position, particle.velocity);
            // Calculation of orbital angular momentum (without mass and in AU^2/day)
            let horb_x = particle.position.y * particle.velocity.z - particle.position.z * particle.velocity.y;
            let horb_y = particle.position.z * particle.velocity.x - particle.position.x * particle.velocity.z;
            let horb_z = particle.position.x * particle.velocity.y - particle.position.y * particle.velocity.x;
            let horbn = (horb_x.powf(2.) + horb_y.powf(2.) + horb_z.powf(2.)).sqrt();
            let output = (
                            semimajor_axis,
                            perihelion_distance,
                            eccentricity,
                            inclination,
                            longitude_of_perihelion,
                            longitude_of_ascending_node,
                            mean_anomaly,
                            horb_x/horbn,
                            horb_y/horbn,
                            horb_z/horbn,
                            horbn,
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
                        );
            bincode::rustc_serialize::encode_into(&output, output_bin, bincode::SizeLimit::Infinite).unwrap();
        }

    }
}
