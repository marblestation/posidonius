use std;
use std::io::{Write, BufWriter};
use super::super::particle::Particles;
use super::super::general::calculate_keplerian_orbital_elements;


pub fn print_output<T: Write>(output_writer: &mut BufWriter<T>, particles: &Particles, current_time: f64, time_step: f64, header: bool) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.
    let debug = true;

    if header {
        // Header
        let _ = output_writer.write(b"current_time\ttime_step\tparticle\t");
        let _ = output_writer.write(b"position_x\tposition_y\tposition_z\t");
        let _ = output_writer.write(b"spin_x\tspin_y\tspin_z\t");
        if debug {
            let _ = output_writer.write(b"velocity_x\tvelocity_y\tvelocity_z\t");
            let _ = output_writer.write(b"acceleration_x\tacceleration_y\tacceleration_z\t");
            let _ = output_writer.write(b"dspin_dt_x\tdspin_dt_y\tdspin_dt_z\t");
            let _ = output_writer.write(b"torque_x\ttorque_y\ttorque_z\t");
            let _ = output_writer.write(b"orthogonal_component_of_the_tidal_force_due_to_stellar_tide\torthogonal_component_of_the_tidal_force_due_to_planetary_tide\t");
            let _ = output_writer.write(b"radial_component_of_the_tidal_force\tradial_component_of_the_tidal_force_conservative_part\tradial_component_of_the_tidal_force_dissipative_part\t");
            let _ = output_writer.write(b"tidal_acceleration_x\ttidal_acceleration_y\ttidal_acceleration_z\t");
            let _ = output_writer.write(b"radial_velocity\tnorm_velocity_vector\tdistance\t");
        }
        let _ = output_writer.write(b"semi-major_axis\tperihelion_distance\teccentricity\tinclination\tlongitude_of_perihelion\tlongitude_of_ascending_node\tmean_anomaly\t");
        let _ = output_writer.write(b"orbital_angular_momentum_x\torbital_angular_momentum_y\torbital_angular_momentum_z\torbital_angular_momentum\t");
        let _ = output_writer.write(b"denergy_dt\t");
        let _ = output_writer.write(b"mass\tradius\t");
        if debug {
            let _ = output_writer.write(b"radius_of_gyration_2\t");
            let _ = output_writer.write(b"dissipation_factor\tlove_number");
        } else {
            let _ = output_writer.write(b"radius_of_gyration_2");
        }

        let _ = output_writer.write(b"\n");
    }

    let current_time_years = current_time/365.25;

    print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
    let _ = std::io::stdout().flush();

    //// SBC
    //for i in 1..N_PARTICLES  {
    //for particle in particles.particles[1..].iter_mut() {
    for (i, particle) in particles.particles.iter().enumerate() {
        let _ = output_writer.write(format!("{}\t{}\t{}\t", current_time_years, time_step, i).as_bytes());
        let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.position.x, particle.position.y, particle.position.z).as_bytes());
        let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.spin.x, particle.spin.y, particle.spin.z).as_bytes());
        if debug {
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.velocity.x, particle.velocity.y, particle.velocity.z).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.acceleration.x, particle.acceleration.y, particle.acceleration.z).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.dspin_dt.x, particle.dspin_dt.y, particle.dspin_dt.z).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.torque.x, particle.torque.y, particle.torque.z).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t", particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.radial_component_of_the_tidal_force, particle.radial_component_of_the_tidal_force_conservative_part, particle.radial_component_of_the_tidal_force_dissipative_part).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.tidal_acceleration.x, particle.tidal_acceleration.y, particle.tidal_acceleration.z).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t", particle.radial_velocity, particle.norm_velocity_vector, particle.distance).as_bytes());
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
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t", semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}\t{:e}\t{:e}\t", horb_x/horbn, horb_y/horbn, horb_z/horbn, horbn).as_bytes());
            let _ = output_writer.write(format!("{:e}\t", particle.denergy_dt).as_bytes());
        } else {
            let _ = output_writer.write(b"\t\t\t\t\t\t\t");
            let _ = output_writer.write(b"\t\t\t\t");
            let _ = output_writer.write(b"\t");
        }
        let _ = output_writer.write(format!("{:e}\t{:e}\t", particle.mass, particle.radius).as_bytes());
        if debug {
            let _ = output_writer.write(format!("{:e}\t", particle.radius_of_gyration_2).as_bytes());
            let _ = output_writer.write(format!("{:e}\t{:e}", particle.dissipation_factor, particle.love_number).as_bytes());
        } else {
            let _ = output_writer.write(format!("{:e}", particle.radius_of_gyration_2).as_bytes());
        }
        let _ = output_writer.write(b"\n");
    }
} 
