use super::{Particle};

pub fn calculate_norm_spin(particles: &mut [Particle]) {
    for particle in particles.iter_mut() {
        // Squared norm of the spin
        particle.norm_spin_vector_2 = (particle.spin.x.powi(2)) 
                            + (particle.spin.y.powi(2))
                            + (particle.spin.z.powi(2));
    }
}

pub fn calculate_distance_and_velocities(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    // Calculation of velocity vv(j), radial velocity vrad(j)
    // velocities in AU/day
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // Norm of the velocity
        particle.norm_velocity_vector_2 = (particle.velocity.x - tidal_host_particle.velocity.x).powi(2) 
                            + (particle.velocity.y - tidal_host_particle.velocity.y).powi(2)
                            + (particle.velocity.z - tidal_host_particle.velocity.z).powi(2);
        particle.norm_velocity_vector = particle.norm_velocity_vector_2.sqrt();

        // (distance to star)^2
        let distance_2 = (particle.position.x - tidal_host_particle.position.x).powi(2) 
                            + (particle.position.y - tidal_host_particle.position.y).powi(2)
                            + (particle.position.z - tidal_host_particle.position.z).powi(2);
        let distance = distance_2.sqrt();
        particle.distance = distance;

        // Radial velocity
        let v_rad = ((particle.position.x - tidal_host_particle.position.x)*(particle.velocity.x - tidal_host_particle.velocity.x) +
                    (particle.position.y - tidal_host_particle.position.y)*(particle.velocity.y - tidal_host_particle.velocity.y) +
                    (particle.position.z - tidal_host_particle.position.z)*(particle.velocity.z - tidal_host_particle.velocity.z)) / distance;
        particle.radial_velocity = v_rad;
        //let tmp = particle.position.x*particle.velocity.x +
                    //particle.position.y*particle.velocity.y +
                    //particle.position.z*particle.velocity.z;
        //println!("{:e} {:e} {:e}", particle.norm_velocity_vector, particle.distance, particle.radial_velocity);
        //println!("{:e} {:e} {:e}", particle.radial_velocity, tmp, tmp/distance); 
        //println!(">p {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z);
        //println!(">v {:e} {:e} {:e}", particle.velocity.x, particle.velocity.y, particle.velocity.z);
    }
}
