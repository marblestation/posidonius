use super::{Particle};

pub fn calculate_norm_spin(particles: &mut [Particle]) {
    for particle in particles.iter_mut() {
        // Squared norm of the spin
        particle.norm_spin_vector_2 = (particle.spin.x.powi(2)) 
                            + (particle.spin.y.powi(2))
                            + (particle.spin.z.powi(2));
    }
}
