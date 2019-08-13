use super::super::constants::{R_SUN};
use super::{Particle};

pub fn calculate_wind_factor(particles: &mut [Particle]) {
    for particle in particles.iter_mut() {
        if particle.wind_k_factor != 0. {
            let threshold = (particle.norm_spin_vector_2).sqrt();
            let factor = - 1. / (particle.moment_of_inertia);
            if threshold >= particle.wind_rotation_saturation {
                particle.wind_factor = factor * particle.wind_k_factor * particle.wind_rotation_saturation_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt()
            } else {
                particle.wind_factor = factor * particle.wind_k_factor * (particle.radius/R_SUN * 1./particle.mass).sqrt()
            }
        };
    }
}
