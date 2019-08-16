use super::super::constants::{R_SUN};
use super::{Particle};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct WindParticleInputParameters {
    pub k_factor: f64,
    pub rotation_saturation: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct WindParticleInternalParameters {
    pub rotation_saturation_2: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct WindParticleOutputParameters {
    pub factor: f64, // Spin related
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct WindParticleParameters {
    pub input: WindParticleInputParameters,
    pub internal: WindParticleInternalParameters,
    pub output: WindParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum WindEffect {
    Interaction,
    None,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Wind {
    pub effect: WindEffect,
    pub parameters: WindParticleParameters,
}

impl Wind {
    pub fn new(effect: WindEffect, k_factor: f64, rotation_saturation: f64) -> Wind {
        Wind {
            effect: effect,
            parameters: WindParticleParameters {
                input: WindParticleInputParameters {
                    k_factor: k_factor,
                    rotation_saturation: rotation_saturation,
                },
                internal: WindParticleInternalParameters {
                    rotation_saturation_2: rotation_saturation.powi(2),
                },
                output: WindParticleOutputParameters {
                    factor: 0.,
                },
            },
        }
    }
}

pub fn calculate_wind_factor(particles: &mut [Particle]) {
    for particle in particles.iter_mut() {
        if particle.wind.effect == WindEffect::Interaction && particle.wind.parameters.input.k_factor != 0. {
            let threshold = (particle.norm_spin_vector_2).sqrt();
            let factor = - 1. / (particle.moment_of_inertia);
            if threshold >= particle.wind.parameters.input.rotation_saturation {
                particle.wind.parameters.output.factor = factor * particle.wind.parameters.input.k_factor * particle.wind.parameters.internal.rotation_saturation_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt()
            } else {
                particle.wind.parameters.output.factor = factor * particle.wind.parameters.input.k_factor * (particle.radius/R_SUN * 1./particle.mass).sqrt()
            }
        };
    }
}
