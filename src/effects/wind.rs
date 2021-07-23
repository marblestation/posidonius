use super::super::constants::{R_SUN};
use super::super::{Particle};
use super::super::{Axes};

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
    pub dangular_momentum_dt: Axes, // Force
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
    Disabled,
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
                    dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
                },
            },
        }
    }
}

pub fn initialize(particles: &mut [Particle]) {
    for particle in particles.iter_mut() {
        if particle.wind.effect == WindEffect::Interaction {
            particle.wind.parameters.output.dangular_momentum_dt.x = 0.;
            particle.wind.parameters.output.dangular_momentum_dt.y = 0.;
            particle.wind.parameters.output.dangular_momentum_dt.z = 0.;
        }
    }
}

pub fn calculate_wind_factor(particles: &mut [Particle]) {
    // TODO: Verify wind factor
    for particle in particles.iter_mut() {
        if particle.wind.effect == WindEffect::Interaction {
            let threshold = particle.norm_spin_vector_2.sqrt();
            // Eq. 14 from Bolmont & Mathis 2016 (adapted to be able to add to the rest of
            // dangular_momentum_dt)
            if threshold >= particle.wind.parameters.input.rotation_saturation {
                // Fast rotator
                particle.wind.parameters.output.dangular_momentum_dt.x = -1. * particle.wind.parameters.input.k_factor * particle.spin.x * particle.wind.parameters.internal.rotation_saturation_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt();
                particle.wind.parameters.output.dangular_momentum_dt.y = -1. * particle.wind.parameters.input.k_factor * particle.spin.y * particle.wind.parameters.internal.rotation_saturation_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt();
                particle.wind.parameters.output.dangular_momentum_dt.z = -1. * particle.wind.parameters.input.k_factor * particle.spin.z * particle.wind.parameters.internal.rotation_saturation_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt();
            } else {
                // Slow rotators
                particle.wind.parameters.output.dangular_momentum_dt.x = -1. * particle.wind.parameters.input.k_factor * particle.spin.x * particle.norm_spin_vector_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt();
                particle.wind.parameters.output.dangular_momentum_dt.y = -1. * particle.wind.parameters.input.k_factor * particle.spin.y * particle.norm_spin_vector_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt();
                particle.wind.parameters.output.dangular_momentum_dt.z = -1. * particle.wind.parameters.input.k_factor * particle.spin.z * particle.norm_spin_vector_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt();
            }
        };
    }
}
