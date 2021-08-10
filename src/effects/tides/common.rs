use super::super::super::{Particle};
use super::super::super::{Axes};
use super::constant_time_lag;
use super::creep_coplanar;


#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleInternalParameters {
    pub distance: f64,
    pub radial_velocity: f64,
    // Constant time lag specific:
    pub scaled_dissipation_factor: f64, // sigma (dissipation_factor_scale*dissipation_factor)
    pub scalar_product_of_vector_position_with_stellar_spin: f64,
    pub scalar_product_of_vector_position_with_planetary_spin: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: f64, // Needed to compute denergy_dt
    //
    // Creep coplanar specific:
    pub shape: Axes,
    //
    pub denergy_dt: f64, // Only for history output
    pub lag_angle: f64, // Used by EvolutionType::BolmontMathis2016, EvolutionType::GalletBolmont2017 and EvolutionType::LeconteChabrier2013(true)
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleParameters {
    pub internal: TidesParticleInternalParameters,
    pub output: TidesParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleCoordinates {
    // Positions/velocities in a heliocentric frame 
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum TidalModel {
    ConstantTimeLag(constant_time_lag::ConstantTimeLagParameters),
    CreepCoplanar(creep_coplanar::CreepCoplanarParameters),
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum TidesEffect {
    CentralBody(TidalModel),
    OrbitingBody(TidalModel),
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Tides {
    pub effect: TidesEffect,
    pub parameters: TidesParticleParameters,
    pub coordinates: TidesParticleCoordinates,
}

impl Tides {
    pub fn new(effect: TidesEffect) -> Tides {
        let scaled_dissipation_factor = match effect {
            TidesEffect::CentralBody(tidal_model) | TidesEffect::OrbitingBody(tidal_model) => {
                match tidal_model {
                    TidalModel::ConstantTimeLag(params) => params.dissipation_factor_scale * params.dissipation_factor,
                    _ => 0.,
                }
            },
            _ => 0.,
        };
        Tides {
            effect: effect,
            parameters: TidesParticleParameters {
                internal: TidesParticleInternalParameters {
                    distance: 0.,
                    radial_velocity: 0.,
                    scaled_dissipation_factor: scaled_dissipation_factor,
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide: 0.,
                    orthogonal_component_of_the_tidal_force_due_to_planetary_tide: 0.,
                    radial_component_of_the_tidal_force: 0.,
                    radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
                    shape: Axes{x: 0., y: 0., z: 0.},
                    denergy_dt: 0., // Only for history output
                    lag_angle: 0., // It will be initialized the first time the evolver is called
                },
                output: TidesParticleOutputParameters {
                    acceleration: Axes{x: 0., y: 0., z: 0.},
                    dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
                },
            },
            coordinates: TidesParticleCoordinates {
                position: Axes{x: 0., y: 0., z: 0.},
                velocity: Axes{x: 0., y: 0., z: 0.},
            },
        }
    }
}

pub fn initialize(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let TidesEffect::CentralBody(_) = host_particle.tides.effect {
        host_particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = 0.;
        host_particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = 0.;
        host_particle.tides.parameters.output.acceleration.x = 0.;
        host_particle.tides.parameters.output.acceleration.y = 0.;
        host_particle.tides.parameters.output.acceleration.z = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
                particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = particle.tides.coordinates.position.x * host_particle.spin.x 
                                + particle.tides.coordinates.position.y * host_particle.spin.y
                                + particle.tides.coordinates.position.z * host_particle.spin.z;
                particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = particle.tides.coordinates.position.x * particle.spin.x 
                                + particle.tides.coordinates.position.y * particle.spin.y
                                + particle.tides.coordinates.position.z * particle.spin.z;
                particle.tides.parameters.output.acceleration.x = 0.;
                particle.tides.parameters.output.acceleration.y = 0.;
                particle.tides.parameters.output.acceleration.z = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let TidesEffect::CentralBody(_) = host_particle.tides.effect {
        // Inertial to Heliocentric positions/velocities
        host_particle.tides.coordinates.position.x = 0.;
        host_particle.tides.coordinates.position.y = 0.;
        host_particle.tides.coordinates.position.z = 0.;
        host_particle.tides.coordinates.velocity.x = 0.;
        host_particle.tides.coordinates.velocity.y = 0.;
        host_particle.tides.coordinates.velocity.z = 0.;
        host_particle.tides.parameters.internal.distance = 0.;
        host_particle.tides.parameters.internal.radial_velocity = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
                particle.tides.coordinates.position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                particle.tides.coordinates.position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                particle.tides.coordinates.position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                particle.tides.coordinates.velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.tides.coordinates.velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.tides.coordinates.velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.tides.parameters.internal.distance = (particle.tides.coordinates.position.x.powi(2) 
                                                                + particle.tides.coordinates.position.y.powi(2)
                                                                + particle.tides.coordinates.position.z.powi(2)).sqrt();
                particle.tides.parameters.internal.radial_velocity = (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x +
                                                                        particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y +
                                                                        particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z) / particle.tides.parameters.internal.distance;
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let TidesEffect::CentralBody(_) = host_particle.tides.effect {
        host_particle.tides.coordinates.position = host_particle.heliocentric_position;
        host_particle.tides.coordinates.velocity = host_particle.heliocentric_velocity;
        host_particle.tides.parameters.internal.distance = host_particle.heliocentric_distance;
        host_particle.tides.parameters.internal.radial_velocity = host_particle.heliocentric_radial_velocity;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
                particle.tides.coordinates.position = particle.heliocentric_position;
                particle.tides.coordinates.velocity = particle.heliocentric_velocity;
                particle.tides.parameters.internal.distance = particle.heliocentric_distance;
                particle.tides.parameters.internal.radial_velocity = particle.heliocentric_radial_velocity;
            }
        }
    }
}



//////////////////////////////////////////////////////////////////////////////
//// TIDES
pub fn calculate_dangular_momentum_dt_due_to_tides(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool) {
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(tidal_model) = particle.tides.effect {
            let torque_due_to_tides = match tidal_model {
                TidalModel::ConstantTimeLag(_) => constant_time_lag::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
                TidalModel::CreepCoplanar(_) => creep_coplanar::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
            };
            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides.x;
                dangular_momentum_dt.y += factor * torque_due_to_tides.y;
                dangular_momentum_dt.z += factor * torque_due_to_tides.z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x = factor * torque_due_to_tides.x;
                particle.tides.parameters.output.dangular_momentum_dt.y = factor * torque_due_to_tides.y;
                particle.tides.parameters.output.dangular_momentum_dt.z = factor * torque_due_to_tides.z;
            }
        }
    }

    if central_body {
        // - Equation 25 from Bolmont et al. 2015
        tidal_host_particle.tides.parameters.output.dangular_momentum_dt.x = dangular_momentum_dt.x;
        tidal_host_particle.tides.parameters.output.dangular_momentum_dt.y = dangular_momentum_dt.y;
        tidal_host_particle.tides.parameters.output.dangular_momentum_dt.z = dangular_momentum_dt.z;
    }

}

pub fn calculate_denergy_dt(particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
            // - Equation 32 from Bolmont et al. 2015
            //// Instantaneous energy loss dE/dt due to tides
            //// in Msun.AU^2.day^(-3)
            //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
            let factor2 = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance;
            particle.tides.parameters.internal.denergy_dt = -((1.0 / particle.tides.parameters.internal.distance * (particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.tides.parameters.internal.radial_velocity))
                        * (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x + particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y + particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z)
                        + factor2 
                        * ((particle.spin.y*particle.tides.coordinates.position.z - particle.spin.z*particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x) * particle.tides.coordinates.velocity.x
                        + (particle.spin.z*particle.tides.coordinates.position.x - particle.spin.x*particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y) * particle.tides.coordinates.velocity.y
                        + (particle.spin.x*particle.tides.coordinates.position.y - particle.spin.y*particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z) * particle.tides.coordinates.velocity.z))
                        - (particle.tides.parameters.output.dangular_momentum_dt.x*particle.spin.x + particle.tides.parameters.output.dangular_momentum_dt.y*particle.spin.y + particle.tides.parameters.output.dangular_momentum_dt.z*particle.spin.z);
        }
    }
}


pub fn calculate_tidal_acceleration(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor2 = 1. / tidal_host_particle.mass;
    let mut sum_tidal_force = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(tidal_model) = particle.tides.effect {
            let tidal_force = match tidal_model {
                TidalModel::ConstantTimeLag(_) => constant_time_lag::calculate_tidal_force(tidal_host_particle, particle),
                TidalModel::CreepCoplanar(_) => creep_coplanar::calculate_tidal_force(tidal_host_particle, particle),
            };
            let factor1 = 1. / particle.mass;
            sum_tidal_force.x += tidal_force.x;
            sum_tidal_force.y += tidal_force.y;
            sum_tidal_force.z += tidal_force.z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * tidal_force.x; 
            particle.tides.parameters.output.acceleration.y = factor1 * tidal_force.y;
            particle.tides.parameters.output.acceleration.z = factor1 * tidal_force.z;
        }
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.tides.parameters.output.acceleration.x += factor2 * sum_tidal_force.x;
        //particle.tides.parameters.output.acceleration.y += factor2 * sum_tidal_force.y;
        //particle.tides.parameters.output.acceleration.z += factor2 * sum_tidal_force.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    tidal_host_particle.tides.parameters.output.acceleration.x = -1.0 * factor2 * sum_tidal_force.x;
    tidal_host_particle.tides.parameters.output.acceleration.y = -1.0 * factor2 * sum_tidal_force.y;
    tidal_host_particle.tides.parameters.output.acceleration.z = -1.0 * factor2 * sum_tidal_force.z;
}


