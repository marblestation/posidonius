use serde::{Serialize, Deserialize};
use super::super::super::{Particle};
use super::super::super::{Axes};
use super::constant_time_lag;
use super::creep_coplanar;
use super::kaula;

// For future reference:
// Within common.rs, "tidal_host_particle" == "host_particle"
// not to be confused with the "host_particle" in kaula.rs, that "host_particle" correspond to the
// object that has tides on itself (tidally perturbed body).

// Future improvement (not in order of priority):
// 1: Read the data once and for all instead of reading it each time step, especially in RUST it is
//    expensive to create or read an array because RUST needs to check the existence of the array/element first.
// 2: "Common.rs" is the starting point of calculation of tidal forces, with you specific choice of
//    tidal model. The original structure of the old tidal models (constant time lag & creep)
//    take into account the both stellar tide and planetary tide with the tidal_model.rs. However,
//    this is not the case for kaula.rs, which is also why we have the bool control parameter called
//    "central_body". This parameter tells the kaula.rs which kind of tides we are calculating (and
//    thus different parameters are needed) and store the forces at "tides.parameters.output.acceleration"
//    under both star (tidal_host_particle) & (particle).
//    It will be better if we can unify the structure of each tidal model such that we can have a neat
//    structure in common.rs.
//    P.S. The idea of having "central_body" comes from the old tidal models, check it for better usage
// 3: The calculation of stellar tide in Kaula.rs requires the planetary tide (for example CTL) to be
//    switched on, thus we set those CTL parameters to be zeros (forces = 0).
//    This is probably due to malfunction of initialization of certain parameters (e.g. coordinates)
//    , so the model calculates a butch of "NaN" forces (essentially zeros).
//    If this can be fixed then for sure we can speed up the code, since we are calculating zeros.
// 4: Currently, one can only include 10 objects in a simulation, due to the limitation of the size
//    of an array in RUST. In RUST, maximum size of array is 32, but we have some arrays that take size
//    of max_number_of_planets*3, thus 10 objects (1 star and 9 planets). Sergi amazingly fixed this
//    in HIS VERSION in a day using some kind of community-developed external module. However, his version
//    needs to be merged as well.
// 5: K2 interpolation: should be a faster way instead of scanning many columns, since the data is
//    split into arrays of sizes in 32 each.
//    Discussion with Tim: Binary Tree data structure or read only certain columns.
//

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
    Kaula(kaula::KaulaParameters),
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
            TidesEffect::CentralBody(ref tidal_model) | TidesEffect::OrbitingBody(ref tidal_model) => {
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
pub fn calculate_dangular_momentum_dt_due_to_tides(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor = -1.0;

    let central_body = false;
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(tidal_model) = &particle.tides.effect {
            let torque_due_to_tides = match tidal_model {
                TidalModel::ConstantTimeLag(_) => constant_time_lag::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
                TidalModel::CreepCoplanar(_) => creep_coplanar::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
                TidalModel::Kaula(_) => kaula::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
            };
            // Integration of the spin (total torque tides):
            particle.tides.parameters.output.dangular_momentum_dt.x = factor * torque_due_to_tides.x;
            particle.tides.parameters.output.dangular_momentum_dt.y = factor * torque_due_to_tides.y;
            particle.tides.parameters.output.dangular_momentum_dt.z = factor * torque_due_to_tides.z;
        }
    }

    let central_body = true;
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::CentralBody(tidal_model) = &tidal_host_particle.tides.effect {
            let torque_due_to_tides = match tidal_model {
                TidalModel::ConstantTimeLag(_) => constant_time_lag::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
                TidalModel::CreepCoplanar(_) => creep_coplanar::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
                TidalModel::Kaula(_) => kaula::calculate_torque_due_to_tides(tidal_host_particle, particle, central_body),
            };
            // Integration of the spin (total torque tides):
            dangular_momentum_dt.x += factor * torque_due_to_tides.x;
            dangular_momentum_dt.y += factor * torque_due_to_tides.y;
            dangular_momentum_dt.z += factor * torque_due_to_tides.z;
        }
    }
    // - Equation 25 from Bolmont et al. 2015
    tidal_host_particle.tides.parameters.output.dangular_momentum_dt.x = dangular_momentum_dt.x;
    tidal_host_particle.tides.parameters.output.dangular_momentum_dt.y = dangular_momentum_dt.y;
    tidal_host_particle.tides.parameters.output.dangular_momentum_dt.z = dangular_momentum_dt.z;

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
    // Leon: I am not sure if we can just use the same formula to calculate the tidal heating for
    // the Kaula model, especially for stellar tide.
    // Nevertheless, assuming we can, we missed some key components such as
    // "radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass".
    // For now we ignored it cause tidal heating within a star is nothing, although I am not sure
    // ignoring the energy part would cause any problem or not. Personally, I would include that for
    // the completeness of the integration and also the conservation of angular momentum and energy
    // (since denergy_dt depends on dangular_momentum_dt too).
    //
    // let s_tide: bool = false;
    // if s_tide{
    //     for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
    //         if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
    //             // - Equation 32 from Bolmont et al. 2015
    //             //// Instantaneous energy loss dE/dt due to tides
    //             //// in Msun.AU^2.day^(-3)
    //             //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
    //             let factor2 = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance;
    //
    //             tidal_host_particle.tides.parameters.internal.denergy_dt = -((1.0 / particle.tides.parameters.internal.distance * (particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.tides.parameters.internal.radial_velocity))
    //                 * (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x + particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y + particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z)
    //                 + factor2
    //                 * ((particle.spin.y*particle.tides.coordinates.position.z - particle.spin.z*particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x) * particle.tides.coordinates.velocity.x
    //                 + (particle.spin.z*particle.tides.coordinates.position.x - particle.spin.x*particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y) * particle.tides.coordinates.velocity.y
    //                 + (particle.spin.x*particle.tides.coordinates.position.y - particle.spin.y*particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z) * particle.tides.coordinates.velocity.z))
    //                 - (particle.tides.parameters.output.dangular_momentum_dt.x*particle.spin.x + particle.tides.parameters.output.dangular_momentum_dt.y*particle.spin.y + particle.tides.parameters.output.dangular_momentum_dt.z*particle.spin.z);
    //         }
    //     }
    // }
}


pub fn calculate_tidal_acceleration(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor2 = 1. / tidal_host_particle.mass;
    let mut sum_tidal_force = Axes{x:0., y:0., z:0.};

    let central_body = false;
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(tidal_model) = &particle.tides.effect {
            let tidal_force = match tidal_model {
                TidalModel::ConstantTimeLag(_) => constant_time_lag::calculate_tidal_force(tidal_host_particle, particle),
                TidalModel::CreepCoplanar(_) => creep_coplanar::calculate_tidal_force(tidal_host_particle, particle),
                TidalModel::Kaula(_) => kaula::calculate_tidal_force(tidal_host_particle, particle, central_body),
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

    //// TODO: Reconsider how to move this stellar tides calculation into calculate_tidal_force while allowing a mixture of tidal models
    // Stellar tides begin
    let central_body = true;

    if matches!(
        tidal_host_particle.tides.effect,
        TidesEffect::CentralBody(TidalModel::Kaula(_))
    ) {
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            let tidal_force = kaula::calculate_tidal_force(particle, tidal_host_particle, central_body);
            let factor1 = 1. / particle.mass;

            // This code is similar to the code above, except that the acceleration is being added up
            // thus the "+=" sign instead of "="
            tidal_host_particle.tides.parameters.output.acceleration.x += factor2 * tidal_force.x;
            tidal_host_particle.tides.parameters.output.acceleration.y += factor2 * tidal_force.y;
            tidal_host_particle.tides.parameters.output.acceleration.z += factor2 * tidal_force.z;

            particle.tides.parameters.output.acceleration.x += -1.0 * factor1 * tidal_force.x;
            particle.tides.parameters.output.acceleration.y += -1.0 * factor1 * tidal_force.y;
            particle.tides.parameters.output.acceleration.z += -1.0 * factor1 * tidal_force.z;
        }
    }
}


