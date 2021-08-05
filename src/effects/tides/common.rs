use std::collections::HashMap;
use super::super::tools;
use super::super::constants::{K2, SMOOTHING_FACTOR_DYN_TIDE_COROTATION};
use super::super::{Particle};
use super::super::{Axes};
use super::{EvolutionType};


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
pub struct ConstantTimeLagParameters {
    pub dissipation_factor: f64,
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub love_number: f64,   // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the 
                            // susceptibility of its shape to change in response to a tidal potential.
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct CreepCoplanarParameters {
    pub uniform_viscosity_coefficient: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum TidalModel {
    ConstantTimeLag(ConstantTimeLagParameters),
    CreepCoplanar(CreepCoplanarParameters),
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


pub fn calculate_planet_dependent_dissipation_factors(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    match tidal_host_particle.evolution {
        EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
            let star_norm_spin_vector = tidal_host_particle.norm_spin_vector_2.sqrt();
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let TidesEffect::OrbitingBody(tidal_model) = particle.tides.effect {
                    if let TidalModel::ConstantTimeLag(_) = tidal_model {
                        let (tidal_host_particle_dissipation_factor_scale, tidal_host_particle_dissipation_factor) = match tidal_host_particle.tides.effect {
                            TidesEffect::CentralBody(tidal_model) => {
                                match tidal_model {
                                    TidalModel::ConstantTimeLag(params) => (params.dissipation_factor_scale, params.dissipation_factor),
                                    _ => (0., 0.)
                                }
                            },
                            _ => (0., 0.)
                        };
                        //
                        //// Excitation frequency needed by the model based on the
                        // instantaneous frequency (using positions, velocities and spins)
                        //let frequency = (particle.tides.coordinates.velocity.x - tidal_host_particle.spin.y*particle.tides.coordinates.position.z + tidal_host_particle.spin.z*particle.tides.coordinates.position.y).powi(2)
                                    //+ (particle.tides.coordinates.velocity.y - tidal_host_particle.spin.z*particle.tides.coordinates.position.x + tidal_host_particle.spin.x*particle.tides.coordinates.position.z).powi(2)
                                    //+ (particle.tides.coordinates.velocity.z - tidal_host_particle.spin.x*particle.tides.coordinates.position.y + tidal_host_particle.spin.y*particle.tides.coordinates.position.x).powi(2);
                        //let inverse_of_half_the_excitation_frequency = particle.tides.parameters.internal.distance / frequency;
                        // NOTE:  two_times_the_inverse_of_the_excitation_frequency: 2/w
                        //        inverse_of_half_the_excitation_frequency : 1/(w/2)
                        //
                        //// Excitation frequency needed by the model based on the
                        // mean frequency (using mean motion and spin). 
                        //
                        // NOTE: The model is already here being used outside the 
                        // validity domain, it seems not justified to use an 
                        // instantaneous frequency.
                        let gm = tidal_host_particle.mass_g+particle.mass_g;
                        let (perihelion_distance, eccentricity) = tools::calculate_perihelion_distance_and_eccentricity(gm, particle.tides.coordinates.position, particle.tides.coordinates.velocity);
                        let mean_motion = gm.sqrt() * (perihelion_distance/(1.0 - eccentricity)).powf(-1.5);
                        let mut half_the_excitation_frequency = (star_norm_spin_vector - mean_motion).abs();
                        // If the dynamical tide is excited, compute the planet-dependent stellar dissipation
                        // If the dynamical tide is not excited (i.e., equilibrium tide), do nothing since the default value corresponds to the planet scaled dissipation factor
                        if half_the_excitation_frequency < star_norm_spin_vector { 
                            // The dynamical tide is excited
                            if half_the_excitation_frequency < SMOOTHING_FACTOR_DYN_TIDE_COROTATION {
                                half_the_excitation_frequency = SMOOTHING_FACTOR_DYN_TIDE_COROTATION;
                            }
                            let inverse_of_half_the_excitation_frequency = 1./half_the_excitation_frequency;

                            // Eq. 4 and Eq.10 of Bolmont & Mathis have a typo, see page 5 of Gallet &
                            // Bolmont 2017. The lag angle has a 1/k2, but here we should *k2, that is why
                            // k2 does not appear here
                            // We add here the dissipation from the equilibrium tide to the dynamical tide one
                            let planet_dependent_dissipation_factor = tidal_host_particle_dissipation_factor_scale
                                * (2.0 * K2 / (3.0*tidal_host_particle.radius.powi(5))
                                * tidal_host_particle.tides.parameters.internal.lag_angle * inverse_of_half_the_excitation_frequency
                                + tidal_host_particle_dissipation_factor);

                            star_planet_dependent_dissipation_factors.insert(particle.id.clone(), planet_dependent_dissipation_factor);
                        }
                    }
                }
            }
            //panic!("Please, contact Posidonius authors before using BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true) evolutionary models. They may not be ready yet for scientific explotation.")
        },
        _ => {},
    }

}

pub fn planet_dependent_dissipation_factor(star_planet_dependent_dissipation_factors: &HashMap<usize, f64>,  id: &usize, evolution: EvolutionType, scaled_dissipation_factor: f64) -> f64 {
    match evolution {
        EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
            match star_planet_dependent_dissipation_factors.get(id) {
                Some(&value) => value,
                _ => scaled_dissipation_factor // This should not happen
            }
        },
        _ => scaled_dissipation_factor,
    }
}


//////////////////////////////////////////////////////////////////////////////
//// TIDES
pub fn calculate_torque_due_to_tides(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool) {
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};
    let mut reference_spin = tidal_host_particle.spin.clone();
    let mut orthogonal_component_of_the_tidal_force: f64;
    let mut reference_rscalspin: f64;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
            if !central_body {
                reference_spin = particle.spin.clone();
                reference_rscalspin = particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin;
                orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
            } else {
                reference_rscalspin = particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin;
                orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
            }

            // distance to star
            let distance = particle.tides.parameters.internal.distance;

            //// Torque calculation (star)
            // - Equation 8-9 from Bolmont et al. 2015
            let torque_due_to_tides_x: f64 = orthogonal_component_of_the_tidal_force 
                            * (distance * reference_spin.x - reference_rscalspin*particle.tides.coordinates.position.x/distance - 1.0/distance
                            * (particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.z - particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.y) );

            let torque_due_to_tides_y :f64 = orthogonal_component_of_the_tidal_force 
                            * (distance * reference_spin.y - reference_rscalspin*particle.tides.coordinates.position.y/distance - 1.0/distance
                            * (particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.x - particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.z) );

            let torque_due_to_tidez_z: f64 = orthogonal_component_of_the_tidal_force 
                            * (distance * reference_spin.z - reference_rscalspin*particle.tides.coordinates.position.z/distance - 1.0/distance
                            * (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.y - particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.x) );

            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += factor * torque_due_to_tidez_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x = factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y = factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z = factor * torque_due_to_tidez_z;
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

pub fn calculate_orthogonal_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    let mut tidal_host_particle = tidal_host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let mut star_planet_dependent_dissipation_factors = star_planet_dependent_dissipation_factors;
    let central_body = true;
    calculate_orthogonal_component_of_the_tidal_force_for(central_body, &mut tidal_host_particle, &mut particles, &mut more_particles, &mut star_planet_dependent_dissipation_factors);
    calculate_orthogonal_component_of_the_tidal_force_for(!central_body, &mut tidal_host_particle, &mut particles, &mut more_particles, &mut star_planet_dependent_dissipation_factors);
}

fn calculate_orthogonal_component_of_the_tidal_force_for(central_body:bool, tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(tidal_model) = particle.tides.effect {
            if let TidalModel::ConstantTimeLag(_) = tidal_model {
                //// Only calculate tides if planet is not in disk
                //if particle.disk_interaction_time == 0.0 {

                // (distance to star)^7
                let distance_7 = particle.tides.parameters.internal.distance.powi(7);

                //// Tidal force calculation (star) :: Only orthogonal component is needed
                if central_body {
                    // - Third line of Equation 5 from Bolmont et al. 2015
                    //   This expression has R**10 (instead of R**5 in Eq. 5) 
                    //   because it uses sigma (i.e., scaled_dissipation_factor) 
                    //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
                    //   there is a R**5 factor as shown in Equation 28)
                    //   - k2 is love number
                    let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(&star_planet_dependent_dissipation_factors, &particle.id, tidal_host_particle.evolution, tidal_host_particle.tides.parameters.internal.scaled_dissipation_factor);
                    particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass.powi(2))
                                                    * (tidal_host_particle.radius.powi(10)) 
                                                    * star_scaled_dissipation_factor / distance_7;
                } else {
                    // - Second line of Equation 5 from Bolmont et al. 2015
                    //   This expression has R**10 (instead of R**5 in Eq. 5) 
                    //   because it uses sigma (i.e., scaled_dissipation_factor) 
                    //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
                    //   there is a R**5 factor as shown in Equation 28)
                    //   - k2 is love number
                    particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (tidal_host_particle.mass.powi(2))
                                                    * (particle.radius.powi(10))
                                                    * particle.tides.parameters.internal.scaled_dissipation_factor / distance_7;

                    // SBC
                    //println!("> {:e} {:e} {:e} {:e}", tidal_host_particle.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
                    //println!("> {:e} {:e} {:e}", particle.tides.coordinates.position.x, particle.tides.coordinates.position.y, particle.tides.coordinates.position.z);
                }
                //} else {
                    //if central_body {
                        //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 0.0
                    //} else{
                        //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0
                    //}
                //}
            }
        }
    }
}

pub fn calculate_radial_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    let star_mass_2 = tidal_host_particle.mass * tidal_host_particle.mass;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(tidal_model) = particle.tides.effect {
            if let TidalModel::ConstantTimeLag(params) = tidal_model {
                let tidal_host_particle_love_number = match tidal_host_particle.tides.effect {
                    TidesEffect::CentralBody(tidal_model) => {
                        match tidal_model {
                            TidalModel::ConstantTimeLag(params) => params.love_number,
                            _ => 0.
                        }
                    },
                    _ => 0.
                };
                let planet_mass_2 = particle.mass * particle.mass;
                // Conservative part of the radial tidal force
                let radial_component_of_the_tidal_force_conservative_part = -3.0 * K2 / particle.tides.parameters.internal.distance.powi(7)
                            * (planet_mass_2 * tidal_host_particle.radius.powi(5) * tidal_host_particle_love_number
                            + star_mass_2 * particle.radius.powi(5) * params.love_number);

                // Dissipative part of the radial tidal force:
                let factor1 = -13.5 * particle.tides.parameters.internal.radial_velocity / particle.tides.parameters.internal.distance.powi(8);
                let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(&star_planet_dependent_dissipation_factors, &particle.id, tidal_host_particle.evolution, tidal_host_particle.tides.parameters.internal.scaled_dissipation_factor);
                let term1 = planet_mass_2
                            * tidal_host_particle.radius.powi(10)
                            * star_scaled_dissipation_factor;
                let term2 = star_mass_2
                            * particle.radius.powi(10)
                            * particle.tides.parameters.internal.scaled_dissipation_factor;
                // If we consider the star as a point mass (used for denergy_dt calculation):
                particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass = factor1 * term2;
                let radial_component_of_the_tidal_force_dissipative_part = particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor1 * term1;

                // Sum of the dissipative and conservative part of the radial force
                // - First line Equation 5 from Bolmont et al. 2015
                particle.tides.parameters.internal.radial_component_of_the_tidal_force = radial_component_of_the_tidal_force_conservative_part + radial_component_of_the_tidal_force_dissipative_part;
            }
        }
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
    let mut sum_total_tidal_force = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
            let factor1 = 1. / particle.mass;

            // - Equation 6 from Bolmont et al. 2015
            let factor3 = particle.tides.parameters.internal.radial_component_of_the_tidal_force
                            + (particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide) * particle.tides.parameters.internal.radial_velocity / particle.tides.parameters.internal.distance;
            let total_tidal_force_x = factor3 * particle.tides.coordinates.position.x / particle.tides.parameters.internal.distance
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance
                                        * (tidal_host_particle.spin.y * particle.tides.coordinates.position.z  - tidal_host_particle.spin.z * particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x)
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance 
                                        * (particle.spin.y * particle.tides.coordinates.position.z  - particle.spin.z * particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x);
            let total_tidal_force_y = factor3 * particle.tides.coordinates.position.y / particle.tides.parameters.internal.distance
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance 
                                        * (tidal_host_particle.spin.z * particle.tides.coordinates.position.x  - tidal_host_particle.spin.x * particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y)
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance 
                                        * (particle.spin.z * particle.tides.coordinates.position.x  - particle.spin.x * particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y);
            let total_tidal_force_z = factor3 * particle.tides.coordinates.position.z / particle.tides.parameters.internal.distance
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance 
                                        * (tidal_host_particle.spin.x * particle.tides.coordinates.position.y  - tidal_host_particle.spin.y * particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z)
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance 
                                        * (particle.spin.x * particle.tides.coordinates.position.y  - particle.spin.y * particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z);
            //println!("factor3 {:e} {:e} {:e}", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
            //println!("d {:e} vrad {:e}", particle.tides.parameters.internal., particle.tides.parameters.internal.radial_velocity);
            //println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);

            sum_total_tidal_force.x += total_tidal_force_x;
            sum_total_tidal_force.y += total_tidal_force_y;
            sum_total_tidal_force.z += total_tidal_force_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x; 
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
        }
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.tides.parameters.output.acceleration.x += factor2 * sum_total_tidal_force.x;
        //particle.tides.parameters.output.acceleration.y += factor2 * sum_total_tidal_force.y;
        //particle.tides.parameters.output.acceleration.z += factor2 * sum_total_tidal_force.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    tidal_host_particle.tides.parameters.output.acceleration.x = -1.0 * factor2 * sum_total_tidal_force.x;
    tidal_host_particle.tides.parameters.output.acceleration.y = -1.0 * factor2 * sum_total_tidal_force.y;
    tidal_host_particle.tides.parameters.output.acceleration.z = -1.0 * factor2 * sum_total_tidal_force.z;
}


