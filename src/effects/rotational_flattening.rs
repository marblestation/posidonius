use super::super::{Particle};
use super::super::{Axes};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleInputParameters {
    pub love_number: f64,   // love number for a completely fluid planet (used for rotational flattening effects)
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleInternalParameters {
    pub distance: f64,
    pub scalar_product_of_vector_position_with_stellar_spin: f64,
    pub scalar_product_of_vector_position_with_planetary_spin: f64,
    pub radial_component_of_the_force_induced_by_rotation: f64,
    pub factor_for_the_force_induced_by_star_rotation: f64,
    pub factor_for_the_force_induced_by_planet_rotation: f64,
    pub orthogonal_component_of_the_force_induced_by_star_rotation: f64,
    pub orthogonal_component_of_the_force_induced_by_planet_rotation: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleParameters {
    pub input: RotationalFlatteningParticleInputParameters,
    pub internal: RotationalFlatteningParticleInternalParameters,
    pub output: RotationalFlatteningParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleCoordinates {
    // Positions/velocities in a heliocentric frame 
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum RotationalFlatteningEffect {
    CentralBody,
    OrbitingBody,
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlattening {
    pub effect: RotationalFlatteningEffect,
    pub parameters: RotationalFlatteningParticleParameters,
    pub coordinates: RotationalFlatteningParticleCoordinates,
}

impl RotationalFlattening {
    pub fn new(effect: RotationalFlatteningEffect, love_number: f64) -> RotationalFlattening {
        RotationalFlattening {
            effect: effect,
            parameters: RotationalFlatteningParticleParameters {
                input: RotationalFlatteningParticleInputParameters {
                    love_number: love_number,
                },
                internal: RotationalFlatteningParticleInternalParameters {
                    distance: 0.,
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    radial_component_of_the_force_induced_by_rotation: 0.,
                    factor_for_the_force_induced_by_star_rotation: 0.,
                    factor_for_the_force_induced_by_planet_rotation: 0.,
                    orthogonal_component_of_the_force_induced_by_star_rotation: 0.,
                    orthogonal_component_of_the_force_induced_by_planet_rotation: 0.,
                },
                output: RotationalFlatteningParticleOutputParameters {
                    acceleration: Axes{x: 0., y: 0., z: 0.},
                    dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
                },
            },
            coordinates: RotationalFlatteningParticleCoordinates {
                position: Axes{x: 0., y: 0., z: 0.},
                velocity: Axes{x: 0., y: 0., z: 0.},
            },
        }
    }
}

pub fn initialize(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let RotationalFlatteningEffect::CentralBody = host_particle.rotational_flattening.effect {
        host_particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = 0.;
        host_particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = 0.;
        host_particle.rotational_flattening.parameters.output.acceleration.x = 0.;
        host_particle.rotational_flattening.parameters.output.acceleration.y = 0.;
        host_particle.rotational_flattening.parameters.output.acceleration.z = 0.;
        host_particle.rotational_flattening.parameters.output.dangular_momentum_dt.x = 0.;
        host_particle.rotational_flattening.parameters.output.dangular_momentum_dt.y = 0.;
        host_particle.rotational_flattening.parameters.output.dangular_momentum_dt.z = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let RotationalFlatteningEffect::OrbitingBody = particle.rotational_flattening.effect {
                particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = particle.rotational_flattening.coordinates.position.x * host_particle.spin.x 
                                                                                    + particle.rotational_flattening.coordinates.position.y * host_particle.spin.y
                                                                                    + particle.rotational_flattening.coordinates.position.z * host_particle.spin.z;
                particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = particle.rotational_flattening.coordinates.position.x * particle.spin.x 
                                                                                    + particle.rotational_flattening.coordinates.position.y * particle.spin.y
                                                                                    + particle.rotational_flattening.coordinates.position.z * particle.spin.z;
                particle.rotational_flattening.parameters.output.acceleration.x = 0.;
                particle.rotational_flattening.parameters.output.acceleration.y = 0.;
                particle.rotational_flattening.parameters.output.acceleration.z = 0.;
                particle.rotational_flattening.parameters.output.dangular_momentum_dt.x = 0.;
                particle.rotational_flattening.parameters.output.dangular_momentum_dt.y = 0.;
                particle.rotational_flattening.parameters.output.dangular_momentum_dt.z = 0.;
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let RotationalFlatteningEffect::CentralBody = host_particle.rotational_flattening.effect {
        // Inertial to Heliocentric positions/velocities
        host_particle.rotational_flattening.coordinates.position.x = 0.;
        host_particle.rotational_flattening.coordinates.position.y = 0.;
        host_particle.rotational_flattening.coordinates.position.z = 0.;
        host_particle.rotational_flattening.coordinates.velocity.x = 0.;
        host_particle.rotational_flattening.coordinates.velocity.y = 0.;
        host_particle.rotational_flattening.coordinates.velocity.z = 0.;
        host_particle.rotational_flattening.parameters.internal.distance = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let RotationalFlatteningEffect::OrbitingBody = particle.rotational_flattening.effect {
                particle.rotational_flattening.coordinates.position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                particle.rotational_flattening.coordinates.position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                particle.rotational_flattening.coordinates.position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                particle.rotational_flattening.coordinates.velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.rotational_flattening.coordinates.velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.rotational_flattening.coordinates.velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.rotational_flattening.parameters.internal.distance = (particle.rotational_flattening.coordinates.position.x.powi(2) 
                                                                + particle.rotational_flattening.coordinates.position.y.powi(2)
                                                                + particle.rotational_flattening.coordinates.position.z.powi(2)).sqrt();
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let RotationalFlatteningEffect::CentralBody = host_particle.rotational_flattening.effect {
        host_particle.rotational_flattening.coordinates.position = host_particle.heliocentric_position;
        host_particle.rotational_flattening.coordinates.velocity = host_particle.heliocentric_velocity;
        host_particle.rotational_flattening.parameters.internal.distance = host_particle.heliocentric_distance;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let RotationalFlatteningEffect::OrbitingBody = particle.rotational_flattening.effect {
                particle.rotational_flattening.coordinates.position = particle.heliocentric_position;
                particle.rotational_flattening.coordinates.velocity = particle.heliocentric_velocity;
                particle.rotational_flattening.parameters.internal.distance = particle.heliocentric_distance;
            }
        }
    }
}

pub fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(rotational_flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let mut rotational_flattening_host_particle = rotational_flattening_host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let central_body = true;
    calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(central_body, &mut rotational_flattening_host_particle, &mut particles, &mut more_particles);
    calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(!central_body, &mut rotational_flattening_host_particle, &mut particles, &mut more_particles);
}

fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(central_body:bool, rotational_flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    // Calculation of the norm square of the spin for the planet
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if central_body {
            // - Star Equation 16 from Bolmont et al. 2015
            particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation = particle.mass * rotational_flattening_host_particle.rotational_flattening.parameters.input.love_number * rotational_flattening_host_particle.norm_spin_vector_2 * rotational_flattening_host_particle.radius.powi(5) / 6.; // Msun.AU^5.day-2
            // - Second part of Equation 15 from Bolmont et al. 2015
            particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_star_rotation = -6. * particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin / (rotational_flattening_host_particle.norm_spin_vector_2 * particle.rotational_flattening.parameters.internal.distance.powi(5));
        } else {
            // - Planet Equation 16 from Bolmont et al. 2015
            particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation = rotational_flattening_host_particle.mass * particle.rotational_flattening.parameters.input.love_number * particle.norm_spin_vector_2 * particle.radius.powi(5) / 6.; // Msun.AU^5.day-2
            // - Second part of Equation 15 from Bolmont et al. 2015
            particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_planet_rotation = -6. * particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin / (particle.norm_spin_vector_2 * particle.rotational_flattening.parameters.internal.distance.powi(5));
        }
    }
}

pub fn calculate_radial_component_of_the_force_induced_by_rotational_flattening(rotational_flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // - First part of Equation 15 from Bolmont et al. 2015
        particle.rotational_flattening.parameters.internal.radial_component_of_the_force_induced_by_rotation = -3./particle.rotational_flattening.parameters.internal.distance.powi(5) * (particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation + particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation)
            + 15./particle.rotational_flattening.parameters.internal.distance.powi(7) * (particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin/rotational_flattening_host_particle.norm_spin_vector_2
                + particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin/particle.norm_spin_vector_2); // Msun.AU.day-1
    }
}

pub fn calculate_torque_induced_by_rotational_flattening(rotational_flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool) {
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};
    let mut reference_spin = rotational_flattening_host_particle.spin.clone();
    let mut orthogonal_component_of_the_force_induced_by_rotation: f64;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if !central_body {
            reference_spin = particle.spin.clone();
            orthogonal_component_of_the_force_induced_by_rotation = particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_planet_rotation;
        } else {
            orthogonal_component_of_the_force_induced_by_rotation = particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_star_rotation;
        }
        
        //// Torque calculation due to rotational flattening
        // - Equation 17-18 from Bolmont et al. 2015
        let torque_induced_by_rotational_flattening_x: f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.rotational_flattening.coordinates.position.y * reference_spin.z - particle.rotational_flattening.coordinates.position.z * reference_spin.y);
        let torque_induced_by_rotational_flattening_y :f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.rotational_flattening.coordinates.position.z * reference_spin.x - particle.rotational_flattening.coordinates.position.x * reference_spin.z);
        let torque_induced_by_rotational_flattening_z: f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.rotational_flattening.coordinates.position.x * reference_spin.y - particle.rotational_flattening.coordinates.position.y * reference_spin.x);

        let factor = -1.0;
        // - Equation 25 from Bolmont et al. 2015
        if central_body {
            // Integration of the spin (total torque rot):
            dangular_momentum_dt.x += factor * torque_induced_by_rotational_flattening_x;
            dangular_momentum_dt.y += factor * torque_induced_by_rotational_flattening_y;
            dangular_momentum_dt.z += factor * torque_induced_by_rotational_flattening_z;
        } else {
            particle.rotational_flattening.parameters.output.dangular_momentum_dt.x = factor * torque_induced_by_rotational_flattening_x;
            particle.rotational_flattening.parameters.output.dangular_momentum_dt.y = factor * torque_induced_by_rotational_flattening_y;
            particle.rotational_flattening.parameters.output.dangular_momentum_dt.z = factor * torque_induced_by_rotational_flattening_z;
        }
    }

    if central_body {
        // - Equation 25 from Bolmont et al. 2015
        rotational_flattening_host_particle.rotational_flattening.parameters.output.dangular_momentum_dt.x = dangular_momentum_dt.x;
        rotational_flattening_host_particle.rotational_flattening.parameters.output.dangular_momentum_dt.y = dangular_momentum_dt.y;
        rotational_flattening_host_particle.rotational_flattening.parameters.output.dangular_momentum_dt.z = dangular_momentum_dt.z;
    }
}

pub fn calculate_acceleration_induced_by_rotational_flattering(rotational_flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor2 = 1. / rotational_flattening_host_particle.mass;

    // Calculation of the norm square of the spin for the star
    let mut sum_total_force_induced_by_rotation = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        let factor1 = 1. / particle.mass;

        // - Equation 15 from Bolmont et al. 2015
        let total_force_induced_by_rotation_x = particle.rotational_flattening.parameters.internal.radial_component_of_the_force_induced_by_rotation * particle.rotational_flattening.coordinates.position.x
            + particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.x
            + particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_star_rotation * rotational_flattening_host_particle.spin.x;
        let total_force_induced_by_rotation_y = particle.rotational_flattening.parameters.internal.radial_component_of_the_force_induced_by_rotation * particle.rotational_flattening.coordinates.position.y
            + particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.y
            + particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_star_rotation * rotational_flattening_host_particle.spin.y;
        let total_force_induced_by_rotation_z = particle.rotational_flattening.parameters.internal.radial_component_of_the_force_induced_by_rotation * particle.rotational_flattening.coordinates.position.z
            + particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.z
            + particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_star_rotation * rotational_flattening_host_particle.spin.z;
        //println!("Rot force = {:e} {:e} {:e}", total_force_induced_by_rotation_x, total_force_induced_by_rotation_y, total_force_induced_by_rotation_z);

        sum_total_force_induced_by_rotation.x += total_force_induced_by_rotation_x;
        sum_total_force_induced_by_rotation.y += total_force_induced_by_rotation_y;
        sum_total_force_induced_by_rotation.z += total_force_induced_by_rotation_z;
          
        // - Equation 19 from Bolmont et al. 2015 (first term)
        particle.rotational_flattening.parameters.output.acceleration.x = factor1 * total_force_induced_by_rotation_x; 
        particle.rotational_flattening.parameters.output.acceleration.y = factor1 * total_force_induced_by_rotation_y;
        particle.rotational_flattening.parameters.output.acceleration.z = factor1 * total_force_induced_by_rotation_z;
    }
    
    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.rotational_flattening.parameters.output.acceleration.x += factor2 * sum_total_force_induced_by_rotation.x;
        //particle.rotational_flattening.parameters.output.acceleration.y += factor2 * sum_total_force_induced_by_rotation.y;
        //particle.rotational_flattening.parameters.output.acceleration.z += factor2 * sum_total_force_induced_by_rotation.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    rotational_flattening_host_particle.rotational_flattening.parameters.output.acceleration.x = -1.0 * factor2 * sum_total_force_induced_by_rotation.x;
    rotational_flattening_host_particle.rotational_flattening.parameters.output.acceleration.y = -1.0 * factor2 * sum_total_force_induced_by_rotation.y;
    rotational_flattening_host_particle.rotational_flattening.parameters.output.acceleration.z = -1.0 * factor2 * sum_total_force_induced_by_rotation.z;
}

