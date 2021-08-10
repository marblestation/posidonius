use super::super::super::{Particle};
use super::super::super::{Axes};
use super::common::RotationalFlatteningEffect;
use super::common::RotationalFlatteningModel;

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct OblateSpheroidParameters {
    pub love_number: f64,   // love number for a completely fluid planet (used for rotational flattening effects)
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
        if let RotationalFlatteningEffect::OrbitingBody(rotational_flattening_model) = particle.rotational_flattening.effect {
            if let RotationalFlatteningModel::OblateSpheroid(params) = rotational_flattening_model {
                let rotational_flattening_host_particle_love_number = match rotational_flattening_host_particle.rotational_flattening.effect {
                    RotationalFlatteningEffect::CentralBody(rotational_flattening_model) => {
                        match rotational_flattening_model {
                            RotationalFlatteningModel::OblateSpheroid(params) => params.love_number,
                            _ => 0.
                        }
                    },
                    _ => 0.
                };
                if central_body {
                    // - Star Equation 16 from Bolmont et al. 2015
                    particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation = particle.mass * rotational_flattening_host_particle_love_number * rotational_flattening_host_particle.norm_spin_vector_2 * rotational_flattening_host_particle.radius.powi(5) / 6.; // Msun.AU^5.day-2
                    // - Second part of Equation 15 from Bolmont et al. 2015
                    particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_star_rotation = -6. * particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin / (rotational_flattening_host_particle.norm_spin_vector_2 * particle.rotational_flattening.parameters.internal.distance.powi(5));
                } else {
                    // - Planet Equation 16 from Bolmont et al. 2015
                    particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation = rotational_flattening_host_particle.mass * params.love_number * particle.norm_spin_vector_2 * particle.radius.powi(5) / 6.; // Msun.AU^5.day-2
                    // - Second part of Equation 15 from Bolmont et al. 2015
                    particle.rotational_flattening.parameters.internal.orthogonal_component_of_the_force_induced_by_planet_rotation = -6. * particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin / (particle.norm_spin_vector_2 * particle.rotational_flattening.parameters.internal.distance.powi(5));
                }
            }
        }
    }
}

pub fn calculate_radial_component_of_the_force_induced_by_rotational_flattening(rotational_flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let RotationalFlatteningEffect::OrbitingBody(_) = particle.rotational_flattening.effect {
            // - First part of Equation 15 from Bolmont et al. 2015
            particle.rotational_flattening.parameters.internal.radial_component_of_the_force_induced_by_rotation = -3./particle.rotational_flattening.parameters.internal.distance.powi(5) * (particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation + particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation)
                + 15./particle.rotational_flattening.parameters.internal.distance.powi(7) * (particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_star_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_stellar_spin/rotational_flattening_host_particle.norm_spin_vector_2
                    + particle.rotational_flattening.parameters.internal.factor_for_the_force_induced_by_planet_rotation * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin * particle.rotational_flattening.parameters.internal.scalar_product_of_vector_position_with_planetary_spin/particle.norm_spin_vector_2); // Msun.AU.day-1
        }
    }
}

pub fn calculate_torque_induced_by_rotational_flattening(rotational_flattening_host_particle: &Particle, particle: &Particle, central_body:bool) -> Axes {
    let mut reference_spin = rotational_flattening_host_particle.spin.clone();
    let orthogonal_component_of_the_force_induced_by_rotation: f64;

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

    Axes{x: torque_induced_by_rotational_flattening_x, y: torque_induced_by_rotational_flattening_y, z: torque_induced_by_rotational_flattening_z}

}

pub fn calculate_acceleration_induced_by_rotational_flattering(rotational_flattening_host_particle: &Particle, particle: &Particle) -> Axes {
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

    Axes{x: total_force_induced_by_rotation_x, y: total_force_induced_by_rotation_y, z: total_force_induced_by_rotation_z}
}
