use super::{Particle};
use super::{Axes};

pub fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let mut flattening_host_particle = flattening_host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let central_body = true;
    calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(central_body, &mut flattening_host_particle, &mut particles, &mut more_particles);
    calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(!central_body, &mut flattening_host_particle, &mut particles, &mut more_particles);
}

fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(central_body:bool, flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    // Calculation of the norm square of the spin for the planet
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if central_body {
            // - Star Equation 16 from Bolmont et al. 2015
            particle.factor_for_the_force_induced_by_star_rotation = particle.mass * flattening_host_particle.fluid_love_number * flattening_host_particle.norm_spin_vector_2 * flattening_host_particle.radius.powi(5) / 6.; // Msun.AU^5.day-2
            // - Second part of Equation 15 from Bolmont et al. 2015
            particle.orthogonal_component_of_the_force_induced_by_star_rotation = -6. * particle.factor_for_the_force_induced_by_star_rotation * particle.scalar_product_of_vector_position_with_stellar_spin / (flattening_host_particle.norm_spin_vector_2 * particle.distance.powi(5));
        } else {
            // - Planet Equation 16 from Bolmont et al. 2015
            particle.factor_for_the_force_induced_by_planet_rotation = flattening_host_particle.mass * particle.fluid_love_number * particle.norm_spin_vector_2 * particle.radius.powi(5) / 6.; // Msun.AU^5.day-2
            // - Second part of Equation 15 from Bolmont et al. 2015
            particle.orthogonal_component_of_the_force_induced_by_planet_rotation = -6. * particle.factor_for_the_force_induced_by_planet_rotation * particle.scalar_product_of_vector_position_with_planetary_spin / (particle.norm_spin_vector_2 * particle.distance.powi(5));
        }
    }
}

pub fn calculate_radial_component_of_the_force_induced_by_rotational_flattening(flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // - First part of Equation 15 from Bolmont et al. 2015
        particle.radial_component_of_the_force_induced_by_rotation = -3./particle.distance.powi(5) * (particle.factor_for_the_force_induced_by_planet_rotation + particle.factor_for_the_force_induced_by_star_rotation)
            + 15./particle.distance.powi(7) * (particle.factor_for_the_force_induced_by_star_rotation * particle.scalar_product_of_vector_position_with_stellar_spin * particle.scalar_product_of_vector_position_with_stellar_spin/flattening_host_particle.norm_spin_vector_2
                + particle.factor_for_the_force_induced_by_planet_rotation * particle.scalar_product_of_vector_position_with_planetary_spin * particle.scalar_product_of_vector_position_with_planetary_spin/particle.norm_spin_vector_2); // Msun.AU.day-1
    }
}

pub fn calculate_torque_induced_by_rotational_flattening(flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool) {
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};
    let mut reference_spin = flattening_host_particle.spin.clone();
    let mut orthogonal_component_of_the_force_induced_by_rotation: f64;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if !central_body {
            reference_spin = particle.spin.clone();
            orthogonal_component_of_the_force_induced_by_rotation = particle.orthogonal_component_of_the_force_induced_by_planet_rotation;
        } else {
            orthogonal_component_of_the_force_induced_by_rotation = particle.orthogonal_component_of_the_force_induced_by_star_rotation;
        }
        
        //// Torque calculation due to rotational flattening
        // - Equation 17-18 from Bolmont et al. 2015
        let torque_induced_by_rotational_flattening_x: f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.position.y * reference_spin.z - particle.position.z * reference_spin.y);
        let torque_induced_by_rotational_flattening_y :f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.position.z * reference_spin.x - particle.position.x * reference_spin.z);
        let torque_induced_by_rotational_flattening_z: f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.position.x * reference_spin.y - particle.position.y * reference_spin.x);

        let factor = -1.0;
        // - Equation 25 from Bolmont et al. 2015
        if central_body {
            // Integration of the spin (total torque rot):
            dangular_momentum_dt.x += factor * torque_induced_by_rotational_flattening_x;
            dangular_momentum_dt.y += factor * torque_induced_by_rotational_flattening_y;
            dangular_momentum_dt.z += factor * torque_induced_by_rotational_flattening_z;
        } else {
            particle.dangular_momentum_dt_induced_by_rotational_flattening.x = factor * torque_induced_by_rotational_flattening_x;
            particle.dangular_momentum_dt_induced_by_rotational_flattening.y = factor * torque_induced_by_rotational_flattening_y;
            particle.dangular_momentum_dt_induced_by_rotational_flattening.z = factor * torque_induced_by_rotational_flattening_z;
        }
    }

    if central_body {
        // - Equation 25 from Bolmont et al. 2015
        flattening_host_particle.dangular_momentum_dt_induced_by_rotational_flattening.x = dangular_momentum_dt.x;
        flattening_host_particle.dangular_momentum_dt_induced_by_rotational_flattening.y = dangular_momentum_dt.y;
        flattening_host_particle.dangular_momentum_dt_induced_by_rotational_flattening.z = dangular_momentum_dt.z;
    }
}

pub fn calculate_acceleration_induced_by_rotational_flattering(flattening_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor2 = 1. / flattening_host_particle.mass;

    // Calculation of the norm square of the spin for the star
    let mut sum_total_force_induced_by_rotation = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        let factor1 = 1. / particle.mass;

        // - Equation 15 from Bolmont et al. 2015
        let total_force_induced_by_rotation_x = particle.radial_component_of_the_force_induced_by_rotation * particle.position.x
            + particle.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.x
            + particle.orthogonal_component_of_the_force_induced_by_star_rotation * flattening_host_particle.spin.x;
        let total_force_induced_by_rotation_y = particle.radial_component_of_the_force_induced_by_rotation * particle.position.y
            + particle.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.y
            + particle.orthogonal_component_of_the_force_induced_by_star_rotation * flattening_host_particle.spin.y;
        let total_force_induced_by_rotation_z = particle.radial_component_of_the_force_induced_by_rotation * particle.position.z
            + particle.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.z
            + particle.orthogonal_component_of_the_force_induced_by_star_rotation * flattening_host_particle.spin.z;
        //println!("Rot force = {:e} {:e} {:e}", total_force_induced_by_rotation_x, total_force_induced_by_rotation_y, total_force_induced_by_rotation_z);

        sum_total_force_induced_by_rotation.x += total_force_induced_by_rotation_x;
        sum_total_force_induced_by_rotation.y += total_force_induced_by_rotation_y;
        sum_total_force_induced_by_rotation.z += total_force_induced_by_rotation_z;
          
        // - Equation 19 from Bolmont et al. 2015 (first term)
        particle.acceleration_induced_by_rotational_flattering.x = factor1 * total_force_induced_by_rotation_x; 
        particle.acceleration_induced_by_rotational_flattering.y = factor1 * total_force_induced_by_rotation_y;
        particle.acceleration_induced_by_rotational_flattering.z = factor1 * total_force_induced_by_rotation_z;
    }
    
    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.acceleration_induced_by_rotational_flattering.x += factor2 * sum_total_force_induced_by_rotation.x;
        //particle.acceleration_induced_by_rotational_flattering.y += factor2 * sum_total_force_induced_by_rotation.y;
        //particle.acceleration_induced_by_rotational_flattering.z += factor2 * sum_total_force_induced_by_rotation.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    flattening_host_particle.acceleration_induced_by_rotational_flattering.x = -1.0 * factor2 * sum_total_force_induced_by_rotation.x;
    flattening_host_particle.acceleration_induced_by_rotational_flattering.y = -1.0 * factor2 * sum_total_force_induced_by_rotation.y;
    flattening_host_particle.acceleration_induced_by_rotational_flattering.z = -1.0 * factor2 * sum_total_force_induced_by_rotation.z;
}

