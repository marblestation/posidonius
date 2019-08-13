use std::collections::HashMap;
use super::super::tools;
use super::super::constants::{K2};
use super::{Particle};
use super::{Axes};
use super::{EvolutionType};


pub fn calculate_planet_dependent_dissipation_factors(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    match tidal_host_particle.evolution_type {
        EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
            let star_norm_spin_vector = tidal_host_particle.norm_spin_vector_2.sqrt();
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                //
                //// Excitation frequency needed by the model based on the
                // instantaneous frequency (using positions, velocities and spins)
                //let frequency = (particle.velocity.x - tidal_host_particle.spin.y*particle.position.z + tidal_host_particle.spin.z*particle.position.y).powi(2)
                            //+ (particle.velocity.y - tidal_host_particle.spin.z*particle.position.x + tidal_host_particle.spin.x*particle.position.z).powi(2)
                            //+ (particle.velocity.z - tidal_host_particle.spin.x*particle.position.y + tidal_host_particle.spin.y*particle.position.x).powi(2);
                //let inverse_of_half_the_excitation_frequency = particle.distance / frequency;
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
                let (perihelion_distance, eccentricity) = tools::calculate_perihelion_distance_and_eccentricity(gm, particle.position, particle.velocity);
                let mean_motion = gm.sqrt() * (perihelion_distance/(1.0 - eccentricity)).powf(-1.5);
                let half_the_excitation_frequency = (star_norm_spin_vector - mean_motion).abs();
                let inverse_of_half_the_excitation_frequency = 1./half_the_excitation_frequency;

                let planet_dependent_dissipation_factor = tidal_host_particle.dissipation_factor_scale * 2.0 * K2
                    * tidal_host_particle.lag_angle * inverse_of_half_the_excitation_frequency / (3.0*tidal_host_particle.radius.powi(5));

                star_planet_dependent_dissipation_factors.insert(particle.id.clone(), planet_dependent_dissipation_factor);
                //println!("Insert {} in {}", planet_dependent_dissipation_factor, particle.id);
            }
            //panic!("Please, contact Posidonius authors before using BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true) evolutionary models. They may not be ready yet for scientific explotation.")
        },
        _ => {},
    }

}

pub fn planet_dependent_dissipation_factor(star_planet_dependent_dissipation_factors: &HashMap<usize, f64>,  id: &usize, evolution_type: EvolutionType, scaled_dissipation_factor: f64) -> f64 {
    match evolution_type {
        EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
            match star_planet_dependent_dissipation_factors.get(id) {
                Some(&value) => value,
                _ => scaled_dissipation_factor // This should not happen
            }
        },
        _ => scaled_dissipation_factor,
    }
}

pub fn calculate_scalar_product_of_vector_position_with_spin (tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        particle.scalar_product_of_vector_position_with_stellar_spin = particle.position.x * tidal_host_particle.spin.x 
                        + particle.position.y * tidal_host_particle.spin.y
                        + particle.position.z * tidal_host_particle.spin.z;
    
        particle.scalar_product_of_vector_position_with_planetary_spin = particle.position.x * particle.spin.x 
                        + particle.position.y * particle.spin.y
                        + particle.position.z * particle.spin.z;
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

        if !central_body {
            reference_spin = particle.spin.clone();
            reference_rscalspin = particle.scalar_product_of_vector_position_with_planetary_spin;
            orthogonal_component_of_the_tidal_force = particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
        } else {
            reference_rscalspin = particle.scalar_product_of_vector_position_with_stellar_spin;
            orthogonal_component_of_the_tidal_force = particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
        }

        // distance to star
        let distance = particle.distance;

        //// Torque calculation (star)
        // - Equation 8-9 from Bolmont et al. 2015
        let torque_due_to_tides_x: f64 = orthogonal_component_of_the_tidal_force 
                        * (distance * reference_spin.x - reference_rscalspin*particle.position.x/distance - 1.0/distance
                        * (particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y) );

        let torque_due_to_tides_y :f64 = orthogonal_component_of_the_tidal_force 
                        * (distance * reference_spin.y - reference_rscalspin*particle.position.y/distance - 1.0/distance
                        * (particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z) );

        let torque_due_to_tidez_z: f64 = orthogonal_component_of_the_tidal_force 
                        * (distance * reference_spin.z - reference_rscalspin*particle.position.z/distance - 1.0/distance
                        * (particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x) );

        let factor = -1.0;
        if central_body {
            // Integration of the spin (total torque tides):
            dangular_momentum_dt.x += factor * torque_due_to_tides_x;
            dangular_momentum_dt.y += factor * torque_due_to_tides_y;
            dangular_momentum_dt.z += factor * torque_due_to_tidez_z;
        } else {
            particle.dangular_momentum_dt_due_to_tides.x = factor * torque_due_to_tides_x;
            particle.dangular_momentum_dt_due_to_tides.y = factor * torque_due_to_tides_y;
            particle.dangular_momentum_dt_due_to_tides.z = factor * torque_due_to_tidez_z;
        }
    }

    if central_body {
        // - Equation 25 from Bolmont et al. 2015
        tidal_host_particle.dangular_momentum_dt_due_to_tides.x = dangular_momentum_dt.x;
        tidal_host_particle.dangular_momentum_dt_due_to_tides.y = dangular_momentum_dt.y;
        tidal_host_particle.dangular_momentum_dt_due_to_tides.z = dangular_momentum_dt.z;
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

        //// Only calculate tides if planet is not in disk
        //if particle.disk_interaction_time == 0.0 {

        // (distance to star)^7
        let distance_7 = particle.distance.powi(7);

        //// Tidal force calculation (star) :: Only orthogonal component is needed
        if central_body {
            // - Third line of Equation 5 from Bolmont et al. 2015
            //   This expression has R**10 (instead of R**5 in Eq. 5) 
            //   because it uses sigma (i.e., scaled_dissipation_factor) 
            //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
            //   there is a R**5 factor as shown in Equation 28)
            //   - k2 is love number
            let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(&star_planet_dependent_dissipation_factors, &tidal_host_particle.id, tidal_host_particle.evolution_type, tidal_host_particle.scaled_dissipation_factor);
            particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass.powi(2))
                                            * (tidal_host_particle.radius.powi(10)) 
                                            * star_scaled_dissipation_factor / distance_7;
        } else {
            // - Second line of Equation 5 from Bolmont et al. 2015
            //   This expression has R**10 (instead of R**5 in Eq. 5) 
            //   because it uses sigma (i.e., scaled_dissipation_factor) 
            //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
            //   there is a R**5 factor as shown in Equation 28)
            //   - k2 is love number
            particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (tidal_host_particle.mass.powi(2))
                                            * (particle.radius.powi(10))
                                            * particle.scaled_dissipation_factor / distance_7;

            // SBC
            //println!("> {:e} {:e} {:e} {:e}", tidal_host_particle.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
            //println!("> {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z);
        }
        //} else {
            //if central_body {
                //particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 0.0
            //} else{
                //particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0
            //}
        //}
    }
}

pub fn calculate_radial_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    let star_mass_2 = tidal_host_particle.mass * tidal_host_particle.mass;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        let planet_mass_2 = particle.mass * particle.mass;
        // Conservative part of the radial tidal force
        let radial_component_of_the_tidal_force_conservative_part = -3.0 * K2 / particle.distance.powi(7)
                    * (planet_mass_2 * tidal_host_particle.radius.powi(5) * tidal_host_particle.love_number 
                    + star_mass_2 * particle.radius.powi(5) * particle.love_number);

        // Dissipative part of the radial tidal force:
        let factor1 = -13.5 * particle.radial_velocity / particle.distance.powi(8);
        let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(&star_planet_dependent_dissipation_factors, &particle.id, tidal_host_particle.evolution_type, tidal_host_particle.scaled_dissipation_factor);
        let term1 = planet_mass_2
                    * tidal_host_particle.radius.powi(10)
                    * star_scaled_dissipation_factor;
        let term2 = star_mass_2
                    * particle.radius.powi(10)
                    * particle.scaled_dissipation_factor;
        // If we consider the star as a point mass (used for denergy_dt calculation):
        particle.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass = factor1 * term2;
        let radial_component_of_the_tidal_force_dissipative_part = particle.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor1 * term1;

        // Sum of the dissipative and conservative part of the radial force
        // - First line Equation 5 from Bolmont et al. 2015
        particle.radial_component_of_the_tidal_force = radial_component_of_the_tidal_force_conservative_part + radial_component_of_the_tidal_force_dissipative_part;
    }
}

//pub fn calculate_denergy_dt(particles: &mut [Particle], more_particles: &mut [Particle]) {
    //for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        //// - Equation 32 from Bolmont et al. 2015
        ////// Instantaneous energy loss dE/dt due to tides
        ////// in Msun.AU^2.day^(-3)
        ////radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
        //let factor2 = particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance;
        //particle.denergy_dt = -((1.0 / particle.distance * (particle.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.radial_velocity))
                    //* (particle.position.x*particle.velocity.x + particle.position.y*particle.velocity.y + particle.position.z*particle.velocity.z)
                    //+ factor2 
                    //* ((particle.spin.y*particle.position.z - particle.spin.z*particle.position.y - particle.velocity.x) * particle.velocity.x
                    //+ (particle.spin.z*particle.position.x - particle.spin.x*particle.position.z - particle.velocity.y) * particle.velocity.y
                    //+ (particle.spin.x*particle.position.y - particle.spin.y*particle.position.x - particle.velocity.z) * particle.velocity.z))
                    //- (particle.dangular_momentum_dt_due_to_tides.x*particle.spin.x + particle.dangular_momentum_dt_due_to_tides.y*particle.spin.y + particle.dangular_momentum_dt_due_to_tides.z*particle.spin.z);
    //}
//}


pub fn calculate_tidal_acceleration(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor2 = 1. / tidal_host_particle.mass;
    let mut sum_total_tidal_force = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        let factor1 = 1. / particle.mass;

        // - Equation 6 from Bolmont et al. 2015
        let factor3 = particle.radial_component_of_the_tidal_force
                        + (particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide) * particle.radial_velocity / particle.distance;
        let total_tidal_force_x = factor3 * particle.position.x / particle.distance
                                + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                    * (tidal_host_particle.spin.y * particle.position.z  - tidal_host_particle.spin.z * particle.position.y - particle.velocity.x)
                                + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                    * (particle.spin.y * particle.position.z  - particle.spin.z * particle.position.y - particle.velocity.x);
        let total_tidal_force_y = factor3 * particle.position.y / particle.distance
                                + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                    * (tidal_host_particle.spin.z * particle.position.x  - tidal_host_particle.spin.x * particle.position.z - particle.velocity.y)
                                + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                    * (particle.spin.z * particle.position.x  - particle.spin.x * particle.position.z - particle.velocity.y);
        let total_tidal_force_z = factor3 * particle.position.z / particle.distance
                                + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                    * (tidal_host_particle.spin.x * particle.position.y  - tidal_host_particle.spin.y * particle.position.x - particle.velocity.z)
                                + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                    * (particle.spin.x * particle.position.y  - particle.spin.y * particle.position.x - particle.velocity.z);
        //println!("factor3 {:e} {:e} {:e}", particle.radial_component_of_the_tidal_force, particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
        //println!("d {:e} vrad {:e}", particle.distance, particle.radial_velocity);
        //println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);

        sum_total_tidal_force.x += total_tidal_force_x;
        sum_total_tidal_force.y += total_tidal_force_y;
        sum_total_tidal_force.z += total_tidal_force_z;

        // - Equation 19 from Bolmont et al. 2015 (first term)
        particle.tidal_acceleration.x = factor1 * total_tidal_force_x; 
        particle.tidal_acceleration.y = factor1 * total_tidal_force_y;
        particle.tidal_acceleration.z = factor1 * total_tidal_force_z;
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.tidal_acceleration.x += factor2 * sum_total_tidal_force.x;
        //particle.tidal_acceleration.y += factor2 * sum_total_tidal_force.y;
        //particle.tidal_acceleration.z += factor2 * sum_total_tidal_force.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    tidal_host_particle.tidal_acceleration.x = -1.0 * factor2 * sum_total_tidal_force.x;
    tidal_host_particle.tidal_acceleration.y = -1.0 * factor2 * sum_total_tidal_force.y;
    tidal_host_particle.tidal_acceleration.z = -1.0 * factor2 * sum_total_tidal_force.z;
}


