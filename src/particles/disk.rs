use super::super::constants::{BOLTZMANN_CONSTANT, MASS_HYDROGEN_ATOM};
use super::{Particle};
use super::{Axes};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiskProperties {
    pub inner_edge_distance: f64,
    pub outer_edge_distance: f64,
    pub lifetime: f64,
    pub alpha: f64,
    pub surface_density_normalization: f64,
    pub mean_molecular_weight: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum Disk {
    Host(DiskProperties),
    Interaction(bool),
    None,
}

impl DiskProperties {
    pub fn calculate_migration_timescale(&self, time: f64, host_particle_mass: f64, host_particle_mass_g: f64, particle_mass: f64, distance: f64, semi_major_axis: f64) -> f64 {
        let gas_radial_velocity = self.calculate_gas_radial_velocity(host_particle_mass, host_particle_mass_g, distance, semi_major_axis);
        let disk_surface_density = self.calculate_disk_surface_density(time, distance);

        let x = 1.0_f64;
        let a_dot = gas_radial_velocity * x.min(2.0 * disk_surface_density * semi_major_axis.powi(2) / particle_mass);
        let timescale = -semi_major_axis / a_dot;
        timescale
    }

    fn calculate_temperature_disk(&self, host_particle_mass: f64, distance: f64) -> f64 {
        // distance in AU, host_particle_mass in M_SUN
        let t_ref = 280.0; //K
        let tdisk = t_ref * distance.powf(-0.5) * host_particle_mass; // K
        tdisk
    }


    fn calculate_disk_surface_density(&self, time: f64, distance: f64) -> f64{
        let initial_disk_surface_density = self.surface_density_normalization * distance.powi(-1)
            * (-1.0 * distance/self.outer_edge_distance).exp() * (1.0 - (self.inner_edge_distance/distance).sqrt()); // Unit of surface_density_normalization

        let mut disk_surface_density = 0.0;
        if distance > self.inner_edge_distance {
            disk_surface_density = initial_disk_surface_density * (-time/self.lifetime).exp();
        }
        disk_surface_density
    }

    fn calculate_keplerian_velocity(&self, host_particle_mass_g: f64, distance: f64) -> f64 {
        let keplerian_frequency = (host_particle_mass_g / (distance).powi(3)).sqrt(); // days^-1
        keplerian_frequency
    }

    fn calculate_disk_viscosity(&self, host_particle_mass: f64, host_particle_mass_g: f64, distance: f64) -> f64{
        let keplerian_frequency = self.calculate_keplerian_velocity(host_particle_mass_g, distance);
        let disk_temperature = self.calculate_temperature_disk(host_particle_mass, distance);
        let speed_of_sound_squared = BOLTZMANN_CONSTANT * disk_temperature / (self.mean_molecular_weight * MASS_HYDROGEN_ATOM);
        let viscosity = self.alpha * speed_of_sound_squared / keplerian_frequency; // AU^2.days^-1
        viscosity
    }

    fn calculate_gas_radial_velocity(&self, host_particle_mass: f64, host_particle_mass_g: f64, distance: f64, semi_major_axis: f64) -> f64{
        let viscosity = self.calculate_disk_viscosity(host_particle_mass, host_particle_mass_g, distance);
        let gas_radial_velocity = -3.0 * viscosity / (2.0 * semi_major_axis); // AU.days^-1
        gas_radial_velocity
    }

}

pub fn calculate_disk_interaction_acceleration(current_time: f64, disk_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let Disk::Host(disk) = disk_host_particle.disk {
        let factor2 = 1. / disk_host_particle.mass;
        let mut sum_total_disk_interaction_force = Axes{x:0., y:0., z:0.};

        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let Disk::Interaction(true) = particle.disk {
                //// Compute norm of the velocity and distance between the body in the
                //// center of the disk and the current body
                let norm_velocity_vector;
                let distance;
                let position;
                let velocity;
                //if self.disk_host_particle_index != self.tidal_host_particle_index {
                    // Norm of the velocity
                    let norm_velocity_vector_2 = (particle.velocity.x - disk_host_particle.velocity.x).powi(2) 
                                                + (particle.velocity.y - disk_host_particle.velocity.y).powi(2)
                                                + (particle.velocity.z - disk_host_particle.velocity.z).powi(2);
                    norm_velocity_vector = norm_velocity_vector_2.sqrt();

                    // (distance to star)^2
                    let distance_2 = (particle.position.x - disk_host_particle.position.x).powi(2) 
                                        + (particle.position.y - disk_host_particle.position.y).powi(2)
                                        + (particle.position.z - disk_host_particle.position.z).powi(2);
                    distance = distance_2.sqrt();
                    position = Axes{
                                    x: particle.inertial_position.x - disk_host_particle.inertial_position.x,
                                    y: particle.inertial_position.y - disk_host_particle.inertial_position.y,
                                    z: particle.inertial_position.z - disk_host_particle.inertial_position.z
                                    };
                    velocity = Axes{
                                    x: particle.inertial_velocity.x - disk_host_particle.inertial_velocity.x,
                                    y: particle.inertial_velocity.y - disk_host_particle.inertial_velocity.y,
                                    z: particle.inertial_velocity.z - disk_host_particle.inertial_velocity.z
                                    };
                //} else {
                    //norm_velocity_vector = particle.norm_velocity_vector;
                    //distance = particle.distance;
                    //position = particle.position;
                    //velocity = particle.velocity;
                //}

                // Semi-major axis
                let sma = (particle.mass_g+disk_host_particle.mass_g) * distance / (2.0 * (particle.mass_g+disk_host_particle.mass_g) -  distance * norm_velocity_vector*norm_velocity_vector);

                particle.migration_timescale = disk.calculate_migration_timescale(current_time, disk_host_particle.mass, disk_host_particle.mass_g, particle.mass, distance, sma);

                let factor1 = 1. / particle.mass;
                // From Alibert et al. 2013 (https://ui.adsabs.harvard.edu/abs/2013A&A...558A.109A)
                let factor_migration = -particle.mass / particle.migration_timescale;

                let eccentricity_damping_timescale = 0.1 * particle.migration_timescale;
                let factor_damping_eccentricity = -particle.mass * 2.0 / eccentricity_damping_timescale;

                let inclination_damping_timescale = eccentricity_damping_timescale;
                let factor_damping_inclination = -particle.mass * 2.0 /inclination_damping_timescale;

                // Force responsible for migration: Eq 8 from Alibert et al. 2013
                let disk_interaction_force_x = factor_migration * velocity.x;
                let disk_interaction_force_y = factor_migration * velocity.y;
                let disk_interaction_force_z = factor_migration * velocity.z;

                // Force responsible for eccentricity damping: Eq 9
                let scalar_product_velocity_radius_over_radius_squared = 1.0/(distance * distance)
                    * (position.x * velocity.x 
                       + position.y * velocity.y
                       + position.z * velocity.z);
                let disk_interaction_eccentricity_damping_force_x = factor_damping_eccentricity * scalar_product_velocity_radius_over_radius_squared * position.x;
                let disk_interaction_eccentricity_damping_force_y = factor_damping_eccentricity * scalar_product_velocity_radius_over_radius_squared * position.y;
                let disk_interaction_eccentricity_damping_force_z = factor_damping_eccentricity * scalar_product_velocity_radius_over_radius_squared * position.z;

                // Force responsible for inclination damping: Eq 10
                let disk_interaction_inclination_damping_force_x = 0.0;
                let disk_interaction_inclination_damping_force_y = 0.0;
                let disk_interaction_inclination_damping_force_z = factor_damping_inclination * velocity.z;

                // Total
                let total_disk_interaction_force_x = disk_interaction_force_x + disk_interaction_eccentricity_damping_force_x + disk_interaction_inclination_damping_force_x;
                let total_disk_interaction_force_y = disk_interaction_force_y + disk_interaction_eccentricity_damping_force_y + disk_interaction_inclination_damping_force_y;
                let total_disk_interaction_force_z = disk_interaction_force_z + disk_interaction_eccentricity_damping_force_z + disk_interaction_inclination_damping_force_z;


                sum_total_disk_interaction_force.x += total_disk_interaction_force_x;
                sum_total_disk_interaction_force.y += total_disk_interaction_force_y;
                sum_total_disk_interaction_force.z += total_disk_interaction_force_z;

                // - As in Equation 19 from Bolmont et al. 2015 (first term) 
                particle.disk_interaction_acceleration.x = factor1 * total_disk_interaction_force_x; 
                particle.disk_interaction_acceleration.y = factor1 * total_disk_interaction_force_y;
                particle.disk_interaction_acceleration.z = factor1 * total_disk_interaction_force_z;
            }
        
            // Instead of the previous code, keep star disk_interaction acceleration separated:
            disk_host_particle.disk_interaction_acceleration.x = -1.0 * factor2 * sum_total_disk_interaction_force.x;
            disk_host_particle.disk_interaction_acceleration.y = -1.0 * factor2 * sum_total_disk_interaction_force.y;
            disk_host_particle.disk_interaction_acceleration.z = -1.0 * factor2 * sum_total_disk_interaction_force.z;
        }
    }
}
