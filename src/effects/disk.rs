use super::super::constants::{BOLTZMANN_CONSTANT, MASS_HYDROGEN_ATOM};
use super::super::{Particle};
use super::super::{Axes};

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
pub struct DiskParticleInternalParameters {
    pub distance: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub migration_timescale: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiskParticleOutputParameters {
    pub acceleration: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiskParticleParameters {
    pub internal: DiskParticleInternalParameters,
    pub output: DiskParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiskParticleCoordinates {
    // Positions/velocities in a heliocentric frame 
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum DiskEffect {
    CentralBody(DiskProperties),
    OrbitingBody,
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Disk {
    pub effect: DiskEffect,
    pub parameters: DiskParticleParameters,
    pub coordinates: DiskParticleCoordinates,
}

impl Disk {
    pub fn new(effect: DiskEffect) -> Disk {
        Disk {
            effect: effect,
            parameters: DiskParticleParameters {
                internal: DiskParticleInternalParameters {
                    distance: 0.,
                    norm_velocity_vector: 0.,
                    norm_velocity_vector_2: 0.,
                    migration_timescale: 0.,
                },
                output: DiskParticleOutputParameters {
                    acceleration: Axes{x: 0., y: 0., z: 0.},
                },
            },
            coordinates: DiskParticleCoordinates {
                position: Axes{x: 0., y: 0., z: 0.},
                velocity: Axes{x: 0., y: 0., z: 0.},
            },
        }
    }
}



pub fn initialize(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let DiskEffect::CentralBody(_disk) = host_particle.disk.effect {
        host_particle.disk.parameters.output.acceleration.x = 0.;
        host_particle.disk.parameters.output.acceleration.y = 0.;
        host_particle.disk.parameters.output.acceleration.z = 0.;
        //host_particle.disk.coordinates.radial_velocity = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let DiskEffect::OrbitingBody = particle.disk.effect {
                particle.disk.parameters.output.acceleration.x = 0.;
                particle.disk.parameters.output.acceleration.y = 0.;
                particle.disk.parameters.output.acceleration.z = 0.;
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let DiskEffect::CentralBody(_disk) = host_particle.disk.effect {
        // Inertial to Heliocentric positions/velocities
        host_particle.disk.coordinates.position.x = 0.;
        host_particle.disk.coordinates.position.y = 0.;
        host_particle.disk.coordinates.position.z = 0.;
        host_particle.disk.coordinates.velocity.x = 0.;
        host_particle.disk.coordinates.velocity.y = 0.;
        host_particle.disk.coordinates.velocity.z = 0.;
        host_particle.disk.parameters.internal.distance = 0.;
        host_particle.disk.parameters.internal.norm_velocity_vector = 0.;
        host_particle.disk.parameters.internal.norm_velocity_vector_2 = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let DiskEffect::OrbitingBody = particle.disk.effect {
                particle.disk.coordinates.position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                particle.disk.coordinates.position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                particle.disk.coordinates.position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                particle.disk.coordinates.velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.disk.coordinates.velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.disk.coordinates.velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.disk.parameters.internal.distance = (particle.disk.coordinates.position.x.powi(2) 
                                                                    + particle.disk.coordinates.position.y.powi(2)
                                                                    + particle.disk.coordinates.position.z.powi(2)).sqrt();
                particle.disk.parameters.internal.norm_velocity_vector_2 = (particle.disk.coordinates.velocity.x 
                                                                                   - host_particle.disk.coordinates.velocity.x).powi(2) 
                                                                                + (particle.disk.coordinates.velocity.y 
                                                                                   - host_particle.disk.coordinates.velocity.y).powi(2)
                                                                                + (particle.disk.coordinates.velocity.z 
                                                                                   - host_particle.disk.coordinates.velocity.z).powi(2);
                particle.disk.parameters.internal.norm_velocity_vector = particle.disk.parameters.internal.norm_velocity_vector_2.sqrt();
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let DiskEffect::CentralBody(_disk) = host_particle.disk.effect {
        host_particle.disk.coordinates.position = host_particle.heliocentric_position;
        host_particle.disk.coordinates.velocity = host_particle.heliocentric_velocity;
        host_particle.disk.parameters.internal.distance = host_particle.heliocentric_distance;
        host_particle.disk.parameters.internal.norm_velocity_vector = host_particle.heliocentric_norm_velocity_vector;
        host_particle.disk.parameters.internal.norm_velocity_vector_2 = host_particle.heliocentric_norm_velocity_vector_2;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let DiskEffect::OrbitingBody = particle.disk.effect {
                particle.disk.coordinates.position = particle.heliocentric_position;
                particle.disk.coordinates.velocity = particle.heliocentric_velocity;
                particle.disk.parameters.internal.distance = particle.heliocentric_distance;
                particle.disk.parameters.internal.norm_velocity_vector = particle.heliocentric_norm_velocity_vector;
                particle.disk.parameters.internal.norm_velocity_vector_2 = particle.heliocentric_norm_velocity_vector_2;
            }
        }
    }
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
    if let DiskEffect::CentralBody(disk) = disk_host_particle.disk.effect {
        let factor2 = 1. / disk_host_particle.mass;
        let mut sum_total_disk_interaction_force = Axes{x:0., y:0., z:0.};

        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let DiskEffect::OrbitingBody = particle.disk.effect {
                // Semi-major axis
                let sma = (particle.mass_g+disk_host_particle.mass_g) * particle.disk.parameters.internal.distance / (2.0 * (particle.mass_g+disk_host_particle.mass_g) -  particle.disk.parameters.internal.distance * particle.disk.parameters.internal.norm_velocity_vector_2);

                particle.disk.parameters.internal.migration_timescale = disk.calculate_migration_timescale(current_time, disk_host_particle.mass, disk_host_particle.mass_g, particle.mass, particle.disk.parameters.internal.distance, sma);

                let factor1 = 1. / particle.mass;
                // From Alibert et al. 2013 (https://ui.adsabs.harvard.edu/abs/2013A&A...558A.109A)
                let factor_migration = -particle.mass / particle.disk.parameters.internal.migration_timescale;

                let eccentricity_damping_timescale = 0.1 * particle.disk.parameters.internal.migration_timescale;
                let factor_damping_eccentricity = -particle.mass * 2.0 / eccentricity_damping_timescale;

                let inclination_damping_timescale = eccentricity_damping_timescale;
                let factor_damping_inclination = -particle.mass * 2.0 /inclination_damping_timescale;

                // Force responsible for migration: Eq 8 from Alibert et al. 2013
                let disk_interaction_force_x = factor_migration * particle.disk.coordinates.velocity.x;
                let disk_interaction_force_y = factor_migration * particle.disk.coordinates.velocity.y;
                let disk_interaction_force_z = factor_migration * particle.disk.coordinates.velocity.z;

                // Force responsible for eccentricity damping: Eq 9
                let scalar_product_velocity_radius_over_radius_squared = 1.0/(particle.disk.parameters.internal.distance * particle.disk.parameters.internal.distance)
                    * (particle.disk.coordinates.position.x * particle.disk.coordinates.velocity.x 
                       + particle.disk.coordinates.position.y * particle.disk.coordinates.velocity.y
                       + particle.disk.coordinates.position.z * particle.disk.coordinates.velocity.z);
                let disk_interaction_eccentricity_damping_force_x = factor_damping_eccentricity * scalar_product_velocity_radius_over_radius_squared * particle.disk.coordinates.position.x;
                let disk_interaction_eccentricity_damping_force_y = factor_damping_eccentricity * scalar_product_velocity_radius_over_radius_squared * particle.disk.coordinates.position.y;
                let disk_interaction_eccentricity_damping_force_z = factor_damping_eccentricity * scalar_product_velocity_radius_over_radius_squared * particle.disk.coordinates.position.z;

                // Force responsible for inclination damping: Eq 10
                let disk_interaction_inclination_damping_force_x = 0.0;
                let disk_interaction_inclination_damping_force_y = 0.0;
                let disk_interaction_inclination_damping_force_z = factor_damping_inclination * particle.disk.coordinates.velocity.z;

                // Total
                let total_disk_interaction_force_x = disk_interaction_force_x + disk_interaction_eccentricity_damping_force_x + disk_interaction_inclination_damping_force_x;
                let total_disk_interaction_force_y = disk_interaction_force_y + disk_interaction_eccentricity_damping_force_y + disk_interaction_inclination_damping_force_y;
                let total_disk_interaction_force_z = disk_interaction_force_z + disk_interaction_eccentricity_damping_force_z + disk_interaction_inclination_damping_force_z;


                sum_total_disk_interaction_force.x += total_disk_interaction_force_x;
                sum_total_disk_interaction_force.y += total_disk_interaction_force_y;
                sum_total_disk_interaction_force.z += total_disk_interaction_force_z;

                // - As in Equation 19 from Bolmont et al. 2015 (first term) 
                particle.disk.parameters.output.acceleration.x = factor1 * total_disk_interaction_force_x; 
                particle.disk.parameters.output.acceleration.y = factor1 * total_disk_interaction_force_y;
                particle.disk.parameters.output.acceleration.z = factor1 * total_disk_interaction_force_z;
            }
        
            // Instead of the previous code, keep star disk_interaction acceleration separated:
            disk_host_particle.disk.parameters.output.acceleration.x = -1.0 * factor2 * sum_total_disk_interaction_force.x;
            disk_host_particle.disk.parameters.output.acceleration.y = -1.0 * factor2 * sum_total_disk_interaction_force.y;
            disk_host_particle.disk.parameters.output.acceleration.z = -1.0 * factor2 * sum_total_disk_interaction_force.z;
        }
    }
}
