use super::super::constants::{BOLTZMANN_CONSTANT, MASS_HYDROGEN_ATOM, G};

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
    Properties(DiskProperties),
    None,
}

impl DiskProperties {
    pub fn calculate_migration_timescale(&self, time: f64, planet_distance: f64, planet_semi_major_axis: f64, planet_mass: f64, star_mass: f64) -> f64 {
        let gas_radial_velocity = self.calculate_gas_radial_velocity(planet_distance, star_mass, planet_semi_major_axis);
        let disk_surface_density = self.calculate_disk_surface_density(time, planet_distance);

        let x = 1.0_f64;
        let a_dot = gas_radial_velocity * x.min(2.0 * disk_surface_density * planet_semi_major_axis.powi(2) / planet_mass);
        let timescale = -planet_semi_major_axis / a_dot;
        timescale
    }

    fn calculate_temperature_disk(&self, planet_distance: f64, star_mass: f64) -> f64 {
        // planet_distance in AU, star_mass in M_SUN
        let t_ref = 280.0; //K
        let tdisk = t_ref * planet_distance.powf(-0.5) * star_mass; // K
        tdisk
    }


    fn calculate_disk_surface_density(&self, time: f64, planet_distance: f64) -> f64{
        let initial_disk_surface_density = self.surface_density_normalization * planet_distance.powi(-1)
            * (-1.0 * planet_distance/self.outer_edge_distance).exp() * (1.0 - (self.inner_edge_distance/planet_distance).sqrt()); // Unit of surface_density_normalization

        let mut disk_surface_density = 0.0;
        if planet_distance > self.inner_edge_distance {
            disk_surface_density = initial_disk_surface_density * (-time/self.lifetime).exp();
        }
        disk_surface_density
    }

    fn calculate_keplerian_velocity(&self, planet_distance: f64, star_mass_g: f64) -> f64 {
        let keplerian_frequency = (star_mass_g / (planet_distance).powi(3)).sqrt(); // days^-1
        keplerian_frequency
    }

    fn calculate_disk_viscosity(&self, planet_distance: f64, star_mass: f64) -> f64{
        let star_mass_g = G * star_mass;
        let keplerian_frequency = self.calculate_keplerian_velocity(planet_distance, star_mass_g);
        let disk_temperature = self.calculate_temperature_disk(planet_distance, star_mass);
        let speed_of_sound_squared = BOLTZMANN_CONSTANT * disk_temperature / (self.mean_molecular_weight * MASS_HYDROGEN_ATOM);
        let viscosity = self.alpha * speed_of_sound_squared / keplerian_frequency; // AU^2.days^-1
        viscosity
    }

    fn calculate_gas_radial_velocity(&self, planet_distance: f64, star_mass: f64, planet_semi_major_axis: f64) -> f64{
        let viscosity = self.calculate_disk_viscosity(planet_distance, star_mass);
        let gas_radial_velocity = -3.0 * viscosity / (2.0 * planet_semi_major_axis); // AU.days^-1
        gas_radial_velocity
    }

}
