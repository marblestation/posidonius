extern crate time;
use std::collections::HashMap;
use super::super::constants::{K2, G, R_SUN, SUN_DYN_FREQ, SPEED_OF_LIGHT_2, MAX_PARTICLES, MAX_DISTANCE_2, DBL_EPSILON_2};
use super::{Evolver, EvolutionType};
use super::{Particle};
use super::{Axes};

#[derive(Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub struct Universe {
    pub initial_time: f64,
    pub time_limit: f64,
    pub particles: [Particle; MAX_PARTICLES],
    pub particles_evolvers: Vec<Evolver>,
    pub n_particles: usize,
    pub evolving_particles_exist: bool,
    pub wind_effects_exist: bool,
    pub consider_tides: bool,
    pub consider_rotational_flattening: bool,
    pub consider_general_relativy: ConsiderGeneralRelativity,
    star_planet_dependent_dissipation_factors : HashMap<usize, f64>, // Central body specific
    temporary_copied_particle_positions: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particle_velocities: [Axes; MAX_PARTICLES], // For optimization purposes (TODO: Delete and adapt python package)
    temporary_copied_particles_masses: [f64; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_radiuses: [f64; MAX_PARTICLES], // For optimization purposes
}

#[derive(Debug, Copy, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub enum IgnoreGravityTerms {
    None,
    WHFastOne,
    WHFastTwo,
}

#[derive(Debug, Copy, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub enum ConsiderGeneralRelativity {
    None,
    Kidder1995, // MercuryT
    Anderson1975, // REBOUNDx gr
    Newhall1983, // REBOUNDx gr full
}

impl Universe {
    pub fn gravity_calculate_acceleration(&mut self, ignore_terms: IgnoreGravityTerms) {
        //let softening = 1e-12;
        // TODO: Try this pattern for gravity calculation
        //for (i, particle_a) in self.particles[..self.n_particles].iter().enumerate() {
            //if let Some((newtonian_acceleration_a, newtonian_accelerations_b)) = newtonian_accelerations[i..self.n_particles].split_first_mut() {
                //for (newtonian_acceleration_b, particle_b) in newtonian_accelerations_b[..self.n_particles-i-1].iter_mut().zip( self.particles[i+1..self.n_particles].iter()) {
                    //let dx = particle_a.inertial_position.x - particle_b.inertial_position.x;
                    //let dy = particle_a.inertial_position.y - particle_b.inertial_position.y;
                    //let dz = particle_a.inertial_position.z - particle_b.inertial_position.z;
                    //let r2 = dx.powi(2) + dy.powi(2) + dz.powi(2);
                    ////println!("r2 {:e} {:e} {:e} {:e}", r2, dx.powi(2), dy.powi(2), dz.powi(2));
                    //let r = r2.sqrt();
                    //let prefac = G/(r2*r);
                    ////println!("prefac {:e} = {:e}/({:e}*{:e})", prefac, G, r2, r);
                    //let prefac_mass_a = prefac*particle_a.mass;
                    //let prefac_mass_b = prefac*particle_b.mass;
                    ////println!("-- mass_a {:e} | mass_b {:e}", particle_a.mass, particle_b.mass);
                    ////println!("* berfore: {:e}", newtonian_acceleration_a.x);
                    //newtonian_acceleration_a.x -= prefac_mass_b*dx;
                    ////println!("* after: {:e}", newtonian_acceleration_a.x);
                    //newtonian_acceleration_a.y -= prefac_mass_b*dy;
                    //newtonian_acceleration_a.z -= prefac_mass_b*dz;
                    ////println!("+ berfore: {:e}", newtonian_acceleration_b.x);
                    //newtonian_acceleration_b.x += prefac_mass_a*dx;
                    ////println!("+ after: {:e}", newtonian_acceleration_b.x);
                    //newtonian_acceleration_b.y += prefac_mass_a*dy;
                    //newtonian_acceleration_b.z += prefac_mass_a*dz;
                //}
            //}
        //}

        for (i, particle) in self.particles[..self.n_particles].iter().enumerate() {
            self.temporary_copied_particle_positions[i].x = particle.inertial_position.x;
            self.temporary_copied_particle_positions[i].y = particle.inertial_position.y;
            self.temporary_copied_particle_positions[i].z = particle.inertial_position.z;
            self.temporary_copied_particles_masses[i] = particle.mass;
            self.temporary_copied_particles_radiuses[i] = particle.radius;
        }

        for (i, particle_a) in self.particles[..self.n_particles].iter_mut().enumerate() {
			particle_a.inertial_acceleration.x = 0.;
			particle_a.inertial_acceleration.y = 0.;
			particle_a.inertial_acceleration.z = 0.;
            for (j, ((particle_b_position, particle_b_radius), particle_b_mass)) in self.temporary_copied_particle_positions[..self.n_particles].iter()
                                                                .zip(self.temporary_copied_particles_radiuses[..self.n_particles].iter())
                                                                .zip(self.temporary_copied_particles_masses[..self.n_particles].iter()).enumerate() {
                if i == j {
                    continue;
                }

                //// Check collisions and ejections //////////////////////////////////
                let dx = particle_a.inertial_position.x - particle_b_position.x;
                let dy = particle_a.inertial_position.y - particle_b_position.y;
                let dz = particle_a.inertial_position.z - particle_b_position.z;
                let distance_2 = dx*dx + dy*dy + dz*dz;
                // Do not check twice the same pair of particles:
                if i < j {
                    // TODO: Optimize roche radius calculation and do it just once
                    // Faber et al, 2005; Pacynski, 1971
                    let roche_radius = (particle_a.radius/0.462)*(particle_a.mass/particle_b_mass).powf(-1./3.);
                    if distance_2 <= roche_radius.powi(2) {
                        println!("\n");
                        panic!("[PANIC {} UTC] Particle {} was destroyed by particle {} due to strong tides (close encounter)!", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), i, j);
                    }
                    // Check if particles are overlapping
                    if distance_2 <= (particle_a.radius + particle_b_radius).powi(2) {
                        println!("\n");
                        panic!("[PANIC {} UTC] Collision between particle {} and {}!", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), i, j);
                    }
                    if i == 0 && distance_2 > MAX_DISTANCE_2 {
                        println!("\n");
                        panic!("[PANIC {} UTC] Particle {} has been ejected!", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), j);
                    }
                }
                //////////////////////////////////////////////////////////////////////

                if ignore_terms == IgnoreGravityTerms::WHFastOne && ((i == 1 && j == 0) || (i == 0 && j == 1)) {
                    // For WHFast Jacobi coordinates, 
                    // ignore interaction between central body and first planet
                    continue
                }
                if ignore_terms == IgnoreGravityTerms::WHFastTwo && (i == 0 || j == 0) {
                    // For WHFast democratic-heliocentric & WHDS coordinates
                    // completely ignore central body
                    continue
                }

                // gravitational force equation for vector quantities
                // Fx = - G * ((m1 * m2) / r^3) * (x2 - x1)
                // Fy = - G * ((m1 * m2) / r^3) * (y2 - y1)
                // Fz = - G * ((m1 * m2) / r^3) * (z2 - z1)
                //
                // acceleration:
                // ax = Fx / m1 = - G * (m2 / r^3) * (x2 - x1)
                //let dx = particle_a.inertial_position.x - particle_b_position.x;
                //let dy = particle_a.inertial_position.y - particle_b_position.y;
                //let dz = particle_a.inertial_position.z - particle_b_position.z;
                //let r = (dx*dx + dy*dy + dz*dz).sqrt();
                //let r = (distance_2 + softening).sqrt();
                let distance = distance_2.sqrt();
                let prefact = -G/(distance_2*distance) * particle_b_mass;

                particle_a.inertial_acceleration.x += prefact * dx;
                particle_a.inertial_acceleration.y += prefact * dy;
                particle_a.inertial_acceleration.z += prefact * dz;
            }
        }

        //if ignore_terms != IgnoreGravityTerms::WHFastOne && ignore_terms != IgnoreGravityTerms::WHFastTwo {
            //// Tides require heliocentric point of reference, the star should continue in the zero point
            //// so we must compensate all the planets (but if WHFast is being used, it is
            //// automatically done by the integrator):
            //if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
                //for particle in particles.iter_mut() {
                    //particle.inertial_acceleration.x -= star.inertial_acceleration.x;
                    //particle.inertial_acceleration.y -= star.inertial_acceleration.y;
                    //particle.inertial_acceleration.z -= star.inertial_acceleration.z;
                //}
                //star.inertial_acceleration.x = 0.;
                //star.inertial_acceleration.y = 0.;
                //star.inertial_acceleration.z = 0.;
            //}
        //}
    }

    pub fn calculate_additional_effects(&mut self, current_time: f64, evolution: bool, dspin_dt: bool, accelerations: bool, ignored_gravity_terms: IgnoreGravityTerms) {
        if (evolution && self.evolving_particles_exist) || 
            ((dspin_dt || accelerations) && self.consider_rotational_flattening) {
            self.calculate_norm_spin(); // Needed for rotational flattening (torque and accelerations) and evolution
        }

        if (dspin_dt && (self.consider_tides || self.consider_rotational_flattening)) ||
            (accelerations && (self.consider_tides || self.consider_rotational_flattening || self.consider_general_relativy != ConsiderGeneralRelativity::None)) {
            self.calculate_distance_and_velocities(); // Needed for tides, rotational flattening, general relativity and evolution
        }

        if evolution && self.evolving_particles_exist {
            self.calculate_particles_evolving_quantities(current_time);
        }

        if dspin_dt && self.wind_effects_exist {
            self.calculate_wind_factor();
        }

        if (dspin_dt && (self.consider_tides || self.consider_rotational_flattening)) ||
            (accelerations && (self.consider_tides || self.consider_rotational_flattening || self.consider_general_relativy != ConsiderGeneralRelativity::None)) {

            self.calculate_orthogonal_components(); // Needed for torques and additional accelerations

            if dspin_dt && (self.consider_tides || self.consider_rotational_flattening) {
                // Not needed for additional accelerations
                self.calculate_torques(); // Needed for dspin_dt
                self.calculate_dspin_dt();
            }

            if accelerations && (self.consider_tides || self.consider_rotational_flattening || self.consider_general_relativy != ConsiderGeneralRelativity::None) {
                if self.consider_tides {
                    self.calculate_radial_component_of_the_tidal_force();  // Needed for calculate_tidal_acceleration
                    self.calculate_tidal_acceleration();
                }

                if self.consider_rotational_flattening {
                    self.calculate_radial_component_of_the_force_induced_by_rotational_flattening();
                    self.calculate_acceleration_induced_by_rotational_flattering();
                }

                match self.consider_general_relativy {
                    ConsiderGeneralRelativity::Kidder1995 => self.calculate_kidder1995_general_relativity_acceleration(),
                    ConsiderGeneralRelativity::Anderson1975 => self.calculate_anderson1975_general_relativity_acceleration(ignored_gravity_terms),
                    ConsiderGeneralRelativity::Newhall1983 => self.calculate_newhall1983_general_relativity_acceleration(ignored_gravity_terms),
                    ConsiderGeneralRelativity::None => {}
                }

                self.apply_acceleration_corrections();

                //for (i, particle) in self.particles[1..self.n_particles].iter().enumerate() {
                //for (i, particle) in self.particles[0..self.n_particles].iter().enumerate() {
                    //println!("{} - Acceleration {:e} {:e} {:e}", i, particle.acceleration.x, particle.acceleration.y, particle.acceleration.z);
                    //println!("{} - Acceleration rot {:e} {:e} {:e}", i, particle.acceleration_induced_by_rotational_flattering.x, particle.acceleration_induced_by_rotational_flattering.y, particle.acceleration_induced_by_rotational_flattering.z);
                //}
            }
        }
    }


    fn calculate_torques(&mut self) {
        let central_body = true;

        if self.consider_tides {
            self.calculate_torque_due_to_tides(central_body);
            self.calculate_torque_due_to_tides(!central_body);
        }

        if self.consider_rotational_flattening {
            self.calculate_torque_induced_by_rotational_flattening(central_body);
            self.calculate_torque_induced_by_rotational_flattening(!central_body);
        }
    }

    fn calculate_orthogonal_components(&mut self) {
        let central_body = true;

        self.calculate_planet_dependent_dissipation_factors(); // Needed by calculate_orthogonal_component_of_the_tidal_force and calculate_orthogonal_component_of_the_tidal_force if BolmontMathis2016/GalletBolmont2017
        self.calculate_scalar_product_of_vector_position_with_spin(); // Needed by tides and rotational flattening
        
        if self.consider_tides {
            self.calculate_orthogonal_component_of_the_tidal_force(central_body);
            self.calculate_orthogonal_component_of_the_tidal_force(!central_body);
        }

        if self.consider_rotational_flattening {
            self.calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(central_body);
            self.calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(!central_body);
        }
    }

    fn apply_acceleration_corrections(&mut self) {
        // Add the tidal+flattening+general relativity accelerations to the gravitational one (already computed)
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                if self.consider_tides {
                    particle.acceleration.x += particle.tidal_acceleration.x - star.tidal_acceleration.x;
                    particle.acceleration.y += particle.tidal_acceleration.y - star.tidal_acceleration.y;
                    particle.acceleration.z += particle.tidal_acceleration.z - star.tidal_acceleration.z;
                }

                if self.consider_rotational_flattening {
                    particle.acceleration.x += particle.acceleration_induced_by_rotational_flattering.x - star.acceleration_induced_by_rotational_flattering.x;
                    particle.acceleration.y += particle.acceleration_induced_by_rotational_flattering.y - star.acceleration_induced_by_rotational_flattering.y;
                    particle.acceleration.z += particle.acceleration_induced_by_rotational_flattering.z - star.acceleration_induced_by_rotational_flattering.z;
                }

                if self.consider_general_relativy != ConsiderGeneralRelativity::None {
                    particle.acceleration.x += particle.general_relativity_acceleration.x - star.general_relativity_acceleration.x;
                    particle.acceleration.y += particle.general_relativity_acceleration.y - star.general_relativity_acceleration.y;
                    particle.acceleration.z += particle.general_relativity_acceleration.z - star.general_relativity_acceleration.z;
                } 
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // TIDES
    fn calculate_torque_due_to_tides(&mut self, central_body:bool) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let mut torque = Axes{x: 0., y: 0., z:0.};
            let mut reference_spin = star.spin.clone();
            let mut orthogonal_component_of_the_tidal_force: f64;
            let mut reference_rscalspin: f64;

            for particle in particles.iter_mut() {
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
                let n_tid_x: f64 = orthogonal_component_of_the_tidal_force * (distance * reference_spin.x - reference_rscalspin*particle.position.x/distance - 1.0/distance
                                * (particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y) );

                let n_tid_y :f64 = orthogonal_component_of_the_tidal_force * (distance * reference_spin.y - reference_rscalspin*particle.position.y/distance - 1.0/distance
                                * (particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z) );

                let n_tid_z: f64 = orthogonal_component_of_the_tidal_force * (distance * reference_spin.z - reference_rscalspin*particle.position.z/distance - 1.0/distance
                                * (particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x) );

                if central_body {
                    // Integration of the spin (total torque tides):
                    //let factor = star.mass_g / (star.mass_g + particle.mass_g);
                    ////let factor = star.mass / (star.mass + particle.mass);
                    let factor = K2 / (star.mass_g + particle.mass_g); // Simplification because later on we will devide by star.mass_g [search SIMPLIFICATION]
                    //let factor = 1. / (star.mass + particle.mass);
                    torque.x += factor * n_tid_x;
                    torque.y += factor * n_tid_y;
                    torque.z += factor * n_tid_z;
                } else {
                    particle.torque_due_to_tides.x = n_tid_x;
                    particle.torque_due_to_tides.y = n_tid_y;
                    particle.torque_due_to_tides.z = n_tid_z;
                }

            }

            if central_body {
                // - Equation 25 from Bolmont et al. 2015
                star.torque_due_to_tides.x = torque.x;
                star.torque_due_to_tides.y = torque.y;
                star.torque_due_to_tides.z = torque.z;
            }
        }

    }


    fn calculate_dspin_dt(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            // - Equation 25 from Bolmont et al. 2015
            let factor = - 1. / (star.radius_of_gyration_2 * star.radius.powi(2)); // Simplification, we dont divide by star.mass_g because we did not mutliplied before (terms cancel out) [search SIMPLIFICATION]
            star.dspin_dt.x = factor * (star.torque_due_to_tides.x + star.torque_induced_by_rotational_flattening.x);
            star.dspin_dt.y = factor * (star.torque_due_to_tides.y + star.torque_induced_by_rotational_flattening.y);
            star.dspin_dt.z = factor * (star.torque_due_to_tides.z + star.torque_induced_by_rotational_flattening.z);

            for particle in particles.iter_mut() {
                //let factor1 = star.mass_g / (star.mass_g + particle.mass_g);
                //let factor2 = - K2 / (particle.mass_g * particle.radius_of_gyration_2 * particle.radius.powi(2));
                //let factor = factor1 * factor2;
                let factor = - K2 * star.mass_g / (particle.mass_g * (particle.mass_g + star.mass_g) 
                                * particle.radius_of_gyration_2 * particle.radius.powi(2));
                // - Equation 25 from Bolmont et al. 2015
                particle.dspin_dt.x = factor * (particle.torque_due_to_tides.x + particle.torque_induced_by_rotational_flattening.x);
                particle.dspin_dt.y = factor * (particle.torque_due_to_tides.y + particle.torque_induced_by_rotational_flattening.y);
                particle.dspin_dt.z = factor * (particle.torque_due_to_tides.z + particle.torque_induced_by_rotational_flattening.z);
            }
        }
    }

    fn calculate_orthogonal_component_of_the_tidal_force(&mut self, central_body:bool) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                // (distance to star)^7
                let distance_7 = particle.distance.powi(7);

                //// Tidal force calculation (star) :: Only orthogonal component is needed
                if central_body {
                    // - Third line of Equation 5 from Bolmont et al. 2015
                    //   This expression has R**10 (instead of R**5 in Eq. 5) 
                    //   because it uses sigma (i.e., scaled_dissipation_factor) 
                    //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
                    //   there is a R**5 factor as shown in Equation 28)
                    let star_scaled_dissipation_factor = Universe::planet_dependent_dissipation_factor(&self.star_planet_dependent_dissipation_factors, &star.id, star.evolution_type, star.scaled_dissipation_factor);
                    particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass_g.powi(2))
                                                    * (star.radius.powi(10)) 
                                                    * star_scaled_dissipation_factor / ( (K2.powi(2)) * distance_7);
                    //particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass.powi(2))
                                                        //* (star.radius.powi(10)) 
                                                        //* star.scaled_dissipation_factor / (distance_7);
                } else {
                    // - Second line of Equation 5 from Bolmont et al. 2015
                    //   This expression has R**10 (instead of R**5 in Eq. 5) 
                    //   because it uses sigma (i.e., scaled_dissipation_factor) 
                    //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
                    //   there is a R**5 factor as shown in Equation 28)
                    particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (star.mass_g.powi(2))
                                                    * (particle.radius.powi(10))
                                                    * particle.scaled_dissipation_factor    / ( (K2.powi(2)) * distance_7);
                    //particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (star.mass.powi(2))
                                                    //* (particle.radius.powi(10))
                                                    //* particle.scaled_dissipation_factor    / (distance_7);

                    // SBC
                    //println!("> {:e} {:e} {:e} {:e}", star.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
                    //println!("> {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z);
                }
            }
        }
    }

    fn calculate_radial_component_of_the_tidal_force(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let star_mass_2 = star.mass_g * star.mass_g;

            for particle in particles.iter_mut() {
                let planet_mass_2 = particle.mass_g * particle.mass_g;
                // Conservative part of the radial tidal force

                let radial_component_of_the_tidal_force_conservative_part = -3.0 / (particle.distance.powi(7) * K2)
                            * (planet_mass_2 * star.radius.powi(5) * star.love_number 
                            + star_mass_2 * particle.radius.powi(5) * particle.love_number);
                //let radial_component_of_the_tidal_force_conservative_part = -3.0 / particle.distance.powi(7)
                            //* (planet_mass_2 * star.radius.powi(5) * star.love_number 
                            //+ star_mass_2 * particle.radius.powi(5) * particle.love_number);

                // Dissipative part of the radial tidal force
                let factor1 = -13.5 * particle.radial_velocity / (particle.distance.powi(8) * K2*K2);
                //let factor1 = -13.5 * particle.radial_velocity / particle.distance.powi(8);
                let star_scaled_dissipation_factor = Universe::planet_dependent_dissipation_factor(&self.star_planet_dependent_dissipation_factors, &particle.id, star.evolution_type, star.scaled_dissipation_factor);
                let term1 = planet_mass_2
                            * star.radius.powi(10)
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
    }

    pub fn calculate_denergy_dt(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                // - Equation 32 from Bolmont et al. 2015
                //// Instantaneous energy loss dE/dt due to tides
                //// in Msun.AU^2.day^(-3)
                //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
                let factor2 = particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance;
                particle.denergy_dt = -((1.0 / particle.distance * (particle.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.radial_velocity))
                            * (particle.position.x*particle.velocity.x + particle.position.y*particle.velocity.y + particle.position.z*particle.velocity.z)
                            + factor2 
                            * ((particle.spin.y*particle.position.z - particle.spin.z*particle.position.y - particle.velocity.x) * particle.velocity.x
                            + (particle.spin.z*particle.position.x - particle.spin.x*particle.position.z - particle.velocity.y) * particle.velocity.y
                            + (particle.spin.x*particle.position.y - particle.spin.y*particle.position.x - particle.velocity.z) * particle.velocity.z))
                            + star.mass_g/(star.mass_g + particle.mass_g) 
                            * (particle.torque_due_to_tides.x*particle.spin.x + particle.torque_due_to_tides.y*particle.spin.y + particle.torque_due_to_tides.z*particle.spin.z);
            }
        }
    }


    fn calculate_tidal_acceleration(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let factor2 = K2 / star.mass_g;
            //let factor2 = 1 / star.mass;
            let mut sum_total_tidal_force = Axes{x:0., y:0., z:0.};

            for particle in particles.iter_mut() {
                let factor1 = K2 / particle.mass_g;
                //let factor1 = 1. / particle.mass;

                // - Equation 6 from Bolmont et al. 2015
                let factor3 = particle.radial_component_of_the_tidal_force
                                + (particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide) * particle.radial_velocity / particle.distance;
                let total_tidal_force_x = factor3 * particle.position.x / particle.distance
                                        + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                            * (star.spin.y * particle.position.z  - star.spin.z * particle.position.y - particle.velocity.x)
                                        + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                            * (particle.spin.y * particle.position.z  - particle.spin.z * particle.position.y - particle.velocity.x);
                let total_tidal_force_y = factor3 * particle.position.y / particle.distance
                                        + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                            * (star.spin.z * particle.position.x  - star.spin.x * particle.position.z - particle.velocity.y)
                                        + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                            * (particle.spin.z * particle.position.x  - particle.spin.x * particle.position.z - particle.velocity.y);
                let total_tidal_force_z = factor3 * particle.position.z / particle.distance
                                        + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                            * (star.spin.x * particle.position.y  - star.spin.y * particle.position.x - particle.velocity.z)
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
            star.tidal_acceleration.x = -1.0 * factor2 * sum_total_tidal_force.x;
            star.tidal_acceleration.y = -1.0 * factor2 * sum_total_tidal_force.y;
            star.tidal_acceleration.z = -1.0 * factor2 * sum_total_tidal_force.z;
        }
    }
    







    fn calculate_kidder1995_general_relativity_acceleration(&mut self) {

        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let factor2 = K2 / star.mass_g;
            //let factor2 = 1 / star.mass;
            let mut sum_total_general_relativity_force = Axes{x:0., y:0., z:0.};

            for particle in particles.iter_mut() {
                let factor1 = K2 / particle.mass_g;
                //let factor1 = 1. / particle.mass;

                // Radial part of the GR force (Kidder 1995, Mardling & Lin 2002)
                // - Equation 11 from Bolmont et al. 2015
                let star_planet_mass_g = star.mass_g + particle.mass_g;
                let distance_2 = particle.distance.powi(2);
                let radial_velocity_2 = particle.radial_velocity.powi(2);
                let radial_component_of_the_general_relativity_force = -star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                    * ( (1.0 + 3.0 * particle.general_relativity_factor) * particle.norm_velocity_vector_2
                    -2.0 * (2.0 + particle.general_relativity_factor) * star_planet_mass_g/particle.distance
                    -1.5 * particle.general_relativity_factor * radial_velocity_2);
                //println!("Radial component GR force {:e}", radial_component_of_the_general_relativity_force);
                // Orthoradial part of the GR force
                // - Equation 11 from Bolmont et al. 2015
                let orthogonal_component_of_the_general_relativity_force = star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                    * 2.0 * (2.0 - particle.general_relativity_factor) * particle.radial_velocity * particle.norm_velocity_vector;
                //println!("Ortho component GR force {:e}", orthogonal_component_of_the_general_relativity_force);
                // Total General Relativity force
                // - Equation 10 from Bolmont et al. 2015
                let total_general_relativity_force_x = particle.mass 
                        * (radial_component_of_the_general_relativity_force * particle.position.x / particle.distance
                        + orthogonal_component_of_the_general_relativity_force * particle.velocity.x / particle.norm_velocity_vector);
                let total_general_relativity_force_y = particle.mass 
                        * (radial_component_of_the_general_relativity_force * particle.position.y / particle.distance
                        + orthogonal_component_of_the_general_relativity_force * particle.velocity.y / particle.norm_velocity_vector);
                let total_general_relativity_force_z = particle.mass 
                        * (radial_component_of_the_general_relativity_force * particle.position.z / particle.distance
                        + orthogonal_component_of_the_general_relativity_force * particle.velocity.z / particle.norm_velocity_vector);
                
                sum_total_general_relativity_force.x += total_general_relativity_force_x;
                sum_total_general_relativity_force.y += total_general_relativity_force_y;
                sum_total_general_relativity_force.z += total_general_relativity_force_z;

                // - Equation 19 from Bolmont et al. 2015 (first term)
                particle.general_relativity_acceleration.x = factor1 * total_general_relativity_force_x;
                particle.general_relativity_acceleration.y = factor1 * total_general_relativity_force_y;
                particle.general_relativity_acceleration.z = factor1 * total_general_relativity_force_z;
                //println!("GR force {:e} {:e} {:e}", total_general_relativity_force_x, total_general_relativity_force_y, total_general_relativity_force_z);
            }
            
            // - Equation 19 from Bolmont et al. 2015 (second term)
            //for particle in particles.iter_mut() {
                //particle.general_relativity_acceleration.x += factor2 * sum_total_general_relativity_force.x;
                //particle.general_relativity_acceleration.y += factor2 * sum_total_general_relativity_force.y;
                //particle.general_relativity_acceleration.z += factor2 * sum_total_general_relativity_force.z;
            //}
            // Instead of the previous code, keep star tidal acceleration separated:
            star.general_relativity_acceleration.x = -1.0 * factor2 * sum_total_general_relativity_force.x;
            star.general_relativity_acceleration.y = -1.0 * factor2 * sum_total_general_relativity_force.y;
            star.general_relativity_acceleration.z = -1.0 * factor2 * sum_total_general_relativity_force.z;
        }
    }

    fn calculate_scalar_product_of_vector_position_with_spin (&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.scalar_product_of_vector_position_with_stellar_spin = particle.position.x * star.spin.x 
                                + particle.position.y * star.spin.y
                                + particle.position.z * star.spin.z;
            
                particle.scalar_product_of_vector_position_with_planetary_spin = particle.position.x * particle.spin.x 
                                + particle.position.y * particle.spin.y
                                + particle.position.z * particle.spin.z;
            }
        }
    }

    /// 
    fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(&mut self, central_body:bool) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            // Calculation of the norm square of the spin for the planet
            for particle in particles.iter_mut() {
                if central_body {
                    // - Star Equation 16 from Bolmont et al. 2015
                    particle.factor_for_the_force_induced_by_star_rotation = particle.mass_g * star.fluid_love_number * star.norm_spin_vector_2 * star.radius.powi(5) / (6. * K2); // Msun.AU^5.day-2
                    // - Second part of Equation 15 from Bolmont et al. 2015
                    particle.orthogonal_component_of_the_force_induced_by_star_rotation = -6. * particle.factor_for_the_force_induced_by_star_rotation * particle.scalar_product_of_vector_position_with_stellar_spin / (star.norm_spin_vector_2 * particle.distance.powi(5));
                } else {
                    // - Planet Equation 16 from Bolmont et al. 2015
                    particle.factor_for_the_force_induced_by_planet_rotation = star.mass_g * particle.fluid_love_number * particle.norm_spin_vector_2 * particle.radius.powi(5) / (6. * K2); // Msun.AU^5.day-2
                    // - Second part of Equation 15 from Bolmont et al. 2015
                    particle.orthogonal_component_of_the_force_induced_by_planet_rotation = -6. * particle.factor_for_the_force_induced_by_planet_rotation * particle.scalar_product_of_vector_position_with_planetary_spin / (particle.norm_spin_vector_2 * particle.distance.powi(5));
                }
            }
        }
    }

    fn calculate_radial_component_of_the_force_induced_by_rotational_flattening(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                // - First part of Equation 15 from Bolmont et al. 2015
                particle.radial_component_of_the_force_induced_by_rotation = -3./particle.distance.powi(5) * (particle.factor_for_the_force_induced_by_planet_rotation + particle.factor_for_the_force_induced_by_star_rotation)
                    + 15./particle.distance.powi(7) * (particle.factor_for_the_force_induced_by_star_rotation * particle.scalar_product_of_vector_position_with_stellar_spin * particle.scalar_product_of_vector_position_with_stellar_spin/star.norm_spin_vector_2
                        + particle.factor_for_the_force_induced_by_planet_rotation * particle.scalar_product_of_vector_position_with_planetary_spin * particle.scalar_product_of_vector_position_with_planetary_spin/particle.norm_spin_vector_2); // Msun.AU.day-1
            }
        }
    }

    fn calculate_torque_induced_by_rotational_flattening(&mut self, central_body:bool) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let mut torque = Axes{x: 0., y: 0., z:0.};
            let mut reference_spin = star.spin.clone();
            let mut orthogonal_component_of_the_force_induced_by_rotation: f64;

            for particle in particles.iter_mut() {
                if !central_body {
                    reference_spin = particle.spin.clone();
                    orthogonal_component_of_the_force_induced_by_rotation = particle.orthogonal_component_of_the_force_induced_by_planet_rotation;
                } else {
                    orthogonal_component_of_the_force_induced_by_rotation = particle.orthogonal_component_of_the_force_induced_by_star_rotation;
                }
                
                //// Torque calculation due to rotational flattening
                // - Equation 17-18 from Bolmont et al. 2015
                let n_rot_x: f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.position.y * reference_spin.z - particle.position.z * reference_spin.y);
                let n_rot_y :f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.position.z * reference_spin.x - particle.position.x * reference_spin.z);
                let n_rot_z: f64 = orthogonal_component_of_the_force_induced_by_rotation * (particle.position.x * reference_spin.y - particle.position.y * reference_spin.x);

                // - Equation 25 from Bolmont et al. 2015
                if central_body {
                    // Integration of the spin (total torque rot):
                    //let factor = star.mass_g / (star.mass_g + particle.mass_g);
                    ////let factor = star.mass / (star.mass + particle.mass);
                    let factor = K2 / (star.mass_g + particle.mass_g); // Simplification because later on we will divide by star.mass_g [search SIMPLIFICATION]
                    //let factor = 1. / (star.mass + particle.mass);
                    torque.x += factor * n_rot_x;
                    torque.y += factor * n_rot_y;
                    torque.z += factor * n_rot_z;
                } else {
                    particle.torque_induced_by_rotational_flattening.x = n_rot_x;
                    particle.torque_induced_by_rotational_flattening.y = n_rot_y;
                    particle.torque_induced_by_rotational_flattening.z = n_rot_z;
                }
            }

            if central_body {
                // - Equation 25 from Bolmont et al. 2015
                star.torque_induced_by_rotational_flattening.x = torque.x;
                star.torque_induced_by_rotational_flattening.y = torque.y;
                star.torque_induced_by_rotational_flattening.z = torque.z;
            }
        }
    }

    fn calculate_acceleration_induced_by_rotational_flattering(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let factor2 = K2 / star.mass_g;
            //let factor2 = 1 / star.mass;

            // Calculation of the norm square of the spin for the star
            let mut sum_total_force_induced_by_rotation = Axes{x:0., y:0., z:0.};

            for particle in particles.iter_mut() {
                let factor1 = K2 / particle.mass_g;

                // - Equation 15 from Bolmont et al. 2015
                let total_force_induced_by_rotation_x = particle.radial_component_of_the_force_induced_by_rotation * particle.position.x
                    + particle.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.x
                    + particle.orthogonal_component_of_the_force_induced_by_star_rotation * star.spin.x;
                let total_force_induced_by_rotation_y = particle.radial_component_of_the_force_induced_by_rotation * particle.position.y
                    + particle.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.y
                    + particle.orthogonal_component_of_the_force_induced_by_star_rotation * star.spin.y;
                let total_force_induced_by_rotation_z = particle.radial_component_of_the_force_induced_by_rotation * particle.position.z
                    + particle.orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.z
                    + particle.orthogonal_component_of_the_force_induced_by_star_rotation * star.spin.z;
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
            star.acceleration_induced_by_rotational_flattering.x = -1.0 * factor2 * sum_total_force_induced_by_rotation.x;
            star.acceleration_induced_by_rotational_flattering.y = -1.0 * factor2 * sum_total_force_induced_by_rotation.y;
            star.acceleration_induced_by_rotational_flattering.z = -1.0 * factor2 * sum_total_force_induced_by_rotation.z;
        }
    }




    pub fn calculate_norm_spin(&mut self) {
        for particle in self.particles[..self.n_particles].iter_mut() {
            // Squared norm of the spin
            particle.norm_spin_vector_2 = (particle.spin.x.powi(2)) 
                                + (particle.spin.y.powi(2))
                                + (particle.spin.z.powi(2));
        }
    }

    fn calculate_distance_and_velocities(&mut self) {
        // Calculation of velocity vv(j), radial velocity vrad(j)
        // velocities in AU/day
        for particle in self.particles[1..self.n_particles].iter_mut() {
            // Norm of the velocity
            particle.norm_velocity_vector_2 = (particle.velocity.x.powi(2)) 
                                + (particle.velocity.y.powi(2))
                                + (particle.velocity.z.powi(2));
            particle.norm_velocity_vector = particle.norm_velocity_vector_2.sqrt();

            // (distance to star)^2
            let distance_2 = (particle.position.x.powi(2)) 
                                + (particle.position.y.powi(2))
                                + (particle.position.z.powi(2));
            let distance = distance_2.sqrt();
            particle.distance = distance;

            // Radial velocity
            let v_rad = (particle.position.x*particle.velocity.x +
                        particle.position.y*particle.velocity.y +
                        particle.position.z*particle.velocity.z)
                        / distance;
            particle.radial_velocity = v_rad;
            //let tmp = particle.position.x*particle.velocity.x +
                        //particle.position.y*particle.velocity.y +
                        //particle.position.z*particle.velocity.z;
            //println!("{:e} {:e} {:e}", particle.norm_velocity_vector, particle.distance, particle.radial_velocity);
            //println!("{:e} {:e} {:e}", particle.radial_velocity, tmp, tmp/distance); 
            //println!(">p {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z);
            //println!(">v {:e} {:e} {:e}", particle.velocity.x, particle.velocity.y, particle.velocity.z);
        }
    }

    fn calculate_planet_dependent_dissipation_factors(&mut self) {
        let star_index = 0; // index
        match self.particles[star_index].evolution_type {
            EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) => {
                if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
                    for particle in particles.iter() {
                        let frequency = (particle.velocity.x - star.spin.y*particle.position.z + star.spin.z*particle.position.y).powi(2)
                                    + (particle.velocity.y - star.spin.z*particle.position.x + star.spin.x*particle.position.z).powi(2)
                                    + (particle.velocity.z - star.spin.x*particle.position.y + star.spin.y*particle.position.x).powi(2);
                        // two_times_the_inverse_of_the_excitation_frequency: 2/w
                        // inverse_of_half_the_excitation_frequency : 1/(w/2)
                        let inverse_of_half_the_excitation_frequency = particle.distance / frequency;
                        let planet_dependent_dissipation_factor = star.dissipation_factor_scale * 2.0 * K2
                            * star.lag_angle * inverse_of_half_the_excitation_frequency / (3.0*star.radius.powi(5));

                        self.star_planet_dependent_dissipation_factors.insert(particle.id.clone(), planet_dependent_dissipation_factor);
                        //println!("Insert {} in {}", planet_dependent_dissipation_factor, particle.id);
                    }
                }
            },
            _ => {},
        }

    }

    pub fn planet_dependent_dissipation_factor(star_planet_dependent_dissipation_factors: &HashMap<usize, f64>,  id: &usize, evolution_type: EvolutionType, scaled_dissipation_factor: f64) -> f64 {
        match evolution_type {
            EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) => {
                match star_planet_dependent_dissipation_factors.get(id) {
                    Some(&value) => value,
                    _ => scaled_dissipation_factor // This should not happen
                }
            },
            _ => scaled_dissipation_factor,
        }
    }

    pub fn compute_total_energy(&self) -> f64 {
        let mut e_kin = 0.;
        let mut e_pot = 0.;
        let e_offset = 0.; // Energy offset due to collisions and ejections

        // Kinectic energy
        for particle in self.particles[..self.n_particles].iter() {
            e_kin += 0.5 * particle.mass * (particle.velocity.x.powi(2) + particle.velocity.y.powi(2) + particle.velocity.z.powi(2));
        }
        // Gravitationl potential energy
        for (i, particle_a) in self.particles[..self.n_particles].iter().enumerate() {
            for particle_b in self.particles[i+1..self.n_particles].iter() {
                let dx = particle_a.position.x - particle_b.position.x;
                let dy = particle_a.position.y - particle_b.position.y;
                let dz = particle_a.position.z - particle_b.position.z;
                e_pot -= particle_b.mass_g*particle_a.mass/(dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
            }
        }
        
        e_kin + e_pot + e_offset
    }

    pub fn compute_total_angular_momentum(&self) -> f64 {
        let mut total_angular_momentum = Axes{x:0., y:0., z:0.}; // L
        for particle in self.particles[..self.n_particles].iter() {
            total_angular_momentum.x += particle.mass*(particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y);
            total_angular_momentum.y += particle.mass*(particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z);
            total_angular_momentum.z += particle.mass*(particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x);
        }
        let total_angular_momentum = (total_angular_momentum.x.powf(2.) + total_angular_momentum.y.powf(2.) + total_angular_momentum.z.powf(2.)).sqrt();
        total_angular_momentum
    }

    pub fn calculate_wind_factor(&mut self) {
        for particle in self.particles[..self.n_particles].iter_mut() {
            if particle.wind_k_factor != 0. {
                let threshold = (particle.norm_spin_vector_2).sqrt();
                let factor = - K2 / (particle.mass_g * particle.radius_of_gyration_2 * particle.radius.powi(2));
                if threshold >= particle.wind_rotation_saturation {
                    // Friendly reminder that m(1) is in solar mass * K2
                    //tmp2 = hdt * K_wind * wsat_wind*wsat_wind * sqrt(Rsth*K2/(Rsun*m(1)))
                    particle.wind_factor = factor * particle.wind_k_factor * particle.wind_rotation_saturation_2 * (particle.radius/R_SUN * 1./particle.mass).sqrt()
                } else {
                    //tmp2 = hdt * K_wind * sqrt(Rsth*K2/(Rsun*m(1)))
                    particle.wind_factor = factor * particle.wind_k_factor * (particle.radius/R_SUN * 1./particle.mass).sqrt()
                }
            };
        }
    }
    
    pub fn calculate_particles_evolving_quantities(&mut self, current_time: f64) {
        for (particle, evolver) in self.particles[..self.n_particles].iter_mut().zip(self.particles_evolvers.iter_mut()) {
            ////////////////////////////////////////////////////////////////////
            // Radius and radius of gyration 2
            ////////////////////////////////////////////////////////////////////
            // If particle radius/radius_of_gyration_2 evolves
            // - Compute ratio between previous and new moment of inertia 
            // - The ratio will be used in the spin integration
            let new_radius = evolver.radius(current_time, particle.radius);
            let new_radius_of_gyration_2 = evolver.radius_of_gyration_2(current_time, particle.radius_of_gyration_2);
            if new_radius != particle.radius || new_radius_of_gyration_2 != particle.radius_of_gyration_2 {
                // Update moment of inertia ratio only if it is not during first initialization
                if current_time > 0. {
                    particle.moment_of_inertia_ratio = (particle.radius_of_gyration_2 * particle.radius.powi(2)) 
                                                        / (new_radius_of_gyration_2 * new_radius.powi(2));
                }
                particle.radius = new_radius;
                particle.radius_of_gyration_2 = new_radius_of_gyration_2;
            } else {
                particle.moment_of_inertia_ratio = 1.;
            }

            ////////////////////////////////////////////////////////////////////
            // Love number
            ////////////////////////////////////////////////////////////////////
            particle.love_number = evolver.love_number(current_time, particle.love_number);

            ////////////////////////////////////////////////////////////////////
            // Lag angle
            ////////////////////////////////////////////////////////////////////
            particle.lag_angle = match evolver.evolution_type {
                EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) => {
                        let inverse_tidal_q_factor = evolver.inverse_tidal_q_factor(current_time, 0.);
                        let epsilon_squared = particle.norm_spin_vector_2/SUN_DYN_FREQ;
                        // Normal formula = 3.d0*epsilon_squared*Q_str_inv/(4.d0*k2s)
                        // but as for sigma it is necessary to divide by k2s, we do not divide here
                        let lag_angle = 3.0*epsilon_squared*inverse_tidal_q_factor/4.0;
                        lag_angle
                    },
                _ => 0.,
            };
            //println!("[{}] Evolve Radius {:e} Gyration {:e} Love {:e} Lag {:e}", current_time, particle.radius, particle.radius_of_gyration_2, particle.love_number, particle.lag_angle);
        }

    }

    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
    // [start] General Relativity based on REBOUNDx gr.c
    fn calculate_anderson1975_general_relativity_acceleration(&mut self, ignored_gravity_terms: IgnoreGravityTerms) {
        // Calculate Newtonian accelerations in the current setup and considering all particles
        let newtonian_inertial_accelerations = self.get_anderson1975_newhall1983_newtonian_inertial_accelerations(ignored_gravity_terms);

        // Transform to Jacobi coordinates
        let (jacobi_star_mass, _jacobi_star_position, _jacobi_star_velocity, _jacobi_star_acceleration, jacobi_particles_positions, jacobi_particles_velocities, mut jacobi_particles_accelerations) = self.anderson1975_general_relativity_inertial_to_jacobi_posvelacc(newtonian_inertial_accelerations);

        //println!("Star {:?}", jacobi_star_mass);
        //println!("---");

        let mu = self.particles[0].mass_g;
        //println!("$$$$$$$$$$ {:e}", mu);
        for ((jacobi_particle_acceleration, jacobi_particle_velocity), jacobi_particle_position) in jacobi_particles_accelerations[..self.n_particles-1].iter_mut()
                                                                                                        .zip(jacobi_particles_velocities[..self.n_particles-1].iter())
                                                                                                        .zip(jacobi_particles_positions[..self.n_particles-1].iter()) {
            //println!("Transformed Planets a {:?}", jacobi_particle_acceleration);
            //println!("Transformed Planets x {:?}", jacobi_particle_position);
            //println!("Transformed Planets v {:?}", jacobi_particle_velocity);
            //println!("...");
            let mut vi = Axes{x: jacobi_particle_velocity.x, y: jacobi_particle_velocity.y, z: jacobi_particle_velocity.z};
            let mut vi2 = jacobi_particle_velocity.x.powi(2) + jacobi_particle_velocity.y.powi(2) + jacobi_particle_velocity.z.powi(2);
            let ri = (jacobi_particle_position.x.powi(2) + jacobi_particle_position.y.powi(2) + jacobi_particle_position.z.powi(2)).sqrt();
            let mut factor_a = (0.5*vi2 + 3.*mu/ri)/SPEED_OF_LIGHT_2;
            let mut old_v = Axes{x:0., y:0., z:0.};

            let max_iterations = 10;
            //println!("$$$$$ a {:e}", factor_a);
            for q in 0..max_iterations {
                old_v.x = vi.x;
                old_v.y = vi.y;
                old_v.z = vi.z;
                vi.x = jacobi_particle_velocity.x/(1.-factor_a);
                vi.y = jacobi_particle_velocity.y/(1.-factor_a);
                vi.z = jacobi_particle_velocity.z/(1.-factor_a);
                vi2 =vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
                factor_a = (0.5*vi2 + 3.*mu/ri)/SPEED_OF_LIGHT_2;
                //println!("$$$$$ a {:e}", factor_a);
                let dvx = vi.x - old_v.x;
                let dvy = vi.y - old_v.y;
                let dvz = vi.z - old_v.z;
                if (dvx*dvx + dvy*dvy + dvz*dvz)/vi2 < DBL_EPSILON_2 {
                    break;
                } else if q == max_iterations {
                    println!("[WARNING {} UTC] {} iterations in general relativity failed to converge. This is typically because the perturbation is too strong for the current implementation.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), max_iterations);
                }
            }

            let factor_b = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/SPEED_OF_LIGHT_2;
            let rdotrdot = jacobi_particle_position.x*jacobi_particle_velocity.x
                            + jacobi_particle_position.y*jacobi_particle_velocity.y
                            + jacobi_particle_position.z*jacobi_particle_velocity.z;
            let vidot = Axes{x: jacobi_particle_acceleration.x + factor_b*jacobi_particle_position.x,
                                y: jacobi_particle_acceleration.y + factor_b*jacobi_particle_position.y,
                                z: jacobi_particle_acceleration.z + factor_b*jacobi_particle_position.z};
            let vdotvdot = vi.x*vidot.x + vi.y*vidot.y + vi.z*vidot.z;
            let factor_d = (vdotvdot - 3.*mu/(ri*ri*ri)*rdotrdot)/SPEED_OF_LIGHT_2;
            //println!("$$$$$ factors: {:e} {:e} {:e}", factor_b, factor_a, factor_d);
            jacobi_particle_acceleration.x = factor_b*(1.-factor_a)*jacobi_particle_position.x - factor_a*jacobi_particle_acceleration.x - factor_d*vi.x;
            jacobi_particle_acceleration.y = factor_b*(1.-factor_a)*jacobi_particle_position.y - factor_a*jacobi_particle_acceleration.y - factor_d*vi.y;
            jacobi_particle_acceleration.z = factor_b*(1.-factor_a)*jacobi_particle_position.z - factor_a*jacobi_particle_acceleration.z - factor_d*vi.z;
        }


        let jacobi_star_acceleration = Axes{x:0., y:0., z:0.};
        let (star_acceleration, particles_accelerations) = self.anderson1975_general_relativity_jacobi_to_inertial_acc(jacobi_star_mass, jacobi_star_acceleration, jacobi_particles_accelerations);


        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            // This algorithm computes general_relativity_acceleration in the inertial frame
            // we need to convert it to the heliocentric frame to be consistent with the rest of
            // additional forces by:
            //
            // inertial_position = heliocentric_position - (1/total_mass) * sum(particle.mass * particle.position)
            // inertial_velocity = heliocentric_velocity - (1/total_mass) * sum(particle.mass * particle.velocity)
            // inertial_acceleration = heliocentric_acceleration - (1/total_mass) * sum(particle.mass * particle.acceleration)
            //
            // heliocentric_position = inertial_position + (1/star.mass) * sum(particle.mass * particle.position)
            // heliocentric_velocity = inertial_velocity + (1/star.mass) * sum(particle.mass * particle.velocity)
            // heliocentric_acceleration = inertial_acceleration + (1/star.mass) * sum(particle.mass * particle.acceleration)
            //
            // 1.- Compute second term using only general relativity accelerations
            let mut reference_gr_acceleration = Axes{x:0., y:0., z:0.};
            reference_gr_acceleration.x  = star_acceleration.x*star.mass;
            reference_gr_acceleration.y  = star_acceleration.y*star.mass;
            reference_gr_acceleration.z  = star_acceleration.z*star.mass;
            for (particle, particle_acceleration) in particles.iter().zip(particles_accelerations.iter()) {
                reference_gr_acceleration.x += particle_acceleration.x*particle.mass;
                reference_gr_acceleration.y += particle_acceleration.y*particle.mass;
                reference_gr_acceleration.z += particle_acceleration.z*particle.mass;
            }
            reference_gr_acceleration.x /= star.mass;
            reference_gr_acceleration.y /= star.mass;
            reference_gr_acceleration.z /= star.mass;
            // 2.- Convert from inertial to heliocentric
            for (particle, particle_acceleration) in particles.iter_mut().zip(particles_accelerations[..self.n_particles-1].iter()){
                particle.general_relativity_acceleration.x = particle_acceleration.x + reference_gr_acceleration.x;
                particle.general_relativity_acceleration.y = particle_acceleration.y + reference_gr_acceleration.y;
                particle.general_relativity_acceleration.z = particle_acceleration.z + reference_gr_acceleration.z;
            }
            star.general_relativity_acceleration.x = star_acceleration.x + reference_gr_acceleration.x;
            star.general_relativity_acceleration.y = star_acceleration.y + reference_gr_acceleration.y;
            star.general_relativity_acceleration.z = star_acceleration.z + reference_gr_acceleration.z;
        }
    }

    fn anderson1975_general_relativity_inertial_to_jacobi_posvelacc(&mut self, newtonian_inertial_accelerations: [Axes; MAX_PARTICLES]) -> (f64, Axes, Axes, Axes, [Axes; MAX_PARTICLES-1], [Axes; MAX_PARTICLES-1], [Axes; MAX_PARTICLES-1]) {
        let mut jacobi_star_mass = 0.;
        let mut jacobi_star_position = Axes{x:0., y:0., z:0. };
        let mut jacobi_star_velocity = Axes{x:0., y:0., z:0. };
        let mut jacobi_star_acceleration = Axes{x:0., y:0., z:0. };
        let mut jacobi_particles_positions = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];
        let mut jacobi_particles_velocities = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];
        let mut jacobi_particles_accelerations = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];

        if let Some((star_newtonian_acceleration, particles_newtonian_accelerations)) = newtonian_inertial_accelerations[..self.n_particles].split_first() {
            if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
                let m0 = star.mass;
                let mut eta = m0;
                let mut s_x = eta * star.inertial_position.x;
                let mut s_y = eta * star.inertial_position.y;
                let mut s_z = eta * star.inertial_position.z;
                let mut s_vx = eta * star.inertial_velocity.x;
                let mut s_vy = eta * star.inertial_velocity.y;
                let mut s_vz = eta * star.inertial_velocity.z;
                let mut s_ax = eta * star_newtonian_acceleration.x;
                let mut s_ay = eta * star_newtonian_acceleration.y;
                let mut s_az = eta * star_newtonian_acceleration.z;
                for ((((jacobi_particle_position, jacobi_particle_velocity), jacobi_particle_acceleration), particle), particle_newtonian_acceleration) in jacobi_particles_positions[..self.n_particles-1].iter_mut()
                                                                        .zip(jacobi_particles_velocities[..self.n_particles-1].iter_mut())
                                                                            .zip(jacobi_particles_accelerations[..self.n_particles-1].iter_mut())
                                                                            .zip(particles.iter())
                                                                            .zip(particles_newtonian_accelerations.iter()) {
                    let ei = 1./eta;
                    eta += particle.mass;
                    let pme = eta*ei;
                    jacobi_particle_position.x = particle.inertial_position.x - s_x*ei;
                    jacobi_particle_position.y = particle.inertial_position.y - s_y*ei;
                    jacobi_particle_position.z = particle.inertial_position.z - s_z*ei;
                    jacobi_particle_velocity.x = particle.inertial_velocity.x - s_vx*ei;
                    jacobi_particle_velocity.y = particle.inertial_velocity.y - s_vy*ei;
                    jacobi_particle_velocity.z = particle.inertial_velocity.z - s_vz*ei;
                    jacobi_particle_acceleration.x = particle_newtonian_acceleration.x - s_ax*ei;
                    jacobi_particle_acceleration.y = particle_newtonian_acceleration.y - s_ay*ei;
                    jacobi_particle_acceleration.z = particle_newtonian_acceleration.z - s_az*ei;
                    s_x  = s_x  * pme + particle.mass*jacobi_particle_position.x ;
                    s_y  = s_y  * pme + particle.mass*jacobi_particle_position.y ;
                    s_z  = s_z  * pme + particle.mass*jacobi_particle_position.z ;
                    s_vx = s_vx * pme + particle.mass*jacobi_particle_velocity.x;
                    s_vy = s_vy * pme + particle.mass*jacobi_particle_velocity.y;
                    s_vz = s_vz * pme + particle.mass*jacobi_particle_velocity.z;
                    s_ax = s_ax * pme + particle.mass*jacobi_particle_acceleration.x;
                    s_ay = s_ay * pme + particle.mass*jacobi_particle_acceleration.y;
                    s_az = s_az * pme + particle.mass*jacobi_particle_acceleration.z;
                }
                let mtot = eta;
                let mtot_i = 1./mtot;
                jacobi_star_mass = mtot;
                jacobi_star_position.x = s_x * mtot_i;
                jacobi_star_position.y = s_y * mtot_i;
                jacobi_star_position.z = s_z * mtot_i;
                jacobi_star_velocity.x = s_vx * mtot_i;
                jacobi_star_velocity.y = s_vy * mtot_i;
                jacobi_star_velocity.z = s_vz * mtot_i;
                jacobi_star_acceleration.x = s_ax * mtot_i;
                jacobi_star_acceleration.y = s_ay * mtot_i;
                jacobi_star_acceleration.z = s_az * mtot_i;
            }
        }
        return(jacobi_star_mass, jacobi_star_position, jacobi_star_velocity, jacobi_star_acceleration, jacobi_particles_positions, jacobi_particles_velocities, jacobi_particles_accelerations)
    }

    fn anderson1975_general_relativity_jacobi_to_inertial_acc(&mut self, jacobi_star_mass: f64, jacobi_star_acceleration: Axes, jacobi_particles_accelerations: [Axes; MAX_PARTICLES-1]) -> (Axes, [Axes; MAX_PARTICLES-1]) {
        let mut star_acceleration = Axes{x:0., y:0., z:0. };
        let mut particles_accelerations = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];

        let mut eta = jacobi_star_mass;
        let mut s_ax = eta * jacobi_star_acceleration.x;
        let mut s_ay = eta * jacobi_star_acceleration.y;
        let mut s_az = eta * jacobi_star_acceleration.z;
        for ((particle_acceleration, particle), jacobi_particle_acceleration) in particles_accelerations[..self.n_particles-1].iter_mut().rev()
                                                                                    .zip(self.particles[1..self.n_particles].iter().rev())
                                                                                    .zip((jacobi_particles_accelerations[..self.n_particles-1].iter().rev())) {
            let ei = 1./eta;
            s_ax = (s_ax - particle.mass*jacobi_particle_acceleration.x) * ei;
            s_ay = (s_ay - particle.mass*jacobi_particle_acceleration.y) * ei;
            s_az = (s_az - particle.mass*jacobi_particle_acceleration.z) * ei;
            particle_acceleration.x = jacobi_particle_acceleration.x + s_ax;
            particle_acceleration.y = jacobi_particle_acceleration.y + s_ay;
            particle_acceleration.z = jacobi_particle_acceleration.z + s_az;
            eta -= particle.mass;
            s_ax *= eta;
            s_ay *= eta;
            s_az *= eta;
        }
        let mtot = eta;
        let mtot_i = 1./mtot;
        star_acceleration.x = s_ax * mtot_i;
        star_acceleration.y = s_ay * mtot_i;
        star_acceleration.z = s_az * mtot_i;
        return(star_acceleration, particles_accelerations)
    }
    // [end] General Relativity based on REBOUNDx gr.c
    //--------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////

    fn get_anderson1975_newhall1983_newtonian_inertial_accelerations(&self, ignored_gravity_terms: IgnoreGravityTerms) -> [Axes; MAX_PARTICLES] {
        let mut newtonian_inertial_accelerations = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; 
        for (newtonian_acceleration, particle) in newtonian_inertial_accelerations.iter_mut().zip(self.particles[..self.n_particles].iter()) {
            newtonian_acceleration.x = particle.inertial_acceleration.x;
            newtonian_acceleration.y = particle.inertial_acceleration.y;
            newtonian_acceleration.z = particle.inertial_acceleration.z;
        }
        // If some terms where ignored by the integrator, they should be added
        if ignored_gravity_terms == IgnoreGravityTerms::WHFastOne || ignored_gravity_terms == IgnoreGravityTerms::WHFastTwo {
            let n_particles;
            if ignored_gravity_terms == IgnoreGravityTerms::WHFastOne {
                n_particles = 2;
            } else {
                n_particles = self.n_particles;
            }
            if let Some((star, particles)) = self.particles[..n_particles].split_first() {
                if let Some((star_newtonian_acceleration, particles_newtonian_accelerations)) = newtonian_inertial_accelerations[..n_particles].split_first_mut() {
                    for (a_newton_particle, particle_newtonian_acceleration) in particles_newtonian_accelerations.iter_mut().zip(particles.iter()) {
                        let dx = star.inertial_position.x - particle_newtonian_acceleration.position.x;
                        let dy = star.inertial_position.y - particle_newtonian_acceleration.position.y;
                        let dz = star.inertial_position.z - particle_newtonian_acceleration.position.z;
                        let r2 = dx.powi(2) + dy.powi(2) + dz.powi(2);
                        //println!("r2 {:e} {:e} {:e} {:e}", r2, dx.powi(2), dy.powi(2), dz.powi(2));
                        let r = r2.sqrt();
                        let prefac = G/(r2*r);
                        //println!("prefac {:e} = {:e}/({:e}*{:e})", prefac, G, r2, r);
                        let prefac_mass_star = prefac*star.mass;
                        let prefac_mass_particle = prefac*particle_newtonian_acceleration.mass;
                        //println!("-- mass_a {:e} | mass_b {:e}", star.mass, particle_newtonian_acceleration.mass);
                        //println!("* berfore: {:e}", star_newtonian_acceleration.x);
                        star_newtonian_acceleration.x -= prefac_mass_particle*dx;
                        //println!("* after: {:e}", star_newtonian_acceleration.x);
                        star_newtonian_acceleration.y -= prefac_mass_particle*dy;
                        star_newtonian_acceleration.z -= prefac_mass_particle*dz;
                        //println!("+ berfore: {:e}", a_newton_particle.x);
                        a_newton_particle.x += prefac_mass_star*dx;
                        //println!("+ after: {:e}", a_newton_particle.x);
                        a_newton_particle.y += prefac_mass_star*dy;
                        a_newton_particle.z += prefac_mass_star*dz;
                    }
                }
            }
        }
        return newtonian_inertial_accelerations;
    }

    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
    // [start] General Relativity FULL based on REBOUNDx gr.c
    fn calculate_newhall1983_general_relativity_acceleration(&mut self, ignored_gravity_terms: IgnoreGravityTerms) {
        let mut a_const = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; // array that stores the value of the constant term
        let mut a_new = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; // stores the newly calculated term
        let mut rs = [[0.; MAX_PARTICLES]; MAX_PARTICLES];
        let mut drs = [[Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; MAX_PARTICLES];
        
        let a_newton = self.get_anderson1975_newhall1983_newtonian_inertial_accelerations(ignored_gravity_terms);

        for (i, ((drs_i, rs_i), particle_i)) in 
                                        drs[..self.n_particles].iter_mut()
                                        .zip(rs[..self.n_particles].iter_mut())
                                        .zip(self.particles[..self.n_particles].iter())
                                        .enumerate() {
            // compute distances
            for (j, (drs_i_j, particle_j)) in drs_i[..self.n_particles].iter_mut()
                                            .zip(self.particles[..self.n_particles].iter())
                                            .enumerate() {
                if j != i{
                    drs_i_j.x = particle_i.inertial_position.x - particle_j.inertial_position.x;
                    drs_i_j.y = particle_i.inertial_position.y - particle_j.inertial_position.y;
                    drs_i_j.z = particle_i.inertial_position.z - particle_j.inertial_position.z;
                    rs_i[j] = (drs_i_j.x.powi(2) + drs_i_j.y.powi(2) + drs_i_j.z.powi(2)).sqrt();
                    //println!("i j: {} {} {:e}", i, j, rs_i[j]);
                }
            }
        }

        for (i, (((a_const_i, drs_i), rs_i), particle_i)) in a_const[..self.n_particles].iter_mut()
                                                .zip(drs[..self.n_particles].iter())
                                                .zip(rs[..self.n_particles].iter())
                                                .zip(self.particles[..self.n_particles].iter())
                                                .enumerate() {

            // then compute the constant terms:
            let mut a_constx = 0.;
            let mut a_consty = 0.;
            let mut a_constz = 0.;
            // 1st constant part
            for (j, (drs_i_j, particle_j)) in drs_i[..self.n_particles].iter()
                                            .zip(self.particles[..self.n_particles].iter())
                                            .enumerate() {
                if j != i {
                    let dxij = drs_i_j.x;
                    let dyij = drs_i_j.y;
                    let dzij = drs_i_j.z;
                    let rij2 = rs_i[j].powi(2);
                    let rij3 = rij2*rs_i[j];

                    let mut a1 = 0.;
                    for (k, (particle_k, rs_i_k)) in self.particles[..self.n_particles].iter()
                                                    .zip(rs_i[..self.n_particles].iter())
                                                        .enumerate() {
                        if k != i {
                            a1 += (4./(SPEED_OF_LIGHT_2)) * G*particle_k.mass/rs_i_k;
                        }
                    }

                    let mut a2 = 0.;
                    for (l, (particle_l, rs_l)) in self.particles[..self.n_particles].iter()
                                                    .zip(rs[..self.n_particles].iter())
                                                        .enumerate() {
                        if l != j {
                            a2 += (1./(SPEED_OF_LIGHT_2)) * G*particle_l.mass/rs_l[j];
                        }
                    }

                    let vi2= particle_i.inertial_velocity.x.powi(2) + particle_i.inertial_velocity.y.powi(2) + particle_i.inertial_velocity.z.powi(2);
                    let a3 = -vi2/SPEED_OF_LIGHT_2;
                    
                    let vj2 = particle_j.inertial_velocity.x.powi(2) + particle_j.inertial_velocity.y.powi(2) + particle_j.inertial_velocity.z.powi(2);
                    let a4 = -2.*vj2/SPEED_OF_LIGHT_2;

                    let a5 = (4./SPEED_OF_LIGHT_2) * (particle_i.inertial_velocity.x*particle_j.inertial_velocity.x + particle_i.inertial_velocity.y*particle_j.inertial_velocity.y + particle_i.inertial_velocity.z*particle_j.inertial_velocity.z);

                    let a6_0 = dxij*particle_j.inertial_velocity.x + dyij*particle_j.inertial_velocity.y + dzij*particle_j.inertial_velocity.z;
                    let a6 = (3./(2.*SPEED_OF_LIGHT_2)) * a6_0.powi(2)/rij2;

                    let factor1 = a1 + a2 + a3 + a4 + a5 + a6;
                    //println!("factors {:e} {:e} {:e} {:e} {:e} {:e}", a1, a2, a3, a4, a5, a6);
                    a_constx += G*particle_j.mass*dxij*factor1/rij3;
                    a_consty += G*particle_j.mass*dyij*factor1/rij3;
                    a_constz += G*particle_j.mass*dzij*factor1/rij3;

                    // 2nd constant part
                    let dvxij = particle_i.inertial_velocity.x - particle_j.inertial_velocity.x; 
                    let dvyij = particle_i.inertial_velocity.y - particle_j.inertial_velocity.y; 
                    let dvzij = particle_i.inertial_velocity.z - particle_j.inertial_velocity.z; 

                    let factor2 = dxij*(4.*particle_i.inertial_velocity.x - 3.*particle_j.inertial_velocity.x) + dyij*(4.*particle_i.inertial_velocity.y -3.*particle_j.inertial_velocity.y) + dzij*(4.*particle_i.inertial_velocity.z - 3.*particle_j.inertial_velocity.z);

                    a_constx += G*particle_j.mass*factor2*dvxij/rij3/SPEED_OF_LIGHT_2;
                    a_consty += G*particle_j.mass*factor2*dvyij/rij3/SPEED_OF_LIGHT_2;
                    a_constz += G*particle_j.mass*factor2*dvzij/rij3/SPEED_OF_LIGHT_2;
                }
            }
            a_const_i.x = a_constx;
            a_const_i.y = a_consty;
            a_const_i.z = a_constz;
            //println!("a_const_i {:?}", a_const_i);
        }


        let dev_limit = 1.0e-30;
        let max_iterations = 10;
        // Now running the substitution again and again through the loop below
        for k in 0..max_iterations {
            let a_old = a_new.clone();
            // now add on the non-constant term
            for (i, (((a_new_i, drs_i), rs_i), a_const_i)) in a_new[..self.n_particles].iter_mut()
                                            .zip(drs[..self.n_particles].iter())
                                            .zip(rs[..self.n_particles].iter())
                                            .zip(a_const[..self.n_particles].iter())
                                            .enumerate() {
                let mut non_constx = 0.;
                let mut non_consty = 0.;
                let mut non_constz = 0.;
                for (j, ((((a_old_j, a_newton_j), drs_i_j), rs_i_j), particle_j)) in a_old[..self.n_particles].iter()
                                                                                            .zip(a_newton[..self.n_particles].iter())
                                                                                            .zip(drs_i[..self.n_particles].iter())
                                                                                            .zip(rs_i[..self.n_particles].iter())
                                                                                            .zip(self.particles[..self.n_particles].iter())
                                                                                            .enumerate() {
                    if j != i {
                        let dxij = drs_i_j.x;
                        let dyij = drs_i_j.y;
                        let dzij = drs_i_j.z;
                        let rij = rs_i_j;
                        let rij2 = rij.powi(2);
                        let rij3 = rij2*rij;
                        non_constx += (G*particle_j.mass*dxij/rij3)*(dxij*(a_newton_j.x+a_old_j.x)+dyij*(a_newton_j.y+a_old_j.y)+
                                    dzij*(a_newton_j.z+a_old_j.z))/(2.*SPEED_OF_LIGHT_2) + (7./(2.*SPEED_OF_LIGHT_2))*G*particle_j.mass*(a_newton_j.x+a_old_j.x)/rij;
                        non_consty += (G*particle_j.mass*dyij/rij3)*(dxij*(a_newton_j.x+a_old_j.x)+dyij*(a_newton_j.y+a_old_j.y)+
                                    dzij*(a_newton_j.z+a_old_j.z))/(2.*SPEED_OF_LIGHT_2) + (7./(2.*SPEED_OF_LIGHT_2))*G*particle_j.mass*(a_newton_j.y+a_old_j.y)/rij;
                        non_constz += (G*particle_j.mass*dzij/rij3)*(dxij*(a_newton_j.x+a_old_j.x)+dyij*(a_newton_j.y+a_old_j.y)+
                                    dzij*(a_newton_j.z+a_old_j.z))/(2.*SPEED_OF_LIGHT_2) + (7./(2.*SPEED_OF_LIGHT_2))*G*particle_j.mass*(a_newton_j.z+a_old_j.z)/rij;
                    }
                }
                a_new_i.x = a_const_i.x + non_constx;
                a_new_i.y = a_const_i.y + non_consty;
                a_new_i.z = a_const_i.z + non_constz;
                //println!("non_constx {:?}", non_constx);
                //println!("non_consty {:?}", non_consty);
                //println!("non_constz {:?}", non_constz);
            }
            
            // break out loop if a_new is converging
            let mut maxdev = 0.;
            let mut dx = 0.;
            let mut dy = 0.;
            let mut dz = 0.;
            for (a_new_i, a_old_i) in a_new[..self.n_particles].iter()
                                    .zip(a_old[..self.n_particles].iter()) {
                if a_new_i.x.abs() < dev_limit {
                    dx = (a_new_i.x - a_old_i.x).abs() / a_new_i.x;
                }
                if a_new_i.y.abs() < dev_limit {
                    dy = (a_new_i.y - a_old_i.y).abs() / a_new_i.y;
                }
                if a_new_i.z.abs() < dev_limit {
                    dz = (a_new_i.z - a_old_i.z).abs() / a_new_i.z;
                }
                if dx > maxdev { maxdev = dx; }
                if dy > maxdev { maxdev = dy; }
                if dz > maxdev { maxdev = dz; }

            }

            if maxdev < dev_limit {
                break;
            } else if k == max_iterations {
                println!("[WARNING {} UTC] {} iterations in general relativity failed to converge.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), max_iterations);
            }

        }
        
        // update acceleration in particles
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            if let Some((a_new_star, a_new_particles)) = a_new[..self.n_particles].split_first_mut() {
                // This algorithm computes general_relativity_acceleration in the inertial frame
                // we need to convert it to the heliocentric frame to be consistent with the rest of
                // additional forces by subtracting star.inertial_acceleration
                //
                // inertial_position = heliocentric_position - (1/total_mass) * sum(particle.mass * particle.position)
                // inertial_velocity = heliocentric_velocity - (1/total_mass) * sum(particle.mass * particle.velocity)
                // inertial_acceleration = heliocentric_acceleration - (1/total_mass) * sum(particle.mass * particle.acceleration)
                //
                // heliocentric_position = inertial_position + (1/star.mass) * sum(particle.mass * particle.position)
                // heliocentric_velocity = inertial_velocity + (1/star.mass) * sum(particle.mass * particle.velocity)
                // heliocentric_acceleration = inertial_acceleration + (1/star.mass) * sum(particle.mass * particle.acceleration)
                //
                // 1.- Compute second term using only general relativity accelerations
                let mut reference_gr_acceleration = Axes{x:0., y:0., z:0.};
                reference_gr_acceleration.x  = a_new_star.x*star.mass;
                reference_gr_acceleration.y  = a_new_star.y*star.mass;
                reference_gr_acceleration.z  = a_new_star.z*star.mass;
                for (particle, particle_acceleration) in particles.iter().zip(a_new_particles.iter()) {
                    reference_gr_acceleration.x += particle_acceleration.x*particle.mass;
                    reference_gr_acceleration.y += particle_acceleration.y*particle.mass;
                    reference_gr_acceleration.z += particle_acceleration.z*particle.mass;
                }
                reference_gr_acceleration.x /= star.mass;
                reference_gr_acceleration.y /= star.mass;
                reference_gr_acceleration.z /= star.mass;
                // 2.- Convert from inertial to heliocentric
                for (particle, a_new_particle) in particles.iter_mut().zip(a_new_particles.iter()){
                    particle.general_relativity_acceleration.x = a_new_particle.x + reference_gr_acceleration.x;
                    particle.general_relativity_acceleration.y = a_new_particle.y + reference_gr_acceleration.y;
                    particle.general_relativity_acceleration.z = a_new_particle.z + reference_gr_acceleration.z;
                }
                star.general_relativity_acceleration.x = a_new_star.x + reference_gr_acceleration.x;
                star.general_relativity_acceleration.y = a_new_star.y + reference_gr_acceleration.y;
                star.general_relativity_acceleration.z = a_new_star.z + reference_gr_acceleration.z;
            }
        }

    }
    // [end] General Relativity FULL based on REBOUNDx gr.c
    //--------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////

}
