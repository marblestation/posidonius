use std::collections::HashMap;
use super::super::constants::{K2, G, R_SUN, SUN_DYN_FREQ, SPEED_OF_LIGHT_2, MAX_PARTICLES, MAX_DISTANCE_2, K_WIND, WSAT_WIND, WSAT_WIND_2};
use super::{Evolver, EvolutionType, SolarEvolutionType};
use super::{Particle};
use super::{Axes};

#[derive(Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub struct Universe {
    pub initial_time: f64,
    pub time_limit: f64,
    pub particles: [Particle; MAX_PARTICLES],
    pub particles_evolvers: Vec<Evolver>,
    pub n_particles: usize,
    pub consider_tides: bool,
    pub consider_rotational_flattening: bool,
    pub consider_general_relativy: bool,
    star_planet_dependent_dissipation_factors : HashMap<usize, f64>, // Central body specific
    temporary_copied_particle_positions: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particle_velocities: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_masses: [f64; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_radiuses: [f64; MAX_PARTICLES], // For optimization purposes
}

impl Universe {
    pub fn new(mut particles: Vec<Particle>, initial_time: f64, time_limit: f64, 
              consider_tides: bool, consider_rotational_flattening:
              bool, consider_general_relativy: bool) -> Universe {
        if consider_general_relativy {
            let star_index = 0; // index
            let local_copy_star_mass_g = particles[star_index].mass_g;
            for particle in particles[1..].iter_mut() {
                particle.general_relativity_factor =  local_copy_star_mass_g*particle.mass_g / (local_copy_star_mass_g + particle.mass_g).powi(2)
            }
        }

        let temporary_copied_particle_positions = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES];
        let temporary_copied_particle_velocities = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES];
        let temporary_copied_particles_masses = [0.; MAX_PARTICLES];
        let temporary_copied_particles_radiuses = [0.; MAX_PARTICLES];
    
        // OPTIMIZATION: Transform vector to array
        // - Arrays are stored in the stack which is faster than the heap (where vectors are allocated)
        // - The array should have a fixed size, thus it should always be equal or greater to the vector passed
        // - The array elements to be considered will be limited by the n_particles value
        let n_particles = particles.len();
        if n_particles > MAX_PARTICLES {
            panic!("Only {} bodies are allowed, you need to increase the MAX_PARTICLE constant.", MAX_PARTICLES);
        }
        let mut transformed_particles = [Particle::new_dummy(); MAX_PARTICLES];
        let mut particles_evolvers : Vec<Evolver> = Vec::with_capacity(n_particles);
        for i in 0..n_particles {
            transformed_particles[i] = particles[i];
            transformed_particles[i].id = i;
            particles_evolvers.push(Evolver::new(transformed_particles[i].evolution_type, initial_time, time_limit));
        }

        let mut universe = Universe {
                    initial_time: initial_time,
                    time_limit: time_limit,
                    particles: transformed_particles,
                    particles_evolvers: particles_evolvers,
                    n_particles: n_particles,
                    consider_tides: consider_tides,
                    consider_rotational_flattening: consider_rotational_flattening,
                    consider_general_relativy: consider_general_relativy,
                    star_planet_dependent_dissipation_factors:HashMap::new(),
                    temporary_copied_particle_positions: temporary_copied_particle_positions,
                    temporary_copied_particle_velocities: temporary_copied_particle_velocities,
                    temporary_copied_particles_masses: temporary_copied_particles_masses,
                    temporary_copied_particles_radiuses: temporary_copied_particles_radiuses,
                    };
        let time_step = 0.;
        let current_time = 0.;
        universe.evolve_particles(current_time, time_step); // Make sure we start with the good initial values
        universe
    }

    pub fn gravity_calculate_acceleration(&mut self, integrator_is_whfasthelio: bool) {

        for (i, particle) in self.particles[..self.n_particles].iter().enumerate() {
            self.temporary_copied_particle_positions[i].x = particle.position.x;
            self.temporary_copied_particle_positions[i].y = particle.position.y;
            self.temporary_copied_particle_positions[i].z = particle.position.z;
            self.temporary_copied_particle_velocities[i].x = particle.velocity.x;
            self.temporary_copied_particle_velocities[i].y = particle.velocity.y;
            self.temporary_copied_particle_velocities[i].z = particle.velocity.z;
            self.temporary_copied_particles_masses[i] = particle.mass;
            self.temporary_copied_particles_radiuses[i] = particle.radius;
        }

        for (i, particle_a) in self.particles[..self.n_particles].iter_mut().enumerate() {
			particle_a.acceleration.x = 0.;
			particle_a.acceleration.y = 0.;
			particle_a.acceleration.z = 0.;
            for (j, (((particle_b_position, particle_b_velocity), particle_b_radius), particle_b_mass)) in self.temporary_copied_particle_positions[..self.n_particles].iter()
                                                                .zip(self.temporary_copied_particle_velocities[..self.n_particles].iter())
                                                                .zip(self.temporary_copied_particles_radiuses[..self.n_particles].iter())
                                                                .zip(self.temporary_copied_particles_masses[..self.n_particles].iter()).enumerate() {
                if i == j {
                    continue;
                }

                //// Check collisions and ejections //////////////////////////////////
                let dx = particle_a.position.x - particle_b_position.x;
                let dy = particle_a.position.y - particle_b_position.y;
                let dz = particle_a.position.z - particle_b_position.z;
                let distance_2 = dx*dx + dy*dy + dz*dz;
                // Do not check twice the same pair of particles:
                if i < j {
                    // Check if particles are overlapping
                    if distance_2 < (particle_a.radius + particle_b_radius).powi(2) {
                        let dvx = particle_a.velocity.x - particle_b_velocity.x;
                        let dvy = particle_a.velocity.y - particle_b_velocity.y;
                        let dvz = particle_a.velocity.z - particle_b_velocity.z;
                        // Check if particles are approaching each other by checking the sign of the non-normalized radial velocity
                        if dvx*dx + dvy*dy + dvz*dz < 0. {
                            panic!("Collision between particle {} and {}!", i, j);
                        }
                    }
                    if i == 0 && distance_2 > MAX_DISTANCE_2 {
                        panic!("Particle {} has been ejected!", j);
                    }
                }
                //////////////////////////////////////////////////////////////////////

                if integrator_is_whfasthelio && (i == 0 || j == 0) {
                    // For WHFastHelio, ignore central body
                    continue
                }

                // gravitational force equation for vector quantities
                // Fx = - G * ((m1 * m2) / r^3) * (x2 - x1)
                // Fy = - G * ((m1 * m2) / r^3) * (y2 - y1)
                // Fz = - G * ((m1 * m2) / r^3) * (z2 - z1)
                //
                // acceleration:
                // ax = Fx / m1 = - G * (m2 / r^3) * (x2 - x1)
                //let dx = particle_a.position.x - particle_b_position.x;
                //let dy = particle_a.position.y - particle_b_position.y;
                //let dz = particle_a.position.z - particle_b_position.z;
                //let r = (dx*dx + dy*dy + dz*dz).sqrt();
                let r = distance_2.sqrt();
                let prefact = -G/(r*r*r) * particle_b_mass;

                particle_a.acceleration.x += prefact * dx;
                particle_a.acceleration.y += prefact * dy;
                particle_a.acceleration.z += prefact * dz;
            }
        }

        if !integrator_is_whfasthelio {
            // Tides require heliocentric point of reference, the star should continue in the zero point
            // so we must compensate all the planets (but if WHFastHelio is being used, it is
            // automatically done by the integrator):
            if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
                for particle in particles.iter_mut() {
                    particle.acceleration.x -= star.acceleration.x;
                    particle.acceleration.y -= star.acceleration.y;
                    particle.acceleration.z -= star.acceleration.z;
                }
                star.acceleration.x = 0.;
                star.acceleration.y = 0.;
                star.acceleration.z = 0.;
            }
        }
    }

    pub fn calculate_additional_forces(&mut self, current_time: f64, time_step: f64, only_dspin_dt: bool) {
        self.evolve_particles(current_time, time_step);

        // Compute this universe
        if self.consider_tides || self.consider_rotational_flattening || self.consider_general_relativy  {
            self.calculate_distance_and_velocities(); // Needed for tides, rotational flattening and general relativity
            if self.consider_tides || self.consider_rotational_flattening  {
                self.calculate_torques(); // Needed for dspin_dt
                self.calculate_dspin_dt(); // Needed for tides and rotational flattening
            }
        }
        if !only_dspin_dt {
            self.calculate_acceleration_corrections();
            self.apply_acceleration_corrections();
        }

        // Recover first original body as frame of reference
        self.center_to_first_particle();
        //for (i, particle) in self.particles[1..self.n_particles].iter().enumerate() {
            //println!("{} - Acceleration {:e} {:e} {:e}", i, particle.acceleration.x, particle.acceleration.y, particle.acceleration.z);
        //}
    }


    pub fn center_to_first_particle(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.position.x  -= star.position.x;
                particle.position.y  -= star.position.y;
                particle.position.z  -= star.position.z;
                particle.velocity.x -= star.velocity.x;
                particle.velocity.y -= star.velocity.y;
                particle.velocity.z -= star.velocity.z;
                particle.acceleration.x -= star.acceleration.x;
                particle.acceleration.y -= star.acceleration.y;
                particle.acceleration.z -= star.acceleration.z;
            }
            star.position.x = 0.;
            star.position.y = 0.;
            star.position.z = 0.;
            star.velocity.x = 0.;
            star.velocity.y = 0.;
            star.velocity.z = 0.;
            star.acceleration.x = 0.;
            star.acceleration.y = 0.;
            star.acceleration.z = 0.;
            star.torque.x = 0.;
            star.torque.y = 0.;
            star.torque.z = 0.;
        }
    }

    fn calculate_torques(&mut self) {
        let central_body = true;

        self.calculate_planet_dependent_dissipation_factors(); // Needed by calculate_orthogonal_component_of_the_tidal_force and calculate_orthogonal_component_of_the_tidal_force if MathisSolarLike
        self.calculate_scalar_product_of_vector_position_with_spin(); // Needed by tides and rotational flattening
        
        if self.consider_tides {
            self.calculate_orthogonal_component_of_the_tidal_force(central_body);
            self.calculate_orthogonal_component_of_the_tidal_force(!central_body);
            self.calculate_torque_due_to_tides(central_body);
            self.calculate_torque_due_to_tides(!central_body);
        }

        if self.consider_rotational_flattening {
            self.calculate_radial_and_orthogonal_components_of_the_force_induced_by_rotational_flattening(); // radial component used later on by calculate_acceleration_induced_by_rotational_flattering
            self.calculate_torque_induced_by_rotational_flattening(central_body);
            self.calculate_torque_induced_by_rotational_flattening(!central_body);
        }
    }

    fn calculate_acceleration_corrections(&mut self) {
        if self.consider_tides {
            self.calculate_radial_component_of_the_tidal_force();  // Needed for calculate_tidal_acceleration
            self.calculate_tidal_acceleration();
        }

        if self.consider_rotational_flattening {
            self.calculate_acceleration_induced_by_rotational_flattering();
        }

        if self.consider_general_relativy {
            self.calculate_general_relativity_acceleration();
        }
    }

    fn apply_acceleration_corrections(&mut self) {
        // Add the tidal+flattening+general relativity accelerations to the gravitational one (already computed)
        //for particle in self.particles[1..self.n_particles].iter_mut() {
        for particle in self.particles[..self.n_particles].iter_mut() {
            if self.consider_tides {
                particle.acceleration.x += particle.tidal_acceleration.x;
                particle.acceleration.y += particle.tidal_acceleration.y;
                particle.acceleration.z += particle.tidal_acceleration.z;
                //println!("Tides acceleration {:e} {:e} {:e}", particle.tidal_acceleration.x, particle.tidal_acceleration.y, particle.tidal_acceleration.z);
            }

            if self.consider_rotational_flattening {
                particle.acceleration.x += particle.acceleration_induced_by_rotational_flattering.x;
                particle.acceleration.y += particle.acceleration_induced_by_rotational_flattering.y;
                particle.acceleration.z += particle.acceleration_induced_by_rotational_flattering.z;
                //println!("Rot acceleration {:e} {:e} {:e}", particle.acceleration_induced_by_rotational_flattering.x, particle.acceleration_induced_by_rotational_flattering.y, particle.acceleration_induced_by_rotational_flattering.z);
            }

            if self.consider_general_relativy {
                particle.acceleration.x += particle.general_relativity_acceleration.x;
                particle.acceleration.y += particle.general_relativity_acceleration.y;
                particle.acceleration.z += particle.general_relativity_acceleration.z;
                //println!("GR acceleration {:e} {:e} {:e}", particle.general_relativity_acceleration.x, particle.general_relativity_acceleration.y, particle.general_relativity_acceleration.z);
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
                    let factor = K2 / (star.mass_g + particle.mass_g);
                    //let factor = 1. / (star.mass + particle.mass);
                    torque.x += factor * n_tid_x;
                    torque.y += factor * n_tid_y;
                    torque.z += factor * n_tid_z;
                } else {
                    particle.torque.x = n_tid_x;
                    particle.torque.y = n_tid_y;
                    particle.torque.z = n_tid_z;
                }

            }

            if central_body {
                // - Equation 25 from Bolmont et al. 2015
                star.torque.x = torque.x;
                star.torque.y = torque.y;
                star.torque.z = torque.z;
            }
        }

    }


    fn calculate_dspin_dt(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            // - Equation 25 from Bolmont et al. 2015
            let factor = - 1. / (star.radius_of_gyration_2 * star.radius.powi(2));
            star.dspin_dt.x = factor * star.torque.x;
            star.dspin_dt.y = factor * star.torque.y;
            star.dspin_dt.z = factor * star.torque.z;

            for particle in particles.iter_mut() {
                let factor = - K2 * star.mass_g / (particle.mass_g * (particle.mass_g + star.mass_g) 
                                * particle.radius_of_gyration_2 * particle.radius.powi(2));
                // - Equation 25 from Bolmont et al. 2015
                particle.dspin_dt.x = factor * particle.torque.x;
                particle.dspin_dt.y = factor * particle.torque.y;
                particle.dspin_dt.z = factor * particle.torque.z;
            }
        }
    }

    fn calculate_orthogonal_component_of_the_tidal_force(&mut self, central_body:bool) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                // (distance to star)^7
                let distance_7 = particle.distance.powi(7);

                //// Tidal force calculation (star) :: Only orthogonal component is needed
                // - Equation 5 from Bolmont et al. 2015
                // - Ftides in Msun.AU.day-1
                if central_body {
                    // - F_tides_ortho_star
                    let star_dissipation_factor = Universe::planet_dependent_dissipation_factor(&self.star_planet_dependent_dissipation_factors, &star.id, star.evolution_type, star.scaled_dissipation_factor);
                    particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass_g.powi(2))
                                                    * (star.radius.powi(10)) 
                                                    * star_dissipation_factor / ( (K2.powi(2)) * distance_7);
                    //particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass.powi(2))
                                                        //* (star.radius.powi(10)) 
                                                        //* star.scaled_dissipation_factor / (distance_7);
                } else {
                    // - F_tides_ortho_plan
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
        // dEdt_tides
        // F_tides_rad
        // F_tides_rad_cons
        // F_tides_rad_diss
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let star_mass_2 = star.mass_g * star.mass_g;

            for particle in particles.iter_mut() {
                let planet_mass_2 = particle.mass_g * particle.mass_g;
                // Conservative part of the radial tidal force
                // - Ftidr_cons

                let radial_component_of_the_tidal_force_conservative_part = -3.0 / (particle.distance.powi(7) * K2)
                            * (planet_mass_2 * star.radius.powi(5) * star.love_number 
                            + star_mass_2 * particle.radius.powi(5) * particle.love_number);
                //let radial_component_of_the_tidal_force_conservative_part = -3.0 / particle.distance.powi(7)
                            //* (planet_mass_2 * star.radius.powi(5) * star.love_number 
                            //+ star_mass_2 * particle.radius.powi(5) * particle.love_number);

                // Dissipative part of the radial tidal force
                // - Ftidr_diss
                let factor1 = -13.5 * particle.radial_velocity / (particle.distance.powi(8) * K2*K2);
                //let factor1 = -13.5 * particle.radial_velocity / particle.distance.powi(8);
                let star_dissipation_factor = Universe::planet_dependent_dissipation_factor(&self.star_planet_dependent_dissipation_factors, &star.id, star.evolution_type, star.scaled_dissipation_factor);
                let term1 = planet_mass_2
                            * star.radius.powi(10)
                            * star_dissipation_factor;
                let term2 = star_mass_2
                            * particle.radius.powi(10)
                            * particle.scaled_dissipation_factor;
                let radial_component_of_the_tidal_force_dissipative_part = factor1 * (term1 + term2 );

                // If we consider the star as a point mass (used for denergy_dt calculation):
                let radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass = factor1 * term2;

                // Sum of the dissipative and conservative part of the radial force
                particle.radial_component_of_the_tidal_force = radial_component_of_the_tidal_force_conservative_part + radial_component_of_the_tidal_force_dissipative_part;

                ////// Instantaneous energy loss dE/dt due to tides
                ////// in Msun.AU^2.day^(-3)
                //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
                let factor2 = particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance;
                particle.denergy_dt = -((1.0 / particle.distance * (radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.radial_velocity))
                            * (particle.position.x*particle.velocity.x + particle.position.y*particle.velocity.y + particle.position.z*particle.velocity.z)
                            + factor2 
                            * ((particle.spin.y*particle.position.z - particle.spin.z*particle.position.y - particle.velocity.x) * particle.velocity.x
                            + (particle.spin.z*particle.position.x - particle.spin.x*particle.position.z - particle.velocity.y) * particle.velocity.y
                            + (particle.spin.x*particle.position.y - particle.spin.y*particle.position.x - particle.velocity.z) * particle.velocity.z))
                            + star.mass_g/(star.mass_g + particle.mass_g) 
                            * (particle.torque.x*particle.spin.x + particle.torque.y*particle.spin.y + particle.torque.z*particle.spin.z) ;
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

    fn calculate_general_relativity_acceleration(&mut self) {

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
    fn calculate_radial_and_orthogonal_components_of_the_force_induced_by_rotational_flattening(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            // Calculation of the norm square of the spin for the planet
            let normspin_2_star = star.spin.x.powi(2) + star.spin.y.powi(2) + star.spin.z.powi(2);
            for particle in particles.iter_mut() {
                //// Star:
                // - Equation 16 from Bolmont et al. 2015
                let csi = particle.mass_g * star.fluid_love_number * normspin_2_star * star.radius.powi(5) / (6. * K2); // Msun.AU^5.day-2
                // - Equation 16 from Bolmont et al. 2015
                particle.orthogonal_component_of_the_force_induced_by_star_rotation = -6. * csi * particle.scalar_product_of_vector_position_with_stellar_spin / (normspin_2_star * particle.distance.powi(5));
                
                //// Planet:
                // Calculation of the norm square of the spin for the planet
                let normspin_2_planet = particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2);
                // - Equation 16 from Bolmont et al. 2015
                let cpi = star.mass_g * particle.fluid_love_number * normspin_2_planet * particle.radius.powi(5) / (6. * K2); // Msun.AU^5.day-2
            
                // - Equation 16 from Bolmont et al. 2015
                particle.orthogonal_component_of_the_force_induced_by_planet_rotation = -6. * cpi * particle.scalar_product_of_vector_position_with_planetary_spin / (normspin_2_planet * particle.distance.powi(5));
                
                ////
                // - Equation 15 from Bolmont et al. 2015
                particle.radial_component_of_the_force_induced_by_rotation = -3./particle.distance.powi(5) * (cpi + csi)
                    + 15./particle.distance.powi(7) * (csi * particle.scalar_product_of_vector_position_with_stellar_spin * particle.scalar_product_of_vector_position_with_stellar_spin/normspin_2_star
                        + cpi * particle.scalar_product_of_vector_position_with_planetary_spin * particle.scalar_product_of_vector_position_with_planetary_spin/normspin_2_planet); // Msun.AU.day-1
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
                    // Integration of the spin (total torque tides):
                    let factor = K2 / (star.mass_g + particle.mass_g);
                    //let factor = 1. / (star.mass + particle.mass);
                    torque.x += factor * n_rot_x;
                    torque.y += factor * n_rot_y;
                    torque.z += factor * n_rot_z;
                } else {
                    particle.torque.x += n_rot_x; // Add to the torque due to tides (if computed)
                    particle.torque.y += n_rot_y;
                    particle.torque.z += n_rot_z;
                }
            }

            if central_body {
                // - Equation 25 from Bolmont et al. 2015
                star.torque.x = torque.x;
                star.torque.y = torque.y;
                star.torque.z = torque.z;
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


    fn calculate_distance_and_velocities(&mut self) {
        // Calculation of velocity vv(j), radial velocity vrad(j)
        // velocities in AU/day
        for particle in self.particles[1..self.n_particles].iter_mut() {
            // Norm of the velocity
            let v_2 = (particle.velocity.x.powi(2)) 
                                + (particle.velocity.y.powi(2))
                                + (particle.velocity.z.powi(2));
            let norm_vel = v_2.sqrt();
            particle.norm_velocity_vector_2 = v_2;
            particle.norm_velocity_vector = norm_vel;

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
            EvolutionType::SolarLike(model) => {
                match model {
                    SolarEvolutionType::EvolvingDissipation(_) => {
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
            },
            _ => {},
        }

    }

    pub fn planet_dependent_dissipation_factor(star_planet_dependent_dissipation_factors: &HashMap<usize, f64>,  id: &usize, evolution_type: EvolutionType, scaled_dissipation_factor: f64) -> f64 {
        match evolution_type {
            EvolutionType::SolarLike(model) => {
                match model {
                    SolarEvolutionType::EvolvingDissipation(_) => {
                        match star_planet_dependent_dissipation_factors.get(id) {
                            Some(&value) => value,
                            _ => scaled_dissipation_factor // This should not happen
                        }
                    },
                    _ => scaled_dissipation_factor,
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
    
    pub fn evolve_particles(&mut self, current_time: f64, time_step: f64) {
        for (particle, evolver) in self.particles[..self.n_particles].iter_mut().zip(self.particles_evolvers.iter_mut()) {
            ////////////////////////////////////////////////////////////////////
            // Wind
            ////////////////////////////////////////////////////////////////////
            // - It requires the radius before being evolved
            // - If time_step is zero, it's the first initialization and the wind factor will be zero
            particle.wind_factor = match evolver.evolution_type {
                EvolutionType::SolarLike(_) => { 
                    let threshold = (particle.spin.x.powi(2) + particle.spin.y.powi(2) +
                                     particle.spin.z.powi(2)).sqrt();
                    if threshold >= WSAT_WIND {
                        // Friendly reminder that m(1) is in solar mass * K2
                        //tmp2 = hdt * K_wind * wsat_wind*wsat_wind * sqrt(Rsth*K2/(Rsun*m(1)))
                        let old_radius = particle.radius; // TODO
                        let factor = time_step * K_WIND * WSAT_WIND_2 * (old_radius/R_SUN * 1./particle.mass).sqrt();
                        Axes{ 
                            x: factor * particle.spin.x,
                            y: factor * particle.spin.y,
                            z: factor * particle.spin.z,
                        }
                    } else {
                        //tmp2 = hdt * K_wind * sqrt(Rsth*K2/(Rsun*m(1)))
                        let old_radius = particle.radius; // TODO
                        let factor = time_step * K_WIND * (old_radius/R_SUN * 1./particle.mass).sqrt();
                        Axes{ 
                            x: factor * particle.spin.x.powi(3),
                            y: factor * particle.spin.y.powi(3),
                            z: factor * particle.spin.z.powi(3),
                        }
                    }
                },
                _ => Axes{x: 0., y:0., z: 0.},
            };

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
                    EvolutionType::SolarLike(model) => {
                        match model {
                            SolarEvolutionType::EvolvingDissipation(_) => {
                                let inverse_tidal_q_factor = evolver.inverse_tidal_q_factor(current_time, 0.);
                                //
                                // Calculation of the norm square of the spin for the star
                                let normspin_2 = particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2);
                                let epsilon_squared = normspin_2/SUN_DYN_FREQ;
                                // Normal formula = 3.d0*epsilon_squared*Q_str_inv/(4.d0*k2s)
                                // but as for sigma it is necessary to divide by k2s, we do not divide here
                                let lag_angle = 3.0*epsilon_squared*inverse_tidal_q_factor/4.0;
                                lag_angle
                            },
                            _ => 0.,
                        }
                    },
                    _ => 0.,
            };
            //println!("[{}] Evolve Radius {:e} Gyration {:e} Love {:e} Lag {:e}", current_time, particle.radius, particle.radius_of_gyration_2, particle.love_number, particle.lag_angle);
        }

    }

}
