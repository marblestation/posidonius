use std::collections::HashMap;
use super::super::constants::{K2, G, SUN_DYN_FREQ, SPEED_OF_LIGHT_2, MAX_PARTICLES, MAX_DISTANCE_2};
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
    pub consider_all_body_interactions: bool,
    star_planet_dependent_dissipation_factors : HashMap<usize, f64>, // Central body specific
    parallel_universes : Vec<Universe>, // If consider all the interactions
    temporary_copied_particle_positions: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particle_velocities: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_masses: [f64; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_radiuses: [f64; MAX_PARTICLES], // For optimization purposes
}

impl Universe {
    pub fn new(mut particles: Vec<Particle>, initial_time: f64, time_limit: f64, 
              consider_tides: bool, consider_rotational_flattening:
              bool, consider_general_relativy: bool, consider_all_body_interactions: bool) -> Universe {
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

        // Parallel universes
        let n_parallel_universes;
        if consider_all_body_interactions {
            n_parallel_universes = n_particles - 2;
        } else {
            n_parallel_universes = 0;
        }
        let parallel_universes : Vec<Universe> = Vec::with_capacity(n_parallel_universes);

        let mut universe = Universe {
                    initial_time: initial_time,
                    time_limit: time_limit,
                    particles: transformed_particles,
                    particles_evolvers: particles_evolvers,
                    n_particles: n_particles,
                    consider_tides: consider_tides,
                    consider_rotational_flattening: consider_rotational_flattening,
                    consider_general_relativy: consider_general_relativy,
                    consider_all_body_interactions: consider_all_body_interactions && n_parallel_universes > 0,
                    star_planet_dependent_dissipation_factors:HashMap::new(),
                    parallel_universes: parallel_universes,
                    temporary_copied_particle_positions: temporary_copied_particle_positions,
                    temporary_copied_particle_velocities: temporary_copied_particle_velocities,
                    temporary_copied_particles_masses: temporary_copied_particles_masses,
                    temporary_copied_particles_radiuses: temporary_copied_particles_radiuses,
                    };
        let current_time = 0.;
        universe.evolve_particles(current_time); // Make sure we start with the good initial values
        for i in 0..n_parallel_universes {
            let mut parallel_universe = universe.clone();
            parallel_universe.n_particles -= i+1;
            // The correct transfer all particles except central body will be executed when needed
            // later on
            universe.parallel_universes.push(parallel_universe);
        }
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

                //// Check collisions //////////////////////////////////////////////
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

    pub fn gravity_calculate_acceleration2(&mut self) {
        // https://www.cs.utoronto.ca/~wayne/research/thesis/msc/node5.html
        // force softening, i.e.,  replacing  tex2html_wrap_inline2549 with  tex2html_wrap_inline2551 in
        // the denominator of the gravitational force computation for some small constant
        // tex2html_wrap_inline2553 , usually chosen to approximate the average inter-particle separation.
        // This is done because it allows a smaller N to approximate a larger N, and also to eliminate the
        // singularity at r=0 [6].
        // [6] James Binney and Scott Tremaine. Galactic Dynamics. Princeton Series in Astrophysics.
        //      Princeton University Press, 1987.
        //
        //  https://en.wikipedia.org/wiki/N-body_problem#Few_bodies
        //  For a small number of bodies, an n-body problem can be solved using direct methods,
        //  also called particle-particle methods. These methods numerically integrate the
        //  differential equations of motion. Numerical integration for this problem can be a
        //  challenge for several reasons. First, the gravitational potential is singular; it goes
        //  to infinity as the distance between two particles goes to zero. The gravitational
        //  potential may be softened to remove the singularity at small distances
        //
        // TODO: Check rebound
    }

    pub fn calculate_additional_forces(&mut self, current_time: f64, only_dspin_dt: bool) {
        ////////////////////////////////////////////////////////////////////////
        // Parallel universes: Copies of the current universe but without the
        // first central body
        //  1.- First parallel universe removes first body and uses the second 
        //      body as point of reference (central body)
        //  2.- Second parallel universe removes first and second body, it uses 
        //      the third body as point of reference (central body)
        //  3.- Third parallel universe...
        ////////////////////////////////////////////////////////////////////////
        // Dependencies for Parallel Computation
        // NOTE: tests show that this parallel approach leads to slower results)
        //
        //extern crate threadpool;
        //use self::threadpool::ThreadPool;
        //use std::sync::mpsc;
        //let n_workers = 10;
        //let thread_pool: ThreadPool::new(n_workers),
        //let (tx, rx) = mpsc::channel();
        ////////////////////////////////////////////////////////////////////////

        self.evolve_particles(current_time);

        if self.consider_all_body_interactions {
            // Build parallel universes if needed
            for (i, parallel_universe) in self.parallel_universes.iter_mut().enumerate() {
                // Transfer all particles except central body
                for (particle, particle_parallel_universe) in self.particles[i+1..self.n_particles].iter()
                                                            .zip(parallel_universe.particles[..parallel_universe.n_particles].iter_mut()) {
                    particle_parallel_universe.position = particle.position;
                    particle_parallel_universe.velocity = particle.velocity;
                    // Forget accelerations from the original universe
                    particle_parallel_universe.acceleration.x = 0.0;
                    particle_parallel_universe.acceleration.y = 0.0;
                    particle_parallel_universe.acceleration.z = 0.0;
                    // Potentially evolving values
                    particle_parallel_universe.radius = particle.radius;
                    particle_parallel_universe.radius_of_gyration_2 = particle.radius_of_gyration_2;
                    particle_parallel_universe.love_number = particle.love_number;
                    particle_parallel_universe.lag_angle = particle.lag_angle;
                }
                // Re-center system over the new first body
                parallel_universe.center_to_first_particle();

                // Compute parallel universes
                parallel_universe.calculate_dspin_dt();
                if !only_dspin_dt {
                    parallel_universe.calculate_acceleration_corrections();
                }
            }

        }

        ////////////////////////////////////////////////////////////////////////
        ////// Execute in parallel
        //if self.consider_all_body_interactions {
            //for parallel_universe in self.parallel_universes.iter_mut() {
                //let tx = tx.clone();
                //let mut parallel_universe = parallel_universe.clone();

                //thread_pool.execute(move || {
                    //parallel_universe.particles[0].acceleration.x = 1.;
                    //parallel_universe.calculate_dspin_dt();
                    //if !only_dspin_dt {
                        //parallel_universe.calculate_acceleration_corrections();
                    //}
                    
                    //tx.send(parallel_universe).unwrap();
                //});
            //}
        //}
        ////////////////////////////////////////////////////////////////////////

        // Compute this universe
        self.calculate_dspin_dt();
        if !only_dspin_dt {
            self.calculate_acceleration_corrections();
        }

        ////////////////////////////////////////////////////////////////////////
        ////// Collect results from parallel threads
        //if self.consider_all_body_interactions {
            //self.parallel_universes.clear();
            //for _ in 0..n_parallel_universes {
                //// It will block the current thread if there no messages available
                //let parallel_universe = rx.recv().unwrap();
                //self.parallel_universes.push(parallel_universe);
            //}
        //}
        ////////////////////////////////////////////////////////////////////////


        if self.consider_all_body_interactions {
            // Integrate results from all parallel universes into this universe
            for parallel_universe in self.parallel_universes.iter() {
                //let mut num = 0;
                let id = self.n_particles - parallel_universe.n_particles - 1;

                for (particle, particle_parallel_universe) in self.particles[1+id..self.n_particles].iter_mut().zip(parallel_universe.particles[..parallel_universe.n_particles].iter()) {
                    //println!("Body {}:", num);
                    if !only_dspin_dt {
                        //println!("      Acceleration {:e} {:e} {:e}", particle.acceleration.x, particle.acceleration.y, particle.acceleration.z);
                        particle.acceleration.x += particle_parallel_universe.acceleration.x;
                        particle.acceleration.y += particle_parallel_universe.acceleration.y;
                        particle.acceleration.z += particle_parallel_universe.acceleration.z; 
                        //println!("      Acceleration {:e} {:e} {:e}", particle_parallel_universe.acceleration.x, particle_parallel_universe.acceleration.y, particle_parallel_universe.acceleration.z);
                        //println!("Added              {:e} {:e} {:e}", particle.acceleration.x, particle.acceleration.y, particle.acceleration.z);
                    }
                    //println!("      dspin_dt {:e} {:e} {:e}", particle.dspin_dt.x, particle.dspin_dt.y, particle.dspin_dt.z);
                    particle.dspin_dt.x += particle_parallel_universe.dspin_dt.x;
                    particle.dspin_dt.y += particle_parallel_universe.dspin_dt.y;
                    particle.dspin_dt.z += particle_parallel_universe.dspin_dt.z;
                    //println!("      dspin_dt {:e} {:e} {:e}", particle_parallel_universe.dspin_dt.x, particle_parallel_universe.dspin_dt.y, particle_parallel_universe.dspin_dt.z);
                    //println!("Added          {:e} {:e} {:e}", particle.dspin_dt.x, particle.dspin_dt.y, particle.dspin_dt.z);
                    //num += 1;
                }
                //println!("--------------------------------------------");
            }
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
        }
    }

    fn calculate_dspin_dt(&mut self) {
        let central_body = true;

        self.calculate_distance_and_velocities(); // Needed for calculate_torque_due_to_tides and calculate_planet_dependent_dissipation_factors
        self.calculate_planet_dependent_dissipation_factors(); // Needed by calculate_orthogonal_component_of_the_tidal_force and calculate_orthogonal_component_of_the_tidal_force if MathisSolarLike

        self.calculate_orthogonal_component_of_the_tidal_force(central_body);   // Needed for calculate_torque_due_to_tides and calculate_tidal_acceleration
        self.calculate_orthogonal_component_of_the_tidal_force(!central_body);  // Needed for calculate_torque_due_to_tides and calculate_tidal_acceleration
        self.calculate_torque_due_to_tides(central_body);   // Needed for spin integration
        self.calculate_torque_due_to_tides(!central_body);  // Needed for spin integration
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

            for particle in particles.iter_mut() {
                if !central_body {
                    reference_spin = particle.spin.clone();
                    orthogonal_component_of_the_tidal_force = particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
                } else {
                    orthogonal_component_of_the_tidal_force = particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
                }
                // Calculation of r scalar spin
                let rscalspin = particle.position.x * reference_spin.x 
                                + particle.position.y * reference_spin.y
                                + particle.position.z * reference_spin.z;

                // distance to star
                let distance = particle.distance;

                //// Torque calculation (star)
                // - Equation 8-9 from Bolmont et al. 2015
                let n_tid_x: f64 = orthogonal_component_of_the_tidal_force * (distance * reference_spin.x - rscalspin*particle.position.x/distance - 1.0/distance
                                * (particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y) );

                let n_tid_y :f64 = orthogonal_component_of_the_tidal_force * (distance * reference_spin.y - rscalspin*particle.position.y/distance - 1.0/distance
                                * (particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z) );

                let n_tid_z: f64 = orthogonal_component_of_the_tidal_force * (distance * reference_spin.z - rscalspin*particle.position.z/distance - 1.0/distance
                                * (particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x) );

                if central_body {
                    // Integration of the spin (total torque tides):
                    let factor = K2 / (star.mass_g + particle.mass_g);
                    //let factor = 1. / (star.mass + particle.mass);
                    torque.x += factor * n_tid_x;
                    torque.y += factor * n_tid_y;
                    torque.z += factor * n_tid_z;
                } else {
                    let factor = - K2 * star.mass_g / (particle.mass_g * (particle.mass_g + star.mass_g) 
                                    * particle.radius_of_gyration_2 * particle.radius.powi(2));
                    //let factor = - 1. * star.mass / (particle.mass * (particle.mass + star.mass) 
                                    //* particle.radius_of_gyration_2 * particle.radius.powi(2));
                    particle.torque.x = n_tid_x;
                    particle.torque.y = n_tid_y;
                    particle.torque.z = n_tid_z;
                    particle.dspin_dt.x = factor * n_tid_x;
                    particle.dspin_dt.y = factor * n_tid_y;
                    particle.dspin_dt.z = factor * n_tid_z;
                }

            }

            if central_body {
                let factor = - 1. / (star.radius_of_gyration_2 * star.radius.powi(2));
                star.torque.x = torque.x;
                star.torque.y = torque.y;
                star.torque.z = torque.z;
                star.dspin_dt.x = factor * torque.x;
                star.dspin_dt.y = factor * torque.y;
                star.dspin_dt.z = factor * torque.z;
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

                //println!("{:e} {:e}", particle.denergy_dt, factor2);
                //thread::sleep_ms(1000);

                //tmp = Ftidop/rr
                //tmp1 = 1.d0/rr*(Ftidr_diss + tmp*v_rad)

                //dEdt = -(tmp1*(xhx*vhx+xhy*vhy+xhz*vhz) &
                        //+ tmp*((spiny*xhz-spinz*xhy-vhx)*vhx &
                              //+(spinz*xhx-spinx*xhz-vhy)*vhy &
                              //+(spinx*xhy-spiny*xhx-vhz)*vhz)) &
                        //+ m(1)/(m(1)+m(i))*(N_tid_px*spinx+N_tid_py*spiny+N_tid_pz*spinz)
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
                let total_general_relativity_force_x = particle.mass_g 
                        * (radial_component_of_the_general_relativity_force * particle.position.x / particle.distance
                        + orthogonal_component_of_the_general_relativity_force * particle.velocity.x / particle.norm_velocity_vector);
                let total_general_relativity_force_y = particle.mass_g 
                        * (radial_component_of_the_general_relativity_force * particle.position.y / particle.distance
                        + orthogonal_component_of_the_general_relativity_force * particle.velocity.y / particle.norm_velocity_vector);
                let total_general_relativity_force_z = particle.mass_g 
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
            star.general_relativity_acceleration.x = -1.0 * factor2 * sum_total_general_relativity_force.x;
            star.general_relativity_acceleration.y = -1.0 * factor2 * sum_total_general_relativity_force.y;
            star.general_relativity_acceleration.z = -1.0 * factor2 * sum_total_general_relativity_force.z;
        }
    }


    fn calculate_acceleration_induced_by_rotational_flattering(&mut self) {
        if let Some((star, particles)) = self.particles[..self.n_particles].split_first_mut() {
            let factor2 = K2 / star.mass_g;
            //let factor2 = 1 / star.mass;

            // Calculation of the norm square of the spin for the star
            let normspin_2_star = star.spin.x.powi(2) + star.spin.y.powi(2) + star.spin.z.powi(2);
            let mut sum_total_force_induced_by_rotation = Axes{x:0., y:0., z:0.};

            for particle in particles.iter_mut() {
                let factor1 = K2 / particle.mass_g;
                
                // Calculation of the norm square of the spin for the planet
                let normspin_2_planet = particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2);

                // - Equation 16 from Bolmont et al. 2015
                let cpi = star.mass_g * particle.fluid_love_number * normspin_2_planet * particle.radius.powi(5) / (6. * K2); // Msun.AU^5.day-2
                let csi = particle.mass_g * star.fluid_love_number * normspin_2_star * star.radius.powi(5) / (6. * K2); // Msun.AU^5.day-2
                //println!("cpi, csi = {:e} {:e}", cpi, csi);
                
                // Calculation of r scalar spin for planet
                let rscalspin_planet = particle.position.x * particle.spin.x 
                                + particle.position.y * particle.spin.y
                                + particle.position.z * particle.spin.z;
                // Calculation of r scalar spin for the star
                let rscalspin_star = particle.position.x * star.spin.x 
                                + particle.position.y * star.spin.y
                                + particle.position.z * star.spin.z;
          
                // - Equation 15 from Bolmont et al. 2015
                let radial_component_of_the_force_induced_by_rotation = -3./particle.distance.powi(5) * (cpi + csi)
                    + 15./particle.distance.powi(7) * (csi * rscalspin_star * rscalspin_star/normspin_2_star
                        + cpi * rscalspin_planet * rscalspin_planet/normspin_2_planet); // Msun.AU.day-1
                    
                // - Equation 15 from Bolmont et al. 2015
                let orthogonal_component_of_the_force_induced_by_star_rotation = -6. * csi * rscalspin_star / (normspin_2_star * particle.distance.powi(5));

                // - Equation 15 from Bolmont et al. 2015
                let orthogonal_component_of_the_force_induced_by_planet_rotation = -6. * cpi * rscalspin_planet / (normspin_2_planet * particle.distance.powi(5));

                // - Equation 15 from Bolmont et al. 2015
                let total_force_induced_by_rotation_x = radial_component_of_the_force_induced_by_rotation * particle.position.x
                    + orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.x
                    + orthogonal_component_of_the_force_induced_by_star_rotation * star.spin.x;
                let total_force_induced_by_rotation_y = radial_component_of_the_force_induced_by_rotation * particle.position.y
                    + orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.y
                    + orthogonal_component_of_the_force_induced_by_star_rotation * star.spin.y;
                let total_force_induced_by_rotation_z = radial_component_of_the_force_induced_by_rotation * particle.position.z
                    + orthogonal_component_of_the_force_induced_by_planet_rotation * particle.spin.z
                    + orthogonal_component_of_the_force_induced_by_star_rotation * star.spin.z;
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
    
    pub fn evolve_particles(&mut self, current_time: f64) {
        for (particle, evolver) in self.particles[..self.n_particles].iter_mut().zip(self.particles_evolvers.iter_mut()) {
            particle.radius = evolver.radius(current_time, particle.radius);
            particle.radius_of_gyration_2 = evolver.radius_of_gyration_2(current_time, particle.radius_of_gyration_2);
            particle.love_number = evolver.love_number(current_time, particle.love_number);
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
