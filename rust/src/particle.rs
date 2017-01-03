use super::constants::{K2, G, N_PARTICLES, TIDES, GENERAL_RELATIVITY, SPEED_OF_LIGHT_2};
use super::integrator::IntegratorType;

#[derive(Debug,Copy, Clone)]
pub struct Axes {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug,Copy, Clone)]
pub struct Particle {
    pub mass: f64,
    pub mass_g: f64,
    pub radius: f64,
    pub dissipation_factor: f64, // sigma
    pub radius_of_gyration_2: f64, // radius of gyration square can be computed in terms of the mass moment of inertia, which 
                               // depends on the shape of the body and determines the torque needed for a desired angular acceleration
    pub love_number: f64,
    //
    pub position: Axes,
    pub velocity: Axes,
    pub acceleration: Axes,
    //
    pub radial_velocity: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub distance: f64,
    // Tides
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub radial_component_of_the_tidal_force_conservative_part: f64,
    pub radial_component_of_the_tidal_force_dissipative_part: f64,
    pub denergy_dt: f64,
    pub torque: Axes, // Force
    pub spin: Axes,
    pub dspin_dt: Axes,
    pub tidal_acceleration: Axes,
    pub general_relativity_acceleration: Axes,
    // General Relativity
    pub general_relativity_factor: f64,
}

impl Particle {
    pub fn new(mass: f64, radius: f64, dissipation_factor: f64, radius_of_gyration_2: f64, love_number: f64, position: Axes, velocity: Axes, acceleration: Axes, spin: Axes) -> Particle {
        let torque = Axes{x: 0., y: 0., z: 0.};
        let dspin_dt = Axes{x: 0., y: 0., z: 0.};
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        let general_relativity_acceleration = Axes{x: 0., y: 0., z: 0.};
        Particle { mass:mass, mass_g: mass*K2, radius:radius, dissipation_factor:dissipation_factor, radius_of_gyration_2:radius_of_gyration_2, love_number:love_number,
                    position:position, velocity:velocity, acceleration:acceleration, spin:spin,
                    radial_velocity: 0., norm_velocity_vector:0., norm_velocity_vector_2:0., distance:0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide:0., orthogonal_component_of_the_tidal_force_due_to_planetary_tide:0., radial_component_of_the_tidal_force:0., radial_component_of_the_tidal_force_conservative_part:0., radial_component_of_the_tidal_force_dissipative_part:0.,
                    denergy_dt:0., torque:torque, dspin_dt:dspin_dt, 
                    tidal_acceleration:tidal_acceleration, general_relativity_acceleration:general_relativity_acceleration,
                    general_relativity_factor: 0.}
    }
}


#[derive(Copy, Clone)]
pub struct Particles {
    pub particles : [Particle; N_PARTICLES],
    integrator_type: IntegratorType,
    integrator_is_whfasthelio: bool,
}

impl Particles {
    pub fn new(mut particles: [Particle; N_PARTICLES], integrator_type: IntegratorType) -> Particles {
        let integrator_is_whfasthelio = match integrator_type {
            IntegratorType::WHFastHelio => { true },
            _ => { false }
        };
        if GENERAL_RELATIVITY {
            let star_index = 0; // index
			let local_copy_particles = particles;
            let local_copy_star = &local_copy_particles[star_index];
            for particle in particles[1..].iter_mut() {
                particle.general_relativity_factor =  local_copy_star.mass_g*particle.mass_g / (local_copy_star.mass_g + particle.mass_g).powi(2)
            }
        }
        Particles {
                    particles: particles,
                    integrator_type: integrator_type,
                    integrator_is_whfasthelio: integrator_is_whfasthelio,
                    }
    }

    pub fn gravity_calculate_acceleration(&mut self) {
		let local_copy_particles = self.particles;

        for (i, particle_a) in self.particles.iter_mut().enumerate() {
			particle_a.acceleration.x = 0.;
			particle_a.acceleration.y = 0.;
			particle_a.acceleration.z = 0.;
            for (j, particle_b) in local_copy_particles.iter().enumerate() {
                if self.integrator_is_whfasthelio && (i == 0 || j == 0) {
                    // For WHFastHelio, ignore central body
                    continue
                }
                if i == j {
                    continue;
                }
                // gravitational force equation for vector quantities
                // Fx = - G * ((m1 * m2) / r^3) * (x2 - x1)
                // Fy = - G * ((m1 * m2) / r^3) * (y2 - y1)
                // Fz = - G * ((m1 * m2) / r^3) * (z2 - z1)
                //
                // acceleration:
                // ax = Fx / m1 = - G * (m2 / r^3) * (x2 - x1)
                let dx = particle_a.position.x - particle_b.position.x;
                let dy = particle_a.position.y - particle_b.position.y;
                let dz = particle_a.position.z - particle_b.position.z;
                let r = (dx*dx + dy*dy + dz*dz).sqrt();
                let prefact = -G/(r*r*r) * particle_b.mass;

                particle_a.acceleration.x += prefact * dx;
                particle_a.acceleration.y += prefact * dy;
                particle_a.acceleration.z += prefact * dz;
            }
        }

        if !self.integrator_is_whfasthelio {
            // Tides require heliocentric point of reference, the star should continue in the zero point
            // so we must compensate all the planets (but if WHFastHelio is being used, it is
            // automatically done by the integrator):
            let star_index = 0; // index
            let local_copy_star = &local_copy_particles[star_index];
            for particle in self.particles.iter_mut() {
                particle.acceleration.x -= local_copy_star.acceleration.x;
                particle.acceleration.y -= local_copy_star.acceleration.y;
                particle.acceleration.z -= local_copy_star.acceleration.z;
            }
            let star_index = 0;
            //let star = (&mut self.particles[star_index..star_index+1]).iter_mut().next().unwrap();
            let star = &mut self.particles[star_index];
            star.acceleration.x = 0.;
            star.acceleration.y = 0.;
            star.acceleration.z = 0.;
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

    pub fn calculate_additional_forces(&mut self, only_dspin_dt: bool) {
        let central_body = true;
        self.calculate_distance_and_velocities(); // Needed for calculate_torque_due_to_tides
        self.calculate_orthogonal_component_of_the_tidal_force(central_body);   // Needed for calculate_torque_due_to_tides and calculate_tidal_acceleration
        self.calculate_orthogonal_component_of_the_tidal_force(!central_body);  // Needed for calculate_torque_due_to_tides and calculate_tidal_acceleration
        self.calculate_torque_due_to_tides(central_body);   // Needed for spin integration
        self.calculate_torque_due_to_tides(!central_body);  // Needed for spin integration


        if !only_dspin_dt {
            if TIDES {
                self.calculate_radial_component_of_the_tidal_force();  // Needed for calculate_tidal_acceleration
                self.calculate_tidal_acceleration();
            }

            if GENERAL_RELATIVITY {
                self.calculate_general_relativity_acceleration()
            }

            // Add the tidal+flattening+general relativity accelerations to the gravitational one (already computed)
            for particle in self.particles[1..].iter_mut() {
                if TIDES {
                    particle.acceleration.x += particle.tidal_acceleration.x;
                    particle.acceleration.y += particle.tidal_acceleration.y;
                    particle.acceleration.z += particle.tidal_acceleration.z;
                }

                if GENERAL_RELATIVITY {
                    particle.acceleration.x += particle.general_relativity_acceleration.x;
                    particle.acceleration.y += particle.general_relativity_acceleration.y;
                    particle.acceleration.z += particle.general_relativity_acceleration.z;
                }
            }
        } 
        //else {
            //self.calculate_radial_component_of_the_tidal_force();  // Needed for calculate_tidal_acceleration
            //self.calculate_tidal_acceleration();
            //println!("atide!  {:e} {:e} {:e}", self.particles[1].tidal_acceleration.x, self.particles[1].tidal_acceleration.y, self.particles[1].tidal_acceleration.z);
        //}
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // TIDES
    fn calculate_torque_due_to_tides(&mut self, central_body:bool) {
		let local_copy_particles = self.particles;

        let star_index = 0; // index
        // Get only the particle that corresponds to the star
        // (using a slice of 1 element):
        //let local_copy_star = local_copy_particles[star_index..star_index+1].iter().next().unwrap();
        //let local_copy_star = local_copy_particles[star_index]; // Less optimal since it copies
        let local_copy_star = &local_copy_particles[star_index];
        let mut torque = Axes{x: 0., y: 0., z:0.};

        let mut spin_index = 0;
        if central_body {
            spin_index = star_index;
        }
        let mut local_copy_particle_as_ref_spin = &local_copy_particles[spin_index];
        let mut orthogonal_component_of_the_tidal_force: f64;
        for (i, particle) in self.particles[1..].iter_mut().enumerate() {
            if !central_body {
                spin_index = i + 1;
                local_copy_particle_as_ref_spin = &local_copy_particles[spin_index];
                orthogonal_component_of_the_tidal_force = particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
            } else {
                orthogonal_component_of_the_tidal_force = particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
            }
            // Calculation of r scalar spin
            let rscalspin = particle.position.x * local_copy_particle_as_ref_spin.spin.x 
                            + particle.position.y * local_copy_particle_as_ref_spin.spin.y
                            + particle.position.z * local_copy_particle_as_ref_spin.spin.z;

            // (distance to star)^7
            let distance = particle.distance;

            //// SBC
            //if i == 1 {
                //println!("Ftdiop {:e}", orthogonal_component_of_the_tidal_force);
            //}

            //// Torque calculation (star)
            // - Equation 8-9 from Bolmont et al. 2015
            let n_tid_x: f64 = orthogonal_component_of_the_tidal_force * (distance * local_copy_particle_as_ref_spin.spin.x - rscalspin*particle.position.x/distance - 1.0/distance
                            * (particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y) );

            let n_tid_y :f64 = orthogonal_component_of_the_tidal_force * (distance * local_copy_particle_as_ref_spin.spin.y - rscalspin*particle.position.y/distance - 1.0/distance
                            * (particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z) );

            let n_tid_z: f64 = orthogonal_component_of_the_tidal_force * (distance * local_copy_particle_as_ref_spin.spin.z - rscalspin*particle.position.z/distance - 1.0/distance
                            * (particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x) );
            // SBC
            //if !central_body {
                //println!("pos  {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z );
                //println!("vel  {:e} {:e} {:e}", particle.velocity.x, particle.velocity.y, particle.velocity.z );
                //println!("spin {:e} {:e} {:e}", local_copy_particle_as_ref_spin.spin.x, local_copy_particle_as_ref_spin.spin.y, local_copy_particle_as_ref_spin.spin.z );
                //println!("Ftidos {:e} ", orthogonal_component_of_the_tidal_force);
                //println!("Distance {:e} ", distance);
                //println!("Rscalspin {:e} ", rscalspin);
                //println!("N_tid_x {:e}", n_tid_x);
                //println!("N_tid_y {:e}", n_tid_y);
                //println!("N_tid_z {:e}", n_tid_z);
                //////println!("Star spin {:e} {:e} {:e}", local_copy_star.spin.x, local_copy_star.spin.y, local_copy_star.spin.z);
                //////println!("{} Planet spin {:e} {:e} {:e}", self.current_time/365.25, self.particles[1].spin.x, self.particles[1].spin.y, self.particles[1].spin.z);
            //}

            if central_body {
                // Integration of the spin (total torque tides):
                let factor = K2 / (local_copy_star.mass_g + particle.mass_g);
                //let factor = 1. / (local_copy_star.mass + particle.mass);
                torque.x += factor * n_tid_x;
                torque.y += factor * n_tid_y;
                torque.z += factor * n_tid_z;
            } else {
                let factor = - K2 * local_copy_star.mass_g / (particle.mass_g * (particle.mass_g + local_copy_star.mass_g) 
                                * particle.radius_of_gyration_2 * particle.radius.powi(2));
                //let factor = - 1. * local_copy_star.mass / (particle.mass * (particle.mass + local_copy_star.mass) 
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
            let factor = - 1. / (local_copy_star.radius_of_gyration_2 * local_copy_star.radius.powi(2));
            // Get only the particle that corresponds to the star in mutable form 
            // (using a slice of 1 element):
            //let star = (&mut self.particles[star_index..star_index+1]).iter_mut().next().unwrap();
            let star = &mut self.particles[star_index];
            star.torque.x = torque.x;
            star.torque.y = torque.y;
            star.torque.z = torque.z;
            star.dspin_dt.x = factor * torque.x;
            star.dspin_dt.y = factor * torque.y;
            star.dspin_dt.z = factor * torque.z;
        }

    }

    fn calculate_orthogonal_component_of_the_tidal_force(&mut self, central_body:bool) {
		let local_copy_particles = self.particles;
        let star_index = 0;
        let local_copy_star = &local_copy_particles[star_index];

        for particle in self.particles[1..].iter_mut() {
            // (distance to star)^7
            let distance_7 = particle.distance.powi(7);

            //// Tidal force calculation (star) :: Only orthogonal component is needed
            // - Equation 5 from Bolmont et al. 2015
            // - Ftides in Msun.AU.day-1
            if central_body {
                // - F_tides_ortho_star
                particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass_g.powi(2))
                                                * (local_copy_star.radius.powi(10)) 
                                                * local_copy_star.dissipation_factor / ( (K2.powi(2)) * distance_7);
                //particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass.powi(2))
                                                    //* (local_copy_star.radius.powi(10)) 
                                                    //* local_copy_star.dissipation_factor / (distance_7);
            } else {
                // - F_tides_ortho_plan
                particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (local_copy_star.mass_g.powi(2))
                                                * (particle.radius.powi(10))
                                                * particle.dissipation_factor    / ( (K2.powi(2)) * distance_7);
                //particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (local_copy_star.mass.powi(2))
                                                //* (particle.radius.powi(10))
                                                //* particle.dissipation_factor    / (distance_7);

                // SBC
                //println!("> {:e} {:e} {:e} {:e}", local_copy_star.mass_g, particle.radius.powi(10), particle.dissipation_factor, distance_7);
                //println!("> {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z);
            }
        }
    }

    fn calculate_radial_component_of_the_tidal_force(&mut self) {
        // dEdt_tides
        // F_tides_rad
        // F_tides_rad_cons
        // F_tides_rad_diss
		let local_copy_particles = self.particles;
        let star_index = 0;
        let local_copy_star = &local_copy_particles[star_index];
        let star_mass_2 = local_copy_star.mass_g * local_copy_star.mass_g;

        for particle in self.particles[1..].iter_mut() {
            let planet_mass_2 = particle.mass_g * particle.mass_g;
            // Conservative part of the radial tidal force
            // - Ftidr_cons

            particle.radial_component_of_the_tidal_force_conservative_part = -3.0 / (particle.distance.powi(7) * K2)
                        * (planet_mass_2 * local_copy_star.radius.powi(5) * local_copy_star.love_number 
                        + star_mass_2 * particle.radius.powi(5) * particle.love_number);
            //particle.radial_component_of_the_tidal_force_conservative_part = -3.0 / particle.distance.powi(7)
                        //* (planet_mass_2 * local_copy_star.radius.powi(5) * local_copy_star.love_number 
                        //+ star_mass_2 * particle.radius.powi(5) * particle.love_number);

            // Dissipative part of the radial tidal force
            // - Ftidr_diss
            let factor1 = -13.5 * particle.radial_velocity / (particle.distance.powi(8) * K2*K2);
            //let factor1 = -13.5 * particle.radial_velocity / particle.distance.powi(8);
            let term1 = planet_mass_2
                        * local_copy_star.radius.powi(10)
                        * local_copy_star.dissipation_factor;
            let term2 = star_mass_2
                        * particle.radius.powi(10)
                        * particle.dissipation_factor;
            particle.radial_component_of_the_tidal_force_dissipative_part = factor1 * (term1 + term2 );

            // If we consider the star as a point mass (used for denergy_dt calculation):
            let radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass = factor1 * term2;

            // Sum of the dissipative and conservative part of the radial force
            particle.radial_component_of_the_tidal_force = particle.radial_component_of_the_tidal_force_conservative_part + particle.radial_component_of_the_tidal_force_dissipative_part;

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
                        + local_copy_star.mass_g/(local_copy_star.mass_g + particle.mass_g) 
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


    fn calculate_tidal_acceleration(&mut self) {
		let local_copy_particles = self.particles;
        let star_index = 0;
        let local_copy_star = &local_copy_particles[star_index];
        let factor2 = K2 / local_copy_star.mass_g;
        //let factor2 = 1 / local_copy_star.mass;

        let mut sum_total_tidal_force = Axes{x:0., y:0., z:0.};
        for particle in self.particles[1..].iter_mut() {
            let factor1 = K2 / particle.mass_g;
            //let factor1 = 1. / particle.mass;

            // - Equation 6 from Bolmont et al. 2015
            let factor3 = particle.radial_component_of_the_tidal_force
                            + (particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide) * particle.radial_velocity / particle.distance;
            let total_tidal_force_x = factor3 * particle.position.x / particle.distance
                                    + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                        * (local_copy_star.spin.y * particle.position.z  - local_copy_star.spin.z * particle.position.y - particle.velocity.x)
                                    + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                        * (particle.spin.y * particle.position.z  - particle.spin.z * particle.position.y - particle.velocity.x);
            let total_tidal_force_y = factor3 * particle.position.y / particle.distance
                                    + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                        * (local_copy_star.spin.z * particle.position.x  - local_copy_star.spin.x * particle.position.z - particle.velocity.y)
                                    + particle.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.distance 
                                        * (particle.spin.z * particle.position.x  - particle.spin.x * particle.position.z - particle.velocity.y);
            let total_tidal_force_z = factor3 * particle.position.z / particle.distance
                                    + particle.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.distance 
                                        * (local_copy_star.spin.x * particle.position.y  - local_copy_star.spin.y * particle.position.x - particle.velocity.z)
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
        for particle in self.particles[1..].iter_mut() {
            particle.tidal_acceleration.x += factor2 * sum_total_tidal_force.x;
            particle.tidal_acceleration.y += factor2 * sum_total_tidal_force.y;
            particle.tidal_acceleration.z += factor2 * sum_total_tidal_force.z;
        }
    }

    fn calculate_general_relativity_acceleration(&mut self) {
		let local_copy_particles = self.particles;
        let star_index = 0;
        let local_copy_star = &local_copy_particles[star_index];
        let factor2 = K2 / local_copy_star.mass_g;
        //let factor2 = 1 / local_copy_star.mass;

        let mut sum_total_general_relativity_force = Axes{x:0., y:0., z:0.};
        for particle in self.particles[1..].iter_mut() {
            let factor1 = K2 / particle.mass_g;
            //let factor1 = 1. / particle.mass;

            // Radial part of the GR force (Kidder 1995, Mardling & Lin 2002)
            // - Equation 11 from Bolmont et al. 2015
            let star_planet_mass_g = local_copy_star.mass_g + particle.mass_g;
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
                    * (radial_component_of_the_general_relativity_force * particle.position.x / particle.norm_velocity_vector 
                    + orthogonal_component_of_the_general_relativity_force * particle.velocity.x / particle.distance);
            let total_general_relativity_force_y = particle.mass_g 
                    * (radial_component_of_the_general_relativity_force * particle.position.y / particle.norm_velocity_vector
                    + orthogonal_component_of_the_general_relativity_force * particle.velocity.y / particle.distance);
            let total_general_relativity_force_z = particle.mass_g 
                    * (radial_component_of_the_general_relativity_force * particle.position.z / particle.norm_velocity_vector
                    + orthogonal_component_of_the_general_relativity_force * particle.velocity.z / particle.distance);
            
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
        for particle in self.particles[1..].iter_mut() {
            particle.general_relativity_acceleration.x += factor2 * sum_total_general_relativity_force.x;
            particle.general_relativity_acceleration.y += factor2 * sum_total_general_relativity_force.y;
            particle.general_relativity_acceleration.z += factor2 * sum_total_general_relativity_force.z;
        }
    }

    fn calculate_distance_and_velocities(&mut self) {
        // Calculation of velocity vv(j), radial velocity vrad(j)
        // velocities in AU/day
        for particle in self.particles[1..].iter_mut() {
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

    pub fn compute_total_energy(&self) -> f64 {
        let mut e_kin = 0.;
        let mut e_pot = 0.;
        let e_offset = 0.; // Energy offset due to collisions and ejections

        // Kinectic energy
        for particle in self.particles.iter() {
            e_kin += 0.5 * particle.mass * (particle.velocity.x.powi(2) + particle.velocity.y.powi(2) + particle.velocity.z.powi(2));
        }
        // Gravitationl potential energy
        for (i, particle_a) in self.particles.iter().enumerate() {
            for particle_b in self.particles[i+1..].iter() {
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
        for particle in self.particles.iter() {
            total_angular_momentum.x += particle.mass*(particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y);
            total_angular_momentum.y += particle.mass*(particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z);
            total_angular_momentum.z += particle.mass*(particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x);
        }
        let total_angular_momentum = (total_angular_momentum.x.powf(2.) + total_angular_momentum.y.powf(2.) + total_angular_momentum.z.powf(2.)).sqrt();
        total_angular_momentum
    }

}
