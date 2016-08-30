use super::constants::{K2, G, N_PARTICLES};

#[derive(Copy, Clone)]
pub struct Axes {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Copy, Clone)]
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
    pub distance: f64,
    // Tides
    pub orthogonal_component_of_the_tidal_force: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub radial_component_of_the_tidal_force_conservative_part: f64,
    pub radial_component_of_the_tidal_force_dissipative_part: f64,
    pub denergy_dt: f64,
    pub torque: Axes, // Force
    pub spin: Axes,
    pub dspin_dt: Axes,
    pub tidal_acceleration: Axes,
}

impl Particle {
    pub fn new(mass: f64, radius: f64, dissipation_factor: f64, radius_of_gyration_2: f64, love_number: f64, position: Axes, velocity: Axes, acceleration: Axes, spin: Axes) -> Particle {
        let torque = Axes{x: 0., y: 0., z: 0.};
        let dspin_dt = Axes{x: 0., y: 0., z: 0.};
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        Particle { mass:mass, mass_g: mass*K2, radius:radius, dissipation_factor:dissipation_factor, radius_of_gyration_2:radius_of_gyration_2, love_number:love_number,
                    position:position, velocity:velocity, acceleration:acceleration, spin:spin,
                    radial_velocity: 0., norm_velocity_vector:0., distance:0.,
                    orthogonal_component_of_the_tidal_force:0., radial_component_of_the_tidal_force:0., radial_component_of_the_tidal_force_conservative_part:0., radial_component_of_the_tidal_force_dissipative_part:0.,
                    denergy_dt:0., torque:torque, dspin_dt:dspin_dt, tidal_acceleration:tidal_acceleration }
    }
}


#[derive(Copy, Clone)]
pub struct Particles {
    pub particles : [Particle; N_PARTICLES],
}

impl Particles {
    pub fn new(particles: [Particle; N_PARTICLES]) -> Particles {
        Particles {
                    particles:particles,
                    }
    }

    pub fn gravity_calculate_acceleration(&mut self) {
		let local_copy_particles = self.particles;

        for (i, particle_a) in self.particles.iter_mut().enumerate() {
			particle_a.acceleration.x = 0.;
			particle_a.acceleration.y = 0.;
			particle_a.acceleration.z = 0.;
            for (j, particle_b) in local_copy_particles.iter().enumerate() {
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

        // Tides require heliocentric point of reference, the star should continue in the zero point
        // so we must compensate all the planets:
        for particle in self.particles.iter_mut() {
            particle.acceleration.x -= particle.acceleration.x;
            particle.acceleration.y -= particle.acceleration.y;
            particle.acceleration.z -= particle.acceleration.z;
        }
        let star_index = 0;
        let star = (&mut self.particles[star_index..star_index+1]).iter_mut().next().unwrap();
        star.acceleration.x = 0.;
        star.acceleration.y = 0.;
        star.acceleration.z = 0.;
    }

    pub fn calculate_additional_forces(&mut self) {
        let central_body = true;
        self.calculate_distance_and_velocities(); // Needed for calculate_torque_due_to_tides
        self.calculate_orthogonal_component_of_the_tidal_force(central_body);   // Needed for calculate_torque_due_to_tides and calculate_tidal_acceleration
        self.calculate_orthogonal_component_of_the_tidal_force(!central_body);  // Needed for calculate_torque_due_to_tides and calculate_tidal_acceleration
        self.calculate_torque_due_to_tides(central_body);   // Needed for spin integration
        self.calculate_torque_due_to_tides(!central_body);  // Needed for spin integration
        self.calculate_radial_component_of_the_tidal_force();  // Needed for calculate_tidal_acceleration
        self.calculate_tidal_acceleration();
        // Add the tidal acceleration to the gravitational one (already computed)
        for particle in self.particles[1..].iter_mut() {
            particle.acceleration.x += particle.tidal_acceleration.x;
            particle.acceleration.y += particle.tidal_acceleration.y;
            particle.acceleration.z += particle.tidal_acceleration.z;
        }
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
        for (i, particle) in self.particles[1..].iter_mut().enumerate() {
            if !central_body {
                spin_index = i;
                local_copy_particle_as_ref_spin = &local_copy_particles[spin_index];
            }
            // Calculation of r scalar spin
            let rscalspin = particle.position.x * local_copy_particle_as_ref_spin.spin.x 
                            + particle.position.y * local_copy_particle_as_ref_spin.spin.y
                            + particle.position.z * local_copy_particle_as_ref_spin.spin.z;

            // (distance to star)^7
            let distance = particle.distance;

            //// SBC
            //if i == 1 {
                //println!("Ftdiop {:e}", particle.orthogonal_component_of_the_tidal_force);
            //}

            //// Torque calculation (star)
            // - Equation 10 in Bolmont et al. 2015
            let n_tid_x: f64 = particle.orthogonal_component_of_the_tidal_force * (distance * local_copy_particle_as_ref_spin.spin.x - rscalspin*particle.position.x/distance - 1.0/distance
                            * (particle.position.y*particle.velocity.z - particle.position.z*particle.velocity.y) );

            let n_tid_y :f64 = particle.orthogonal_component_of_the_tidal_force * (distance * local_copy_particle_as_ref_spin.spin.y - rscalspin*particle.position.y/distance - 1.0/distance
                            * (particle.position.z*particle.velocity.x - particle.position.x*particle.velocity.z) );

            let n_tid_z: f64 = particle.orthogonal_component_of_the_tidal_force * (distance * local_copy_particle_as_ref_spin.spin.z - rscalspin*particle.position.z/distance - 1.0/distance
                            * (particle.position.x*particle.velocity.y - particle.position.y*particle.velocity.x) );
            //// SBC
            //if i == 1 && !central_body && self.current_time > 0. {
                ////println!("Ftidos {:e} ", particle.orthogonal_component_of_the_tidal_force);
                //////println!("Distance {:e} ", distance);
                //////println!("Rscalspin {:e} ", rscalspin);
                //////println!("N_tid_z {:e} {:e} {:e} {:e}", n_tid_z, particle.orthogonal_component_of_the_tidal_force, distance, rscalspin);
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
                                * particle.radius_of_gyration_2 * particle.radius.powf(2.));
                //let factor = - 1. * local_copy_star.mass / (particle.mass * (particle.mass + local_copy_star.mass) 
                                //* particle.radius_of_gyration_2 * particle.radius.powf(2.));
                particle.torque.x = n_tid_x;
                particle.torque.y = n_tid_y;
                particle.torque.z = n_tid_z;
                particle.dspin_dt.x = factor * n_tid_x;
                particle.dspin_dt.y = factor * n_tid_y;
                particle.dspin_dt.z = factor * n_tid_z;
            }

        }

        if central_body {
            let factor = - 1. / (local_copy_star.radius_of_gyration_2 * local_copy_star.radius.powf(2.));
            // Get only the particle that corresponds to the star in mutable form 
            // (using a slice of 1 element):
            let star = (&mut self.particles[star_index..star_index+1]).iter_mut().next().unwrap();
            star.torque.x = torque.x;
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
            let distance_7 = particle.distance.powf(7.);

            //// Tidal force calculation (star) :: Only orthogonal component is needed
            // - Equation 6 from Bolmont et al. 2015
            // - Ftides in Msun.AU.day-1
            if central_body {
                // - F_tides_ortho_star
                particle.orthogonal_component_of_the_tidal_force = 4.5 * (particle.mass_g.powf(2.))
                                                * (local_copy_star.radius.powf(10.)) 
                                                * local_copy_star.dissipation_factor / ( (K2.powf(2.)) * distance_7);
                //particle.orthogonal_component_of_the_tidal_force = 4.5 * (particle.mass.powf(2))
                                                    //* (local_copy_star.radius.powf(10.)) 
                                                    //* local_copy_star.dissipation_factor / (distance_7);
            } else {
                // - F_tides_ortho_plan
                particle.orthogonal_component_of_the_tidal_force = 4.5 * (local_copy_star.mass_g.powf(2.))
                                                * (particle.radius.powf(10.))
                                                * particle.dissipation_factor    / ( (K2.powf(2.)) * distance_7);
                //particle.orthogonal_component_of_the_tidal_force = 4.5 * (local_copy_star.mass.powf(2.))
                                                //* (particle.radius.powf(10.))
                                                //* particle.dissipation_factor    / (distance_7);

                // SBC
                //println!("\t {:e} {:e} {:e} {:e}", local_copy_star.mass_g, particle.radius.powf(10.), particle.dissipation_factor, distance_7);
                //println!("**AKI {:e} {:e} {:e}", particle.position.x, particle.position.y, particle.position.z);
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

            particle.radial_component_of_the_tidal_force_conservative_part = -3.0 / (particle.distance.powf(8.) * K2*K2)
                        * (planet_mass_2 * local_copy_star.radius.powf(5.) * local_copy_star.love_number 
                        + star_mass_2 * particle.radius.powf(5.) * particle.love_number);
            //particle.radial_component_of_the_tidal_force_conservative_part = -3.0 / particle.distance.powf(8.)
                        //* (planet_mass_2 * local_copy_star.radius.powf(5.) * local_copy_star.love_number 
                        //+ star_mass_2 * particle.radius.powf(5.) * particle.love_number);

            // Dissipative part of the radial tidal force
            // - Ftidr_diss
            let factor1 = -13.5 * particle.radial_velocity / (particle.distance.powf(8.) * K2*K2);
            //let factor1 = -13.5 * particle.radial_velocity / particle.distance.powf(8.);
            let term1 = planet_mass_2
                        * local_copy_star.radius.powf(10.)
                        * local_copy_star.dissipation_factor;
            let term2 = star_mass_2
                        * particle.radius.powf(10.)
                        * particle.dissipation_factor;
            particle.radial_component_of_the_tidal_force_dissipative_part = factor1 * (term1 + term2 );

            // If we consider the star as a point mass (used for denergy_dt calculation):
            let radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass = factor1 * term2;

            // Sum of the dissipative and conservative part of the radial force
            particle.radial_component_of_the_tidal_force = particle.radial_component_of_the_tidal_force_conservative_part + particle.radial_component_of_the_tidal_force_dissipative_part;

            ////// Instantaneous energy loss dE/dt due to tides
            ////// in Msun.AU^2.day^(-3)
            //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
            let factor2 = particle.orthogonal_component_of_the_tidal_force / particle.distance;
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
        //// acc_tides
		let local_copy_particles = self.particles;
        let star_index = 0;
        let local_copy_star = &local_copy_particles[star_index];

        for particle in self.particles[1..].iter_mut() {
            let factor1 = K2 / particle.mass_g;
            //let factor1 = 1. / particle.mass;
            let factor2 = particle.radial_component_of_the_tidal_force
                            + (local_copy_star.orthogonal_component_of_the_tidal_force + particle.orthogonal_component_of_the_tidal_force) * particle.radial_velocity / particle.distance;
            particle.tidal_acceleration.x = factor1 * (factor2 * particle.position.x / particle.distance
                                    + local_copy_star.orthogonal_component_of_the_tidal_force / particle.distance 
                                        * (local_copy_star.spin.y * local_copy_star.position.z  - local_copy_star.spin.z * local_copy_star.position.y - local_copy_star.velocity.x)
                                    + particle.orthogonal_component_of_the_tidal_force / particle.distance 
                                        * (particle.spin.y * local_copy_star.position.z  - particle.spin.z * local_copy_star.position.y - local_copy_star.velocity.x)
                                    );
            particle.tidal_acceleration.y = factor1 * (factor2 * particle.position.y / particle.distance
                                    + local_copy_star.orthogonal_component_of_the_tidal_force / particle.distance 
                                        * (local_copy_star.spin.x * local_copy_star.position.x  - local_copy_star.spin.x * local_copy_star.position.z - local_copy_star.velocity.y)
                                    + particle.orthogonal_component_of_the_tidal_force / particle.distance 
                                        * (particle.spin.x * local_copy_star.position.x  - particle.spin.x * local_copy_star.position.z - local_copy_star.velocity.y)
                                    );
            particle.tidal_acceleration.z = factor1 * (factor2 * particle.position.z / particle.distance
                                    + local_copy_star.orthogonal_component_of_the_tidal_force / particle.distance 
                                        * (local_copy_star.spin.x * local_copy_star.position.y  - local_copy_star.spin.y * local_copy_star.position.x - local_copy_star.velocity.z)
                                    + particle.orthogonal_component_of_the_tidal_force / particle.distance 
                                        * (particle.spin.x * local_copy_star.position.y  - particle.spin.y * local_copy_star.position.x - local_copy_star.velocity.z)
                                    );
        }

     }


    fn calculate_distance_and_velocities(&mut self) {
        // Calculation of velocity vv(j), radial velocity vrad(j)
        // velocities in AU/day
        for particle in self.particles[1..].iter_mut() {
            // Norm of the velocity
            let v_2 = (particle.velocity.x.powf(2.)) 
                                + (particle.velocity.y.powf(2.))
                                + (particle.velocity.z.powf(2.));
            let norm_vel = v_2.sqrt();
            particle.norm_velocity_vector = norm_vel;

            // (distance to star)^2
            let distance_2 = (particle.position.x.powf(2.)) 
                                + (particle.position.y.powf(2.))
                                + (particle.position.z.powf(2.));
            let distance = distance_2.sqrt();
            particle.distance = distance;

            // Radial velocity
            let v_rad = (particle.position.x*particle.velocity.x +
                        particle.position.y*particle.velocity.y +
                        particle.position.z*particle.velocity.z)
                        / distance;
            particle.radial_velocity = v_rad;
        }
    }

}
