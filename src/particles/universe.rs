extern crate time;
use std::collections::HashMap;
use super::super::constants::{G, MAX_PARTICLES, MAX_DISTANCE_2};
use super::{Evolver, EvolutionType};
use super::{Particle, Disk};
use super::{Axes};
use super::{tides, flattening, general_relativity, evolver, wind, disk, common};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Universe {
    pub initial_time: f64,
    pub time_limit: f64,
    pub particles: [Particle; MAX_PARTICLES],
    pub particles_evolvers: Vec<Evolver>,
    pub n_particles: usize,
    pub evolving_particles_exist: bool,
    pub wind_effects_exist: bool,
    pub consider_tides: bool,
    pub consider_disk_interaction: bool,
    pub consider_rotational_flattening: bool,
    pub consider_general_relativity: ConsiderGeneralRelativity,
    host_particle_index: usize, // Most massive particle
    tidal_host_particle_index: usize, // Particle that is the main one for tidal effects
    disk_host_particle_index: usize, // Particle with disk
    star_planet_dependent_dissipation_factors : HashMap<usize, f64>, // Central body specific
    temporary_copied_particle_positions: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particle_velocities: [Axes; MAX_PARTICLES], // For optimization purposes (TODO: Delete and adapt python package)
    temporary_copied_particles_masses: [f64; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_radiuses: [f64; MAX_PARTICLES], // For optimization purposes
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum IgnoreGravityTerms {
    None,
    WHFastOne,
    WHFastTwo,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum ConsiderGeneralRelativity {
    None,
    Kidder1995, // MercuryT
    Anderson1975, // REBOUNDx gr
    Newhall1983, // REBOUNDx gr full
}

impl Universe {
    pub fn new(initial_time: f64, time_limit: f64, mut particles: Vec<Particle>, consider_tides: bool, consider_rotational_flattening: bool, consider_disk_interaction: bool, consider_general_relativity: ConsiderGeneralRelativity) -> Universe {

        let mut evolving_particles_exist = false;
        let mut wind_effects_exist = false;
        for particle in particles.iter(){
            if particle.wind_k_factor != 0. {
                wind_effects_exist = true;
            }
            match particle.evolution_type {
                EvolutionType::NonEvolving => {},
                _ => {
                    evolving_particles_exist = true;
                },
            }
            if wind_effects_exist && evolving_particles_exist {
                break;
            }
        }

        // Most massive particle
        let host_particle_index = 0;
        
        // Particle that is the main one for tidal effects
        let tidal_host_particle_index = 0;


        // Find the position of the particle that has the disk if any
        let mut disk_host_particle_index = MAX_PARTICLES+1;
        if consider_disk_interaction {
            for (i, particle) in particles.iter().enumerate() {
                if let Disk::Host( disk ) = particle.disk {
                    if disk.inner_edge_distance != 0. || disk.outer_edge_distance != 0. {
                        if disk_host_particle_index == MAX_PARTICLES+1 {
                            disk_host_particle_index = i;
                        } else {
                            panic!("Only one body with a disk is allowed!");
                        }
                    }
                }
            }
        }

        if consider_general_relativity != ConsiderGeneralRelativity::None {
            let star_index = 0; // index
            let local_copy_star_mass_g = particles[star_index].mass_g;
            for particle in particles[1..].iter_mut() {
                particle.general_relativity_factor =  local_copy_star_mass_g*particle.mass_g / (local_copy_star_mass_g + particle.mass_g).powi(2)
            }
        }

        // Populate inertial positions/velocities from heliocentric positions
        let (center_of_mass_position, center_of_mass_velocity) = calculate_center_of_mass(&particles);
        for particle in particles.iter_mut() {
            particle.inertial_position.x = particle.position.x - center_of_mass_position.x;
            particle.inertial_position.y = particle.position.y - center_of_mass_position.y;
            particle.inertial_position.z = particle.position.z - center_of_mass_position.z;
            particle.inertial_velocity.x = particle.velocity.x - center_of_mass_velocity.x;
            particle.inertial_velocity.y = particle.velocity.y - center_of_mass_velocity.y;
            particle.inertial_velocity.z = particle.velocity.z - center_of_mass_velocity.z;
        }


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
        for _i in n_particles..MAX_PARTICLES {
            // For dummy particles
            particles_evolvers.push(Evolver::new(EvolutionType::NonEvolving, initial_time, time_limit));
        }


        let temporary_copied_particle_positions = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES];
        let temporary_copied_particle_velocities = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES];
        let temporary_copied_particles_masses = [0.; MAX_PARTICLES];
        let temporary_copied_particles_radiuses = [0.; MAX_PARTICLES];

        let universe = Universe {
                    initial_time: initial_time,
                    time_limit: time_limit,
                    particles: transformed_particles,
                    particles_evolvers: particles_evolvers,
                    n_particles: n_particles,
                    evolving_particles_exist: evolving_particles_exist,
                    wind_effects_exist: wind_effects_exist,
                    consider_tides: consider_tides,
                    consider_rotational_flattening: consider_rotational_flattening,
                    consider_general_relativity: consider_general_relativity,
                    consider_disk_interaction: consider_disk_interaction,
                    host_particle_index: host_particle_index,
                    tidal_host_particle_index: tidal_host_particle_index,
                    disk_host_particle_index: disk_host_particle_index,
                    star_planet_dependent_dissipation_factors:HashMap::new(),
                    temporary_copied_particle_positions: temporary_copied_particle_positions,
                    temporary_copied_particle_velocities: temporary_copied_particle_velocities,
                    temporary_copied_particles_masses: temporary_copied_particles_masses,
                    temporary_copied_particles_radiuses: temporary_copied_particles_radiuses,
                    };
        universe
    }

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
                let prefact = -G/(distance*distance*distance) * particle_b_mass;

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

    pub fn calculate_particles_evolving_quantities(&mut self, current_time: f64) {
        if self.evolving_particles_exist {
            let (mut particles, _) = self.particles.split_at_mut(self.n_particles);
            evolver::calculate_particles_evolving_quantities(current_time, &mut particles, &mut self.particles_evolvers);
        }
    }

    pub fn calculate_norm_spin(&mut self) {
        if self.evolving_particles_exist || self.consider_rotational_flattening {
            let (mut particles, _) = self.particles.split_at_mut(self.n_particles);
            common::calculate_norm_spin(&mut particles); // Needed for rotational flattening (torque and accelerations) and evolution
        }
    }

    pub fn calculate_additional_effects(&mut self, current_time: f64, evolution: bool, dangular_momentum_dt_per_moment_of_inertia: bool, accelerations: bool, ignored_gravity_terms: IgnoreGravityTerms) {
        let (mut particles, _) = self.particles.split_at_mut(self.n_particles);
        if (evolution && self.evolving_particles_exist) || 
            ((dangular_momentum_dt_per_moment_of_inertia || accelerations) && self.consider_rotational_flattening) {
            common::calculate_norm_spin(&mut particles); // Needed for rotational flattening (torque and accelerations) and evolution
        }
        
        {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.tidal_host_particle_index);
            if let Some((tidal_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                if (dangular_momentum_dt_per_moment_of_inertia && (self.consider_tides || self.consider_rotational_flattening)) ||
                    (accelerations && (self.consider_tides  || self.consider_disk_interaction || self.consider_rotational_flattening || 
                                       self.consider_general_relativity != ConsiderGeneralRelativity::None)) {
                        common::calculate_distance_and_velocities(tidal_host_particle, &mut particles_left, &mut particles_right); // Needed for tides, rotational flattening, general relativity and evolution
                }
            }
        }
       
        if evolution && self.evolving_particles_exist {
            evolver::calculate_particles_evolving_quantities(current_time, &mut particles, &mut self.particles_evolvers);
        }

        if dangular_momentum_dt_per_moment_of_inertia && self.wind_effects_exist {
            wind::calculate_wind_factor(&mut particles);
        }

        {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.tidal_host_particle_index);
            if let Some((mut tidal_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                if (dangular_momentum_dt_per_moment_of_inertia && (self.consider_tides || self.consider_rotational_flattening)) ||
                    (accelerations && (self.consider_tides || self.consider_disk_interaction || self.consider_rotational_flattening || 
                                       self.consider_general_relativity != ConsiderGeneralRelativity::None)) {

                    {
                        //// calculate_orthogonal_components

                        tides::calculate_planet_dependent_dissipation_factors(&mut tidal_host_particle, &mut particles_left, &mut particles_right, &mut self.star_planet_dependent_dissipation_factors); // Needed by calculate_orthogonal_component_of_the_tidal_force and calculate_orthogonal_component_of_the_tidal_force if BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true)
                        tides::calculate_scalar_product_of_vector_position_with_spin(&mut tidal_host_particle, &mut particles_left, &mut particles_right); // Needed by tides and rotational flattening
                        
                        if self.consider_tides {
                            tides::calculate_orthogonal_component_of_the_tidal_force(&mut tidal_host_particle, &mut particles_left, &mut particles_right, &mut self.star_planet_dependent_dissipation_factors);
                        }

                        if self.consider_rotational_flattening {
                            flattening::calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(&mut tidal_host_particle, &mut particles_left, &mut particles_right);
                        }
                    }

                    if dangular_momentum_dt_per_moment_of_inertia && (self.consider_tides || self.consider_rotational_flattening) {
                        // Not needed for additional accelerations
                        {
                            //// calculate_torques // Needed for dangular_momentum_dt_per_moment_of_inertia
                            let central_body = true;

                            if self.consider_tides {
                                tides::calculate_torque_due_to_tides(&mut tidal_host_particle, &mut particles_left, &mut particles_right, central_body);
                                tides::calculate_torque_due_to_tides(&mut tidal_host_particle, &mut particles_left, &mut particles_right, !central_body);
                            }

                            if self.consider_rotational_flattening {
                                flattening::calculate_torque_induced_by_rotational_flattening(&mut tidal_host_particle, &mut particles_left, &mut particles_right, central_body);
                                flattening::calculate_torque_induced_by_rotational_flattening(&mut tidal_host_particle, &mut particles_left, &mut particles_right, !central_body);
                            }
                        }
                    }

                    if accelerations && (self.consider_tides || self.consider_rotational_flattening) {
                        if self.consider_tides {
                            tides::calculate_radial_component_of_the_tidal_force(&mut tidal_host_particle, &mut particles_left, &mut particles_right, &mut self.star_planet_dependent_dissipation_factors);  // Needed for calculate_tidal_acceleration
                            tides::calculate_tidal_acceleration(&mut tidal_host_particle, &mut particles_left, &mut particles_right);
                        }

                        if self.consider_rotational_flattening {
                            flattening::calculate_radial_component_of_the_force_induced_by_rotational_flattening(&mut tidal_host_particle, &mut particles_left, &mut particles_right);
                            flattening::calculate_acceleration_induced_by_rotational_flattering(&mut tidal_host_particle, &mut particles_left, &mut particles_right);
                        }


                    }
                }
            }
        }
        if accelerations && self.consider_general_relativity != ConsiderGeneralRelativity::None {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.host_particle_index);
            if let Some((mut host_particle, mut particles_right)) = particles_right.split_first_mut() {
                match self.consider_general_relativity {
                    ConsiderGeneralRelativity::Kidder1995 => {
                        general_relativity::calculate_kidder1995_general_relativity_acceleration(&mut host_particle, &mut particles_left, &mut particles_right);
                    },
                    ConsiderGeneralRelativity::Anderson1975 => {
                        general_relativity::calculate_anderson1975_general_relativity_acceleration(&mut host_particle, &mut particles_left, &mut particles_right, 
                                                                                                   ignored_gravity_terms);
                    },
                    ConsiderGeneralRelativity::Newhall1983 => {
                        general_relativity::calculate_newhall1983_general_relativity_acceleration(&mut host_particle, &mut particles_left, &mut particles_right, 
                                                                                                  ignored_gravity_terms);
                    },
                    ConsiderGeneralRelativity::None => {}
                }
            }
        }
        if accelerations && self.consider_disk_interaction {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.disk_host_particle_index);
            if let Some((mut disk_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                if let Disk::Host(_) = disk_host_particle.disk {
                    disk::calculate_disk_interaction_acceleration(current_time, &mut disk_host_particle, &mut particles_left, &mut particles_right);
                }
            }
        }
        if dangular_momentum_dt_per_moment_of_inertia {
            if  self.consider_tides || self.consider_rotational_flattening || ConsiderGeneralRelativity::Kidder1995 == self.consider_general_relativity {
                self.calculate_dangular_momentum_dt_per_moment_of_inertia();
            }
        }

        if accelerations {
            self.apply_acceleration_corrections();
        }

    }


    fn apply_acceleration_corrections(&mut self) {
        // Add the tidal+flattening+general relativity accelerations to the gravitational one (already computed)
        let (particles_left, particles_right) = self.particles[..self.n_particles].split_at_mut(self.tidal_host_particle_index);
        if let Some((tidal_host_particle, particles_right)) = particles_right.split_first_mut() {
            for particle in particles_left.iter_mut().chain(particles_right.iter_mut()) {
                if self.consider_tides {
                    particle.inertial_acceleration.x += particle.tidal_acceleration.x;
                    particle.inertial_acceleration.y += particle.tidal_acceleration.y;
                    particle.inertial_acceleration.z += particle.tidal_acceleration.z;
                }

                if self.consider_disk_interaction {
                    particle.inertial_acceleration.x += particle.disk_interaction_acceleration.x;
                    particle.inertial_acceleration.y += particle.disk_interaction_acceleration.y;
                    particle.inertial_acceleration.z += particle.disk_interaction_acceleration.z;
                }

                if self.consider_rotational_flattening {
                    particle.inertial_acceleration.x += particle.acceleration_induced_by_rotational_flattering.x;
                    particle.inertial_acceleration.y += particle.acceleration_induced_by_rotational_flattering.y;
                    particle.inertial_acceleration.z += particle.acceleration_induced_by_rotational_flattering.z;
                }

                if self.consider_general_relativity != ConsiderGeneralRelativity::None {
                    particle.inertial_acceleration.x += particle.general_relativity_acceleration.x;
                    particle.inertial_acceleration.y += particle.general_relativity_acceleration.y;
                    particle.inertial_acceleration.z += particle.general_relativity_acceleration.z;
                } 
            }
            if self.consider_tides {
                tidal_host_particle.inertial_acceleration.x += tidal_host_particle.tidal_acceleration.x;
                tidal_host_particle.inertial_acceleration.y += tidal_host_particle.tidal_acceleration.y;
                tidal_host_particle.inertial_acceleration.z += tidal_host_particle.tidal_acceleration.z;
            }

            if self.consider_disk_interaction {
                tidal_host_particle.inertial_acceleration.x += tidal_host_particle.disk_interaction_acceleration.x;
                tidal_host_particle.inertial_acceleration.y += tidal_host_particle.disk_interaction_acceleration.y;
                tidal_host_particle.inertial_acceleration.z += tidal_host_particle.disk_interaction_acceleration.z;
            }

            if self.consider_rotational_flattening {
                tidal_host_particle.inertial_acceleration.x += tidal_host_particle.acceleration_induced_by_rotational_flattering.x;
                tidal_host_particle.inertial_acceleration.y += tidal_host_particle.acceleration_induced_by_rotational_flattering.y;
                tidal_host_particle.inertial_acceleration.z += tidal_host_particle.acceleration_induced_by_rotational_flattering.z;
            }

            if self.consider_general_relativity != ConsiderGeneralRelativity::None {
                tidal_host_particle.inertial_acceleration.x += tidal_host_particle.general_relativity_acceleration.x;
                tidal_host_particle.inertial_acceleration.y += tidal_host_particle.general_relativity_acceleration.y;
                tidal_host_particle.inertial_acceleration.z += tidal_host_particle.general_relativity_acceleration.z;
            } 
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // TIDES
    fn calculate_dangular_momentum_dt_per_moment_of_inertia(&mut self) {
        for particle in self.particles[..self.n_particles].iter_mut() {
            // - Equation 25 from Bolmont et al. 2015
            particle.dangular_momentum_dt.x = particle.dangular_momentum_dt_due_to_tides.x + particle.dangular_momentum_dt_induced_by_rotational_flattening.x + particle.dangular_momentum_dt_due_to_general_relativity.x;
            particle.dangular_momentum_dt.y = particle.dangular_momentum_dt_due_to_tides.y + particle.dangular_momentum_dt_induced_by_rotational_flattening.y + particle.dangular_momentum_dt_due_to_general_relativity.y;
            particle.dangular_momentum_dt.z = particle.dangular_momentum_dt_due_to_tides.z + particle.dangular_momentum_dt_induced_by_rotational_flattening.z + particle.dangular_momentum_dt_due_to_general_relativity.z;
            //
            let factor = 1. / (particle.moment_of_inertia);
            particle.dangular_momentum_dt_per_moment_of_inertia.x = factor * particle.dangular_momentum_dt.x;
            particle.dangular_momentum_dt_per_moment_of_inertia.y = factor * particle.dangular_momentum_dt.y;
            particle.dangular_momentum_dt_per_moment_of_inertia.z = factor * particle.dangular_momentum_dt.z;
        }
    }

    pub fn calculate_denergy_dt(&mut self) {
        let (particles_left, particles_right) = self.particles[..self.n_particles].split_at_mut(self.tidal_host_particle_index);
        if let Some((_, particles_right)) = particles_right.split_first_mut() {
            for particle in particles_left.iter_mut().chain(particles_right.iter_mut()) {
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
                            - (particle.dangular_momentum_dt_due_to_tides.x*particle.spin.x + particle.dangular_momentum_dt_due_to_tides.y*particle.spin.y + particle.dangular_momentum_dt_due_to_tides.z*particle.spin.z);
            }
        }
    }


    /// 
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
}



fn get_center_of_mass_of_pair(center_of_mass_position: &mut Axes, center_of_mass_velocity: &mut Axes, center_of_mass_acceleration: &mut Axes, center_of_mass_mass: f64, particle: &Particle) -> f64 {
    center_of_mass_position.x      = center_of_mass_position.x*center_of_mass_mass + particle.position.x*particle.mass;
    center_of_mass_position.y      = center_of_mass_position.y*center_of_mass_mass + particle.position.y*particle.mass;
    center_of_mass_position.z      = center_of_mass_position.z*center_of_mass_mass + particle.position.z*particle.mass;
    center_of_mass_velocity.x      = center_of_mass_velocity.x*center_of_mass_mass + particle.velocity.x*particle.mass;
    center_of_mass_velocity.y      = center_of_mass_velocity.y*center_of_mass_mass + particle.velocity.y*particle.mass;
    center_of_mass_velocity.z      = center_of_mass_velocity.z*center_of_mass_mass + particle.velocity.z*particle.mass;
    center_of_mass_acceleration.x  = center_of_mass_acceleration.x*center_of_mass_mass + particle.acceleration.x*particle.mass;
    center_of_mass_acceleration.y  = center_of_mass_acceleration.y*center_of_mass_mass + particle.acceleration.y*particle.mass;
    center_of_mass_acceleration.z  = center_of_mass_acceleration.z*center_of_mass_mass + particle.acceleration.z*particle.mass;
    
    let new_center_of_mass_mass = center_of_mass_mass + particle.mass;
    if new_center_of_mass_mass > 0. {
        center_of_mass_position.x     /= new_center_of_mass_mass;
        center_of_mass_position.y     /= new_center_of_mass_mass;
        center_of_mass_position.z     /= new_center_of_mass_mass;
        center_of_mass_velocity.x     /= new_center_of_mass_mass;
        center_of_mass_velocity.y     /= new_center_of_mass_mass;
        center_of_mass_velocity.z     /= new_center_of_mass_mass;
        center_of_mass_acceleration.x /= new_center_of_mass_mass;
        center_of_mass_acceleration.y /= new_center_of_mass_mass;
        center_of_mass_acceleration.z /= new_center_of_mass_mass;
    }
    new_center_of_mass_mass
}

pub fn calculate_center_of_mass(particles: &Vec<Particle>) -> (Axes, Axes) {
    // Compute center of mass
    let mut center_of_mass_position = Axes{x:0., y:0., z:0.};
    let mut center_of_mass_velocity = Axes{x:0., y:0., z:0.};
    let mut center_of_mass_acceleration = Axes{x:0., y:0., z:0.};
    let mut center_of_mass_mass = 0.;

    for particle in particles.iter() {
        center_of_mass_mass = get_center_of_mass_of_pair(&mut center_of_mass_position, 
                                                                &mut center_of_mass_velocity, 
                                                                &mut center_of_mass_acceleration,
                                                                center_of_mass_mass,
                                                                &particle);
    }

    (center_of_mass_position, center_of_mass_velocity)
}
