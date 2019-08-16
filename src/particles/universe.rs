extern crate time;
use std::collections::HashMap;
use super::super::constants::{G, MAX_PARTICLES, MAX_DISTANCE_2};
use super::{Evolver, EvolutionType};
use super::{Particle};
use super::{Axes};
use super::{tides, flattening, general_relativity, evolver, wind, disk, common};
use super::{TidesEffect, RotationalFlatteningEffect, DiskEffect, WindEffect};
use super::{GeneralRelativityImplementation, GeneralRelativityEffect};

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
    pub consider_general_relativity: bool,
    pub general_relativity_implementation: GeneralRelativityImplementation, // Optimization: fast access to GR implementation for integrators
    pub most_massive_particle_index: usize, // Most massive particle
    general_relativity_host_particle_index: usize, // Central particle for general relativity effects
    tidal_host_particle_index: usize, // Particle that is the main one for tidal effects
    rotational_flattening_host_particle_index: usize, // Particle that is the main one for rotational flattenning effects
    disk_host_particle_index: usize, // Particle with disk
    all_same_particle_indices: bool, // Optimization
    star_planet_dependent_dissipation_factors : HashMap<usize, f64>, // Central body specific
    temporary_copied_particle_positions: [Axes; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_masses: [f64; MAX_PARTICLES], // For optimization purposes
    temporary_copied_particles_radiuses: [f64; MAX_PARTICLES], // For optimization purposes
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum IgnoreGravityTerms {
    None,
    WHFastOne,
    WHFastTwo,
}

impl Universe {
    pub fn new(initial_time: f64, time_limit: f64, mut particles: Vec<Particle>, consider_tides: bool, consider_rotational_flattening: bool, consider_disk_interaction: bool, consider_general_relativity: bool) -> Universe {

        let mut evolving_particles_exist = false;
        let mut wind_effects_exist = false;
        for particle in particles.iter(){
            if particle.wind.effect != WindEffect::None {
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

        let (most_massive_particle_index, tidal_host_particle_index, rotational_flattening_host_particle_index, general_relativity_host_particle_index, disk_host_particle_index) = find_indices(&particles, consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity);
        let all_same_particle_indices = ((consider_tides && most_massive_particle_index == tidal_host_particle_index) || !consider_tides)
            && ((consider_rotational_flattening && most_massive_particle_index == rotational_flattening_host_particle_index) || !consider_rotational_flattening)
            && ((consider_disk_interaction && most_massive_particle_index == disk_host_particle_index) || !consider_disk_interaction)
            && ((consider_general_relativity && most_massive_particle_index == general_relativity_host_particle_index) || !consider_general_relativity);

        // Initialize general relativity factor
        let mut general_relativity_implementation = GeneralRelativityImplementation::Kidder1995;
        if consider_general_relativity {
            let (particles_left, particles_right) = particles.split_at_mut(general_relativity_host_particle_index);
            if let Some((general_relativity_host_particle, particles_right)) = particles_right.split_first_mut() {
                if let GeneralRelativityEffect::CentralBody(implementation) = general_relativity_host_particle.general_relativity.effect {
                    general_relativity_implementation = implementation;
                    let local_copy_star_mass_g = general_relativity_host_particle.mass_g;
                    for particle in particles_left.iter_mut().chain(particles_right.iter_mut()) {
                        particle.general_relativity.parameters.internal.factor =  local_copy_star_mass_g*particle.mass_g / (local_copy_star_mass_g + particle.mass_g).powi(2)
                    }
                }
            }
        }

        // Re-compute inertial positions/velocities (if the user introduced heliocentric
        // positions/velocities, they are transformed to barycentric)
        let (center_of_mass_position, center_of_mass_velocity) = calculate_center_of_mass(&particles);
        for particle in particles.iter_mut() {
            particle.inertial_position.x = particle.heliocentric_position.x - center_of_mass_position.x;
            particle.inertial_position.y = particle.heliocentric_position.y - center_of_mass_position.y;
            particle.inertial_position.z = particle.heliocentric_position.z - center_of_mass_position.z;
            particle.inertial_velocity.x = particle.heliocentric_velocity.x - center_of_mass_velocity.x;
            particle.inertial_velocity.y = particle.heliocentric_velocity.y - center_of_mass_velocity.y;
            particle.inertial_velocity.z = particle.heliocentric_velocity.z - center_of_mass_velocity.z;
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
                    general_relativity_implementation: general_relativity_implementation,
                    most_massive_particle_index: most_massive_particle_index,
                    tidal_host_particle_index: tidal_host_particle_index,
                    rotational_flattening_host_particle_index: rotational_flattening_host_particle_index,
                    general_relativity_host_particle_index: general_relativity_host_particle_index,
                    disk_host_particle_index: disk_host_particle_index,
                    all_same_particle_indices: all_same_particle_indices,
                    star_planet_dependent_dissipation_factors:HashMap::new(),
                    temporary_copied_particle_positions: temporary_copied_particle_positions,
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
                    let roche_radius;
                    if particle_a.mass > *particle_b_mass {
                        // particle a is the most massive of both
                        roche_radius = (particle_a.radius/0.462)*(particle_a.mass/particle_b_mass).powf(-1./3.);
                    } else {
                        // particle b is the most massive of both
                        roche_radius = (particle_b_radius/0.462)*(particle_b_mass/particle_a.mass).powf(-1./3.);
                    }
                    if distance_2 <= roche_radius.powi(2) {
                        println!("\n");
                        panic!("[PANIC {} UTC] Particle {} was destroyed by particle {} due to close encounter!", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), i, j);
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

                if ignore_terms == IgnoreGravityTerms::WHFastOne && 
                        ((i == self.most_massive_particle_index && self.most_massive_particle_index == 0 && j == 1)
                        || (i == self.most_massive_particle_index && self.most_massive_particle_index > 0 && j == 0)
                        || (j == self.most_massive_particle_index && self.most_massive_particle_index == 0 && i == 1)
                        || (j == self.most_massive_particle_index && self.most_massive_particle_index > 0 && i == 0))
                    {
                    // For WHFast Jacobi coordinates, 
                    // ignore interaction between central body and first planet
                    continue
                }
                if ignore_terms == IgnoreGravityTerms::WHFastTwo && (i == self.most_massive_particle_index || j == self.most_massive_particle_index) {
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

    pub fn inertial_to_heliocentric(&mut self) {
        let (particles_left, particles_right) = self.particles.split_at_mut(self.most_massive_particle_index);
        if let Some((host_particle, particles_right)) = particles_right.split_first_mut() {
            for particle in particles_left.iter_mut().chain(particles_right.iter_mut()) {
                particle.heliocentric_position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                particle.heliocentric_position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                particle.heliocentric_position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                particle.heliocentric_velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.heliocentric_velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.heliocentric_velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.heliocentric_distance = (particle.heliocentric_position.x.powi(2) 
                                                                    + particle.heliocentric_position.y.powi(2)
                                                                    + particle.heliocentric_position.z.powi(2)).sqrt();
                particle.heliocentric_radial_velocity = (particle.heliocentric_position.x*particle.heliocentric_velocity.x +
                                                                        particle.heliocentric_position.y*particle.heliocentric_velocity.y +
                                                                        particle.heliocentric_position.z*particle.heliocentric_velocity.z) 
                                                                        / particle.heliocentric_distance;
                particle.heliocentric_norm_velocity_vector_2 = (particle.heliocentric_velocity.x - host_particle.heliocentric_velocity.x).powi(2) 
                                                                + (particle.heliocentric_velocity.y - host_particle.heliocentric_velocity.y).powi(2)
                                                                + (particle.heliocentric_velocity.z - host_particle.heliocentric_velocity.z).powi(2);
                particle.heliocentric_norm_velocity_vector = particle.heliocentric_norm_velocity_vector_2.sqrt();
            }
            host_particle.heliocentric_position.x = 0.;
            host_particle.heliocentric_position.y = 0.;
            host_particle.heliocentric_position.z = 0.;
            host_particle.heliocentric_velocity.x = 0.;
            host_particle.heliocentric_velocity.y = 0.;
            host_particle.heliocentric_velocity.z = 0.;
            host_particle.heliocentric_distance = 0.;
            host_particle.heliocentric_radial_velocity = 0.;
            host_particle.heliocentric_norm_velocity_vector_2 = 0.;
            host_particle.heliocentric_norm_velocity_vector = 0.;
        }
    }

    pub fn initialize(&mut self, dangular_momentum_dt_per_moment_of_inertia: bool, accelerations: bool) {
        let (particles, _) = self.particles.split_at_mut(self.n_particles);
        //
        let (mut particles_left, particles_right) = particles.split_at_mut(self.most_massive_particle_index);
        if let Some((mut host_particle, mut particles_right)) = particles_right.split_first_mut() {
            if self.tidal_host_particle_index == self.most_massive_particle_index {
                tides::copy_heliocentric_coordinates(&mut host_particle, &mut particles_left, &mut particles_right);
                tides::initialize(&mut host_particle, &mut particles_left, &mut particles_right);
            }
            if self.rotational_flattening_host_particle_index == self.most_massive_particle_index {
                flattening::copy_heliocentric_coordinates(&mut host_particle, &mut particles_left, &mut particles_right);
                flattening::initialize(&mut host_particle, &mut particles_left, &mut particles_right);
            }
            if self.general_relativity_host_particle_index == self.most_massive_particle_index {
                general_relativity::copy_heliocentric_coordinates(&mut host_particle, &mut particles_left, &mut particles_right);
                general_relativity::initialize(&mut host_particle, &mut particles_left, &mut particles_right);
            }
            if self.disk_host_particle_index == self.most_massive_particle_index {
                disk::copy_heliocentric_coordinates(&mut host_particle, &mut particles_left, &mut particles_right);
                disk::initialize(&mut host_particle, &mut particles_left, &mut particles_right);
            }
        }
        //
        if !self.all_same_particle_indices {
            let initialize_tides = (dangular_momentum_dt_per_moment_of_inertia && self.consider_tides) || (accelerations && self.consider_tides);
            let initialize_rotational_flattening = (dangular_momentum_dt_per_moment_of_inertia && self.consider_rotational_flattening) || (accelerations && self.consider_rotational_flattening);
            let initialize_general_relativity = (dangular_momentum_dt_per_moment_of_inertia && self.consider_general_relativity && self.general_relativity_implementation == GeneralRelativityImplementation::Kidder1995) || (accelerations && self.consider_general_relativity);
            let initialize_disk = accelerations && self.consider_disk_interaction;
            //
            if initialize_tides && self.tidal_host_particle_index != self.most_massive_particle_index {
                let (mut particles_left, particles_right) = particles.split_at_mut(self.tidal_host_particle_index);
                if let Some((mut tidal_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                    tides::inertial_to_heliocentric_coordinates(&mut tidal_host_particle, &mut particles_left, &mut particles_right);
                    tides::initialize(&mut tidal_host_particle, &mut particles_left, &mut particles_right);
                }
            }
            if initialize_rotational_flattening && self.rotational_flattening_host_particle_index != self.most_massive_particle_index {
                let (mut particles_left, particles_right) = particles.split_at_mut(self.rotational_flattening_host_particle_index);
                if let Some((mut rotational_flattening_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                    flattening::inertial_to_heliocentric_coordinates(&mut rotational_flattening_host_particle, &mut particles_left, &mut particles_right);
                    flattening::initialize(&mut rotational_flattening_host_particle, &mut particles_left, &mut particles_right);
                }
            }
            if initialize_general_relativity && self.general_relativity_host_particle_index != self.most_massive_particle_index {
                let (mut particles_left, particles_right) = particles.split_at_mut(self.general_relativity_host_particle_index);
                if let Some((mut general_relativity_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                    general_relativity::inertial_to_heliocentric_coordinates(&mut general_relativity_host_particle, &mut particles_left, &mut particles_right);
                    general_relativity::initialize(&mut general_relativity_host_particle, &mut particles_left, &mut particles_right);
                }
            }
            if initialize_disk && self.disk_host_particle_index != self.most_massive_particle_index {
                let (mut particles_left, particles_right) = particles.split_at_mut(self.disk_host_particle_index);
                if let Some((mut disk_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                    disk::inertial_to_heliocentric_coordinates(&mut disk_host_particle, &mut particles_left, &mut particles_right);
                    disk::initialize(&mut disk_host_particle, &mut particles_left, &mut particles_right);
                }
            }
        }
    }

    pub fn calculate_additional_effects(&mut self, current_time: f64, evolution: bool, dangular_momentum_dt_per_moment_of_inertia: bool, accelerations: bool, ignored_gravity_terms: IgnoreGravityTerms) {
        self.initialize(dangular_momentum_dt_per_moment_of_inertia, accelerations);

        let (mut particles, _) = self.particles.split_at_mut(self.n_particles);
        if (evolution && self.evolving_particles_exist) || 
            ((dangular_momentum_dt_per_moment_of_inertia || accelerations) && self.consider_rotational_flattening) {
            common::calculate_norm_spin(&mut particles); // Needed for rotational flattening (torque and accelerations) and evolution
        }
       
        if evolution && self.evolving_particles_exist {
            evolver::calculate_particles_evolving_quantities(current_time, &mut particles, &mut self.particles_evolvers);
        }

        let (mut particles, _) = self.particles.split_at_mut(self.n_particles);
        if dangular_momentum_dt_per_moment_of_inertia && self.wind_effects_exist {
            wind::calculate_wind_factor(&mut particles);
        }

        if self.consider_tides || self.consider_rotational_flattening {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.tidal_host_particle_index);
            if let Some((mut tidal_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                if (dangular_momentum_dt_per_moment_of_inertia && (self.consider_tides || self.consider_rotational_flattening)) ||
                    (accelerations && (self.consider_tides || self.consider_disk_interaction || self.consider_rotational_flattening || 
                                       self.consider_general_relativity)) {

                    {
                        //// calculate_orthogonal_components

                        tides::calculate_planet_dependent_dissipation_factors(&mut tidal_host_particle, &mut particles_left, &mut particles_right, &mut self.star_planet_dependent_dissipation_factors); // Needed by calculate_orthogonal_component_of_the_tidal_force and calculate_orthogonal_component_of_the_tidal_force if BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true)
                        
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
        if accelerations && self.consider_general_relativity {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.general_relativity_host_particle_index);
            if let Some((mut general_relativity_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) = general_relativity_host_particle.general_relativity.effect {
                    match general_relativity_implementation {
                        GeneralRelativityImplementation::Kidder1995 => {
                            general_relativity::calculate_kidder1995_general_relativity_acceleration(&mut general_relativity_host_particle, 
                                                                                                     &mut particles_left, &mut particles_right);
                        },
                        GeneralRelativityImplementation::Anderson1975 => {
                            general_relativity::calculate_anderson1975_general_relativity_acceleration(&mut general_relativity_host_particle, 
                                                                                                       &mut particles_left, &mut particles_right, 
                                                                                                       ignored_gravity_terms);
                        },
                        GeneralRelativityImplementation::Newhall1983 => {
                            general_relativity::calculate_newhall1983_general_relativity_acceleration(&mut general_relativity_host_particle, 
                                                                                                      &mut particles_left, &mut particles_right, 
                                                                                                      ignored_gravity_terms);
                        },
                    }
                }
            }
        }
        if accelerations && self.consider_disk_interaction {
            let (mut particles_left, particles_right) = particles.split_at_mut(self.disk_host_particle_index);
            if let Some((mut disk_host_particle, mut particles_right)) = particles_right.split_first_mut() {
                if let DiskEffect::CentralBody(_) = disk_host_particle.disk.effect {
                    disk::calculate_disk_interaction_acceleration(current_time, &mut disk_host_particle, &mut particles_left, &mut particles_right);
                }
            }
        }
        if dangular_momentum_dt_per_moment_of_inertia {
            if  self.consider_tides || self.consider_rotational_flattening || (self.consider_general_relativity && self.general_relativity_implementation == GeneralRelativityImplementation::Kidder1995) {
                self.calculate_dangular_momentum_dt_per_moment_of_inertia();
            }
        }

        if accelerations {
            self.apply_acceleration_corrections();
            //println!("{:?}", self.particles[1].tides.parameters.output.acceleration);
            //println!("{:?}", self.particles[1].rotational_flattening.parameters.output.acceleration);
            //println!("{:?}", self.particles[1].general_relativity.parameters.output.acceleration);
            //println!("{:?}", self.particles[1].disk.parameters.output.acceleration);
            //panic!("*** {:?}", self.particles[1].inertial_acceleration);
        }

    }


    fn apply_acceleration_corrections(&mut self) {
        // Add the tidal+flattening+general relativity accelerations to the gravitational one (already computed)
        for particle in self.particles[..self.n_particles].iter_mut() {
            if self.consider_tides {
                particle.inertial_acceleration.x += particle.tides.parameters.output.acceleration.x;
                particle.inertial_acceleration.y += particle.tides.parameters.output.acceleration.y;
                particle.inertial_acceleration.z += particle.tides.parameters.output.acceleration.z;
            }

            if self.consider_disk_interaction {
                particle.inertial_acceleration.x += particle.disk.parameters.output.acceleration.x;
                particle.inertial_acceleration.y += particle.disk.parameters.output.acceleration.y;
                particle.inertial_acceleration.z += particle.disk.parameters.output.acceleration.z;
            }

            if self.consider_rotational_flattening {
                particle.inertial_acceleration.x += particle.rotational_flattening.parameters.output.acceleration.x;
                particle.inertial_acceleration.y += particle.rotational_flattening.parameters.output.acceleration.y;
                particle.inertial_acceleration.z += particle.rotational_flattening.parameters.output.acceleration.z;
            }

            if self.consider_general_relativity {
                particle.inertial_acceleration.x += particle.general_relativity.parameters.output.acceleration.x;
                particle.inertial_acceleration.y += particle.general_relativity.parameters.output.acceleration.y;
                particle.inertial_acceleration.z += particle.general_relativity.parameters.output.acceleration.z;
            } 
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // TIDES
    fn calculate_dangular_momentum_dt_per_moment_of_inertia(&mut self) {
        for particle in self.particles[..self.n_particles].iter_mut() {
            // - Equation 25 from Bolmont et al. 2015
            particle.dangular_momentum_dt.x = particle.tides.parameters.output.dangular_momentum_dt.x + particle.rotational_flattening.parameters.output.dangular_momentum_dt.x + particle.general_relativity.parameters.output.dangular_momentum_dt.x;
            particle.dangular_momentum_dt.y = particle.tides.parameters.output.dangular_momentum_dt.y + particle.rotational_flattening.parameters.output.dangular_momentum_dt.y + particle.general_relativity.parameters.output.dangular_momentum_dt.y;
            particle.dangular_momentum_dt.z = particle.tides.parameters.output.dangular_momentum_dt.z + particle.rotational_flattening.parameters.output.dangular_momentum_dt.z + particle.general_relativity.parameters.output.dangular_momentum_dt.z;
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
            tides::calculate_denergy_dt(particles_left, particles_right);
        }
    }


    /// 
    pub fn compute_total_energy(&self) -> f64 {
        let mut e_kin = 0.;
        let mut e_pot = 0.;
        let e_offset = 0.; // Energy offset due to collisions and ejections

        // Kinectic energy
        for particle in self.particles[..self.n_particles].iter() {
            e_kin += 0.5 * particle.mass * (particle.heliocentric_velocity.x.powi(2) + particle.heliocentric_velocity.y.powi(2) + particle.heliocentric_velocity.z.powi(2));
        }
        // Gravitationl potential energy
        for (i, particle_a) in self.particles[..self.n_particles].iter().enumerate() {
            for particle_b in self.particles[i+1..self.n_particles].iter() {
                let dx = particle_a.heliocentric_position.x - particle_b.heliocentric_position.x;
                let dy = particle_a.heliocentric_position.y - particle_b.heliocentric_position.y;
                let dz = particle_a.heliocentric_position.z - particle_b.heliocentric_position.z;
                e_pot -= particle_b.mass_g*particle_a.mass/(dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
            }
        }
        
        e_kin + e_pot + e_offset
    }

    pub fn compute_total_angular_momentum(&self) -> f64 {
        let mut total_angular_momentum = Axes{x:0., y:0., z:0.}; // L
        for particle in self.particles[..self.n_particles].iter() {
            total_angular_momentum.x += particle.mass*(particle.heliocentric_position.y*particle.heliocentric_velocity.z - particle.heliocentric_position.z*particle.heliocentric_velocity.y);
            total_angular_momentum.y += particle.mass*(particle.heliocentric_position.z*particle.heliocentric_velocity.x - particle.heliocentric_position.x*particle.heliocentric_velocity.z);
            total_angular_momentum.z += particle.mass*(particle.heliocentric_position.x*particle.heliocentric_velocity.y - particle.heliocentric_position.y*particle.heliocentric_velocity.x);
        }
        let total_angular_momentum = (total_angular_momentum.x.powi(2) + total_angular_momentum.y.powi(2) + total_angular_momentum.z.powi(2)).sqrt();
        total_angular_momentum
    }
}



fn get_center_of_mass_of_pair(center_of_mass_position: &mut Axes, center_of_mass_velocity: &mut Axes, center_of_mass_mass: f64, particle: &Particle) -> f64 {
    center_of_mass_position.x      = center_of_mass_position.x*center_of_mass_mass + particle.heliocentric_position.x*particle.mass;
    center_of_mass_position.y      = center_of_mass_position.y*center_of_mass_mass + particle.heliocentric_position.y*particle.mass;
    center_of_mass_position.z      = center_of_mass_position.z*center_of_mass_mass + particle.heliocentric_position.z*particle.mass;
    center_of_mass_velocity.x      = center_of_mass_velocity.x*center_of_mass_mass + particle.heliocentric_velocity.x*particle.mass;
    center_of_mass_velocity.y      = center_of_mass_velocity.y*center_of_mass_mass + particle.heliocentric_velocity.y*particle.mass;
    center_of_mass_velocity.z      = center_of_mass_velocity.z*center_of_mass_mass + particle.heliocentric_velocity.z*particle.mass;
    
    let new_center_of_mass_mass = center_of_mass_mass + particle.mass;
    if new_center_of_mass_mass > 0. {
        center_of_mass_position.x     /= new_center_of_mass_mass;
        center_of_mass_position.y     /= new_center_of_mass_mass;
        center_of_mass_position.z     /= new_center_of_mass_mass;
        center_of_mass_velocity.x     /= new_center_of_mass_mass;
        center_of_mass_velocity.y     /= new_center_of_mass_mass;
        center_of_mass_velocity.z     /= new_center_of_mass_mass;
    }
    new_center_of_mass_mass
}

pub fn calculate_center_of_mass(particles: &Vec<Particle>) -> (Axes, Axes) {
    // Compute center of mass
    let mut center_of_mass_position = Axes{x:0., y:0., z:0.};
    let mut center_of_mass_velocity = Axes{x:0., y:0., z:0.};
    let mut center_of_mass_mass = 0.;

    for particle in particles.iter() {
        center_of_mass_mass = get_center_of_mass_of_pair(&mut center_of_mass_position, 
                                                                &mut center_of_mass_velocity, 
                                                                center_of_mass_mass,
                                                                &particle);
    }

    (center_of_mass_position, center_of_mass_velocity)
}

fn find_indices(particles: &Vec<Particle>, consider_tides: bool, consider_rotational_flattening: bool, consider_disk_interaction: bool, consider_general_relativity: bool) -> (usize, usize, usize, usize, usize) {
    // Most massive particle
    let mut most_massive_particle_index = MAX_PARTICLES+1;
    let mut max_mass_found = 0.;
    for (i, particle) in particles.iter().enumerate() {
        if particle.mass > max_mass_found {
            max_mass_found = particle.mass;
            most_massive_particle_index = i;
        }
    }

    // Particle that is the main one for general relativity effects and WHFast symplectic integration
    let mut general_relativity_host_particle_index = MAX_PARTICLES+1;
    if consider_general_relativity {
        for (i, particle) in particles.iter().enumerate() {
            if let GeneralRelativityEffect::CentralBody(_implementation) = particle.general_relativity.effect {
                if general_relativity_host_particle_index == MAX_PARTICLES+1 {
                    general_relativity_host_particle_index = i;
                } else {
                    panic!("Only one central body is allowed for general relativity effects!");
                }
            }
        }
        if most_massive_particle_index != general_relativity_host_particle_index {
            panic!("The central body for General Relativity should be the most massive one!");
        }
    } else {
        general_relativity_host_particle_index = most_massive_particle_index;
    }

    
    // Particle that is the main one for tidal effects
    let mut tidal_host_particle_index = MAX_PARTICLES+1;
    if consider_tides {
        for (i, particle) in particles.iter().enumerate() {
            if let TidesEffect::CentralBody = particle.tides.effect {
                if tidal_host_particle_index == MAX_PARTICLES+1 {
                    tidal_host_particle_index = i;
                } else {
                    panic!("Only one central body is allowed for tidal effects!");
                }
            }
        }
    }

    // Particle that is the main one for rotational flatting effects
    let mut rotational_flattening_host_particle_index = MAX_PARTICLES+1;
    if consider_rotational_flattening {
        for (i, particle) in particles.iter().enumerate() {
            if let RotationalFlatteningEffect::CentralBody = particle.rotational_flattening.effect {
                if rotational_flattening_host_particle_index == MAX_PARTICLES+1 {
                    rotational_flattening_host_particle_index = i;
                } else {
                    panic!("Only one central body is allowed for rotational flattening effects!");
                }
            }
        }
    }

    if consider_tides && consider_rotational_flattening && tidal_host_particle_index != rotational_flattening_host_particle_index {
        panic!("The central body for tidal & rotational flattening effects needs to be the same!");
    } else if !consider_tides && consider_rotational_flattening {
        tidal_host_particle_index = rotational_flattening_host_particle_index;
    }

    // Find the position of the particle that has the disk if any
    let mut disk_host_particle_index = MAX_PARTICLES+1;
    if consider_disk_interaction {
        for (i, particle) in particles.iter().enumerate() {
            if let DiskEffect::CentralBody( disk ) = particle.disk.effect {
                if disk.inner_edge_distance != 0. || disk.outer_edge_distance != 0. {
                    if disk_host_particle_index == MAX_PARTICLES+1 {
                        disk_host_particle_index = i;
                    } else {
                        panic!("Only one central body with a disk is allowed!");
                    }
                }
            }
        }
    }

    (most_massive_particle_index, tidal_host_particle_index, rotational_flattening_host_particle_index, general_relativity_host_particle_index, disk_host_particle_index)
}
