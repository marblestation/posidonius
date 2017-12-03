extern crate time;
use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use super::Integrator;
use super::super::particles::Universe;
use super::super::particles::IgnoreGravityTerms;
use super::super::particles::Axes;
use super::output::{write_historic_snapshot};
use rustc_serialize::json;
use bincode::rustc_serialize::{decode_from};
use bincode::SizeLimit;
use std::io::{BufReader};
use std::io::Read;
use std::path::Path;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

/// LeapFrog is a second order symplectic integrator
/// http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/9710043
///
/// http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1411.6671
/// Time-reversible, symplectic integrators should in principle conserve energy to a better level
/// than non-symplectic integrators, since there is no drift present in the energy error.
/// Therefore, by using a symplectic integrator, the number of simulations with large energy error
/// could be reduced. Using a Leapfrog integrator with constant time-steps, we tested this
/// assumption and we find that for resonant 3-body interactions, it is challenging to obtain
/// accurate solutions. The main reason is that, contrary to regular systems like, for example, the
/// Solar System, resonant 3-body interactions often include very close encounters, which need a
/// very small time-step size to be resolved accurately. 
///
///https://arxiv.org/pdf/1110.4876v2.pdf
///Leap-frog is a second-order accurate and symplectic integrator
///for non-rotating frames. Here, the Hamiltonian is split into the
///kinetic part and the potential part. Both
///the drift and kick sub-steps are simple Euler steps. First the positions
///of all particles are advanced for half a time-step while keeping
///the velocities fixed. Then the velocities are advanced for one
///time-step while keeping the positions fixed. In the last sub-step
///the velocities are again advanced for half a time-step. L


#[derive(Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub struct LeapFrog {
    time_step: f64,
    half_time_step: f64,
    pub universe: Universe,
    current_time: f64,
    current_iteration: u32,
    recovery_snapshot_period: f64,
    historic_snapshot_period: f64,
    last_recovery_snapshot_time: f64,
    last_historic_snapshot_time: f64,
    pub n_historic_snapshots: usize,
    pub hash: u64,
}

impl Hash for LeapFrog {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash only works with certain types (for instance, it does not work with f64 values)
        // thus we convert the whole integrator to a string thanks to the debug trait
        // and we hash that value
        format!("{:?}", self).hash(state); 
    }
}

impl LeapFrog {
    pub fn new(time_step: f64, recovery_snapshot_period: f64, historic_snapshot_period: f64, universe: Universe) -> LeapFrog {
        let mut universe_integrator = LeapFrog {
                    time_step:time_step,
                    half_time_step:0.5*time_step,
                    recovery_snapshot_period:recovery_snapshot_period,
                    historic_snapshot_period:historic_snapshot_period,
                    last_recovery_snapshot_time:-1.,
                    last_historic_snapshot_time:-1.,
                    n_historic_snapshots:0,
                    hash: 0,
                    universe:universe,
                    current_time:0.,
                    current_iteration:0,
                    };
        // Initialize physical values
        let current_time = 0.;
        universe_integrator.universe.calculate_norm_spin(); // Needed for evolution
        universe_integrator.universe.calculate_particles_evolving_quantities(current_time); // Make sure we start with the good initial values
        universe_integrator
    }

    pub fn restore_snapshot(universe_integrator_snapshot_path: &Path, verify_integrity: bool) -> Result<LeapFrog, String> {
        let mut universe_integrator: LeapFrog;
        if universe_integrator_snapshot_path.exists() {
            // Open the path in read-only mode, returns `io::Result<File>`
            let mut snapshot_file = match File::open(&universe_integrator_snapshot_path) {
                // The `description` method of `io::Error` returns a string that
                // describes the error
                Err(why) => return Err(format!("Couldn't open {}: {}", universe_integrator_snapshot_path.display(), why)),
                Ok(file) => file,
            };


            if universe_integrator_snapshot_path.extension().unwrap() == "json" { 
                //// Deserialize using `json::decode`
                let mut json_encoded = String::new();
                match snapshot_file.read_to_string(&mut json_encoded) {
                    Err(why) => return Err(format!("Couldn't read {}: {}", universe_integrator_snapshot_path.display(), why)),
                    Ok(_) => {}
                }
                universe_integrator = json::decode(&json_encoded).unwrap();
            } else {
                let mut reader = BufReader::new(snapshot_file);
                universe_integrator = decode_from(&mut reader, SizeLimit::Infinite).unwrap();
            }
            if universe_integrator.hash == 0 {
                if verify_integrity && universe_integrator.current_time != 0. {
                    panic!("[PANIC {} UTC] File '{}' has a zeroed hash (i.e., new simulation) but a current time different from zero ({:})", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display(), universe_integrator.current_time)
                }
                println!("[INFO {} UTC] Created new simulation based on '{}'", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display());
                // Initialize physical values
                let current_time = 0.;
                universe_integrator.universe.calculate_norm_spin(); // Needed for evolution
                universe_integrator.universe.calculate_particles_evolving_quantities(current_time); // Make sure we start with the good initial values
            } else {
                // Verify hash for this universe at this moment of time
                let mut s = DefaultHasher::new();
                let restored_hash = universe_integrator.hash;
                universe_integrator.hash = 0;
                universe_integrator.hash(&mut s);
                let computed_hash = s.finish();
                if verify_integrity && restored_hash != computed_hash {
                    panic!("[PANIC {} UTC] File '{}' seems corrupted because computed hash '{}' does not match restored hash '{}'", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display(), computed_hash, restored_hash)
                }
                universe_integrator.hash = restored_hash;
                println!("[INFO {} UTC] Restored previous simulation from '{}'", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display());
            }
            return Ok(universe_integrator);
        } else {
            return Err(format!("File does not exist"));
        }
    }

}

impl Integrator for LeapFrog {


    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            self.inertial_to_heliocentric_posvelacc();
            if self.universe.consider_tides {
                self.universe.calculate_denergy_dt();
            }
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.n_historic_snapshots as f64*self.historic_snapshot_period; // Instead of self.current_time to avoid small deviations
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            if ! silent_mode {
                print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
                let _ = std::io::stdout().flush();
            }
        }

        let ignore_gravity_terms = IgnoreGravityTerms::None;
        let ignored_gravity_terms = ignore_gravity_terms;

        self.inertial_to_heliocentric_posvelacc();
        // Calculate non-gravity accelerations.
        let evolution = true;
        let dspin_dt = true;
        let accelerations = false;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations, ignored_gravity_terms);
        self.heliocentric_to_inertial_posvelacc();


        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        self.universe.gravity_calculate_acceleration(ignore_gravity_terms);
        self.inertial_to_heliocentric_posvelacc();

        // Calculate non-gravity accelerations.
        let evolution = true;
        let dspin_dt = true;
        let accelerations = true;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations, ignored_gravity_terms);
        self.heliocentric_to_inertial_posvelacc();


        // A 'DKD'-like integrator will do the 'KD' part.
        self.integrator_part2();
        self.current_iteration += 1;

        // Return
        if self.current_time+self.time_step > self.universe.time_limit {
            Err("reached maximum time limit.".to_string())
        } else {
            Ok(first_snapshot_trigger || recovery_snapshot_time_trigger)
        }
    }

    fn prepare_for_recovery_snapshot(&mut self, universe_history_writer: &mut BufWriter<File>) {
        self.last_recovery_snapshot_time = self.current_time;
        universe_history_writer.flush().unwrap();
        // Compute hash for this universe at this moment of time
        let mut s = DefaultHasher::new();
        self.hash = 0;
        self.hash(&mut s);
        self.hash = s.finish();
    }

}

impl LeapFrog {
    // Leapfrog integrator (Drift-Kick-Drift)
    // for non-rotating frame.
    #[allow(dead_code)]
    fn integrator_part1(&mut self) {
        for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            particle.inertial_position.x += self.half_time_step * particle.inertial_velocity.x;
            particle.inertial_position.y += self.half_time_step * particle.inertial_velocity.y;
            particle.inertial_position.z += self.half_time_step * particle.inertial_velocity.z;

            if particle.moment_of_inertia_ratio != 1. {
                particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.half_time_step * particle.dspin_dt.x;
                particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.half_time_step * particle.dspin_dt.y;
                particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.half_time_step * particle.dspin_dt.z;
            } else {
                particle.spin.x = particle.spin.x + self.half_time_step * particle.dspin_dt.x;
                particle.spin.y = particle.spin.y + self.half_time_step * particle.dspin_dt.y;
                particle.spin.z = particle.spin.z + self.half_time_step * particle.dspin_dt.z;
            }
            if particle.wind_factor != 0. {
                // TODO: Verify wind factor
                particle.spin.x += self.half_time_step * particle.wind_factor * particle.spin.x;
                particle.spin.y += self.half_time_step * particle.wind_factor * particle.spin.y;
                particle.spin.z += self.half_time_step * particle.wind_factor * particle.spin.z;
            }
        }
        self.current_time += self.half_time_step;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            particle.inertial_velocity.x += self.time_step * particle.inertial_acceleration.x;
            particle.inertial_velocity.y += self.time_step * particle.inertial_acceleration.y;
            particle.inertial_velocity.z += self.time_step * particle.inertial_acceleration.z;
            particle.inertial_position.x += self.half_time_step * particle.inertial_velocity.x;
            particle.inertial_position.y += self.half_time_step * particle.inertial_velocity.y;
            particle.inertial_position.z += self.half_time_step * particle.inertial_velocity.z;

            if particle.moment_of_inertia_ratio != 1. {
                particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.half_time_step * particle.dspin_dt.x;
                particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.half_time_step * particle.dspin_dt.y;
                particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.half_time_step * particle.dspin_dt.z;
            } else {
                particle.spin.x = particle.spin.x + self.half_time_step * particle.dspin_dt.x;
                particle.spin.y = particle.spin.y + self.half_time_step * particle.dspin_dt.y;
                particle.spin.z = particle.spin.z + self.half_time_step * particle.dspin_dt.z;
            }
            if particle.wind_factor != 0. {
                // TODO: Verify wind factor
                particle.spin.x += self.half_time_step * particle.wind_factor * particle.spin.x;
                particle.spin.y += self.half_time_step * particle.wind_factor * particle.spin.y;
                particle.spin.z += self.half_time_step * particle.wind_factor * particle.spin.z;
            }
        }
        self.current_time += self.half_time_step;
    }

    fn heliocentric_to_inertial_posvelacc(&mut self) {
        // Formulation:
        //
        // inertial_position = heliocentric_position - (1/total_mass) * sum(particle.mass * particle.position)
        // inertial_velocity = heliocentric_velocity - (1/total_mass) * sum(particle.mass * particle.velocity)
        // inertial_acceleration = heliocentric_acceleration - (1/total_mass) * sum(particle.mass * particle.acceleration)
        //
        // Compute center of mass
        let mut reference_position = Axes{x:0., y:0., z:0.};
        let mut reference_velocity = Axes{x:0., y:0., z:0.};
        let mut reference_acceleration = Axes{x:0., y:0., z:0.};
        let mut total_mass = 0.;

        for particle in self.universe.particles[..self.universe.n_particles].iter() {
            reference_position.x      += particle.position.x*particle.mass;
            reference_position.y      += particle.position.y*particle.mass;
            reference_position.z      += particle.position.z*particle.mass;
            reference_velocity.x      += particle.velocity.x*particle.mass;
            reference_velocity.y      += particle.velocity.y*particle.mass;
            reference_velocity.z      += particle.velocity.z*particle.mass;
            reference_acceleration.x  += particle.acceleration.x*particle.mass;
            reference_acceleration.y  += particle.acceleration.y*particle.mass;
            reference_acceleration.z  += particle.acceleration.z*particle.mass;
            total_mass += particle.mass;
        }
        reference_position.x /= total_mass;
        reference_position.y /= total_mass;
        reference_position.z /= total_mass;
        reference_velocity.x /= total_mass;
        reference_velocity.y /= total_mass;
        reference_velocity.z /= total_mass;
        reference_acceleration.x /= total_mass;
        reference_acceleration.y /= total_mass;
        reference_acceleration.z /= total_mass;

        // Convert heliocentric to inertial
        for particle in self.universe.particles[0..self.universe.n_particles].iter_mut() {
            particle.inertial_position.x = particle.position.x - reference_position.x;
            particle.inertial_position.y = particle.position.y - reference_position.y;
            particle.inertial_position.z = particle.position.z - reference_position.z;
            particle.inertial_velocity.x = particle.velocity.x - reference_velocity.x;
            particle.inertial_velocity.y = particle.velocity.y - reference_velocity.y;
            particle.inertial_velocity.z = particle.velocity.z - reference_velocity.z;
            particle.inertial_acceleration.x = particle.acceleration.x - reference_acceleration.x;
            particle.inertial_acceleration.y = particle.acceleration.y - reference_acceleration.y;
            particle.inertial_acceleration.z = particle.acceleration.z - reference_acceleration.z;
        }
    }

    fn inertial_to_heliocentric_posvelacc(&mut self) {
        //// Equivalent formulation to substractinga stellar inertial position/velocity/acceleration:
        ////
        //// heliocentric_position = inertial_position + (1/star.mass) * sum(particle.mass * particle.position)
        //// heliocentric_velocity = inertial_velocity + (1/star.mass) * sum(particle.mass * particle.velocity)
        //// heliocentric_acceleration = inertial_acceleration + (1/star.mass) * sum(particle.mass * particle.acceleration)
        self.inertial_to_heliocentric_posvel();
        self.inertial_to_heliocentric_acc();
    }

    fn inertial_to_heliocentric_posvel(&mut self) {
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.position.x = particle.inertial_position.x - star.inertial_position.x;
                particle.position.y = particle.inertial_position.y - star.inertial_position.y;
                particle.position.z = particle.inertial_position.z - star.inertial_position.z;
                particle.velocity.x = particle.inertial_velocity.x - star.inertial_velocity.x;
                particle.velocity.y = particle.inertial_velocity.y - star.inertial_velocity.y;
                particle.velocity.z = particle.inertial_velocity.z - star.inertial_velocity.z;
            }
            star.position.x = 0.;
            star.position.y = 0.;
            star.position.z = 0.;
            star.velocity.x = 0.;
            star.velocity.y = 0.;
            star.velocity.z = 0.;
        }
    }

    fn inertial_to_heliocentric_acc(&mut self) {
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.acceleration.x = particle.inertial_acceleration.x - star.inertial_acceleration.x;
                particle.acceleration.y = particle.inertial_acceleration.y - star.inertial_acceleration.y;
                particle.acceleration.z = particle.inertial_acceleration.z - star.inertial_acceleration.z;
            }
            star.acceleration.x = 0.;
            star.acceleration.y = 0.;
            star.acceleration.z = 0.;
        }
    }

}
