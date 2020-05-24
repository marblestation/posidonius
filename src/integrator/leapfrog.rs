use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use super::Integrator;
use super::super::particles::Universe;
use super::super::particles::IgnoreGravityTerms;
use super::super::effects::WindEffect;
use super::super::effects::EvolutionType;
use super::output::{write_recovery_snapshot, write_historic_snapshot};
use time;
use std::path::Path;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::any::Any;

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


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LeapFrog {
    time_step: f64,
    half_time_step: f64,
    pub universe: Universe,
    pub current_time: f64,
    current_iteration: u32,
    pub recovery_snapshot_period: f64,
    pub historic_snapshot_period: f64,
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
        let universe_integrator = LeapFrog {
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
        universe_integrator
    }

}

impl Integrator for LeapFrog {

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn get_n_historic_snapshots(&self) -> usize {
        self.n_historic_snapshots
    }
    
    fn get_n_particles(&self) -> usize {
        self.universe.n_particles
    }
    
    fn get_current_time(&self) -> f64 {
        self.current_time
    }

    fn set_time_limit(&mut self, time_limit: f64) {
        if time_limit > 0. && self.universe.time_limit != time_limit {
            if time_limit > self.universe.time_limit && self.universe.consider_effects.evolution {
                // Check if the new time is in the range of the evolutionary model
                for (i, evolver) in self.universe.particles_evolvers.iter().enumerate() {
                    let is_an_evolving_body = match evolver.evolution {
                        EvolutionType::NonEvolving => { false },
                        _ => { true }
                    };
                    if is_an_evolving_body && evolver.time[evolver.time.len()-1] < time_limit {
                        panic!("Your new time limit ({} days) is greater than the maximum allowed age of the evolving body #{} ({} days)", time_limit, i+1, evolver.time[evolver.time.len()-1]);
                    };
                }
            } else if time_limit < self.universe.time_limit {
                if time_limit < self.current_time {
                    panic!("Your new time limit ({} days) is smaller than the current time ({} days)", time_limit, self.current_time);
                }
            }
            println!("[INFO {} UTC] The time limit changed from {} to {} days", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), self.universe.time_limit, time_limit);
            self.universe.time_limit = time_limit;
        }
    }

    fn set_snapshot_periods(&mut self, historic_snapshot_period: f64, recovery_snapshot_period: f64) {
        if historic_snapshot_period > 0. && self.historic_snapshot_period != historic_snapshot_period {
            println!("[INFO {} UTC] The historic snapshot period changed from {} to {} days", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), self.historic_snapshot_period, historic_snapshot_period);
            self.historic_snapshot_period = historic_snapshot_period;
        } else {
            println!("[INFO {} UTC] A historic snapshot will be saved every {} days", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), self.historic_snapshot_period);
        }
        
        if recovery_snapshot_period > 0. && self.recovery_snapshot_period != recovery_snapshot_period {
            println!("[INFO {} UTC] The recovery snapshot period changed from {} to {} days", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), self.recovery_snapshot_period, recovery_snapshot_period);
            self.recovery_snapshot_period = recovery_snapshot_period;
        } else {
            println!("[INFO {} UTC] A recovery snapshot will be saved every {} days", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), self.recovery_snapshot_period);
        }
    }

    fn initialize_physical_values(&mut self) {
        if self.current_time != 0. {
            panic!("Physical values cannot be initialized on a resumed simulation");
        }
        self.universe.calculate_norm_spin(); // Needed for evolution
        self.universe.calculate_particles_evolving_quantities(self.current_time); // Make sure we start with the good initial values
        self.universe.calculate_roche_radiuses(); // Needed for collision detection
    }

    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            self.universe.inertial_to_heliocentric();
            if self.universe.consider_effects.tides {
                self.universe.calculate_denergy_dt();
            }
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.n_historic_snapshots as f64*self.historic_snapshot_period; // Instead of self.current_time to avoid small deviations
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            if ! silent_mode {
                print!("Year: {:0.0} ({:0.1e}) | Time step: {:0.3} days                    \r", current_time_years, current_time_years, self.time_step);
                let _ = std::io::stdout().flush();
            }
        }

        let ignore_gravity_terms = IgnoreGravityTerms::None;
        let ignored_gravity_terms = ignore_gravity_terms;

        self.universe.inertial_to_heliocentric();
        // Calculate non-gravity accelerations.
        let evolution = true;
        let dangular_momentum_dt_per_moment_of_inertia = true;
        let accelerations = false;
        self.universe.calculate_additional_effects(self.current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms);


        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        self.universe.gravity_calculate_acceleration(ignore_gravity_terms);
        self.universe.inertial_to_heliocentric();

        // Calculate non-gravity accelerations.
        let evolution = true;
        let dangular_momentum_dt_per_moment_of_inertia = true;
        let accelerations = true;
        self.universe.calculate_additional_effects(self.current_time, evolution, dangular_momentum_dt_per_moment_of_inertia, accelerations, ignored_gravity_terms);


        // A 'DKD'-like integrator will do the 'KD' part.
        self.integrator_part2();
        self.current_iteration += 1;

        // Return
        if self.current_time+self.time_step > self.universe.time_limit {
            Err("Simulation completed".to_string())
        } else {
            Ok(first_snapshot_trigger || recovery_snapshot_time_trigger)
        }
    }

    fn write_recovery_snapshot(&mut self, snapshot_path: &Path, universe_history_writer: &mut BufWriter<File>) {
        self.last_recovery_snapshot_time = self.current_time;
        universe_history_writer.flush().unwrap();
        // Compute hash for this universe at this moment of time
        let mut s = DefaultHasher::new();
        self.hash = 0;
        self.hash(&mut s);
        self.hash = s.finish();
        write_recovery_snapshot(&snapshot_path, &self);
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
                particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.x;
                particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.y;
                particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.z;
            } else {
                particle.spin.x = particle.spin.x + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.x;
                particle.spin.y = particle.spin.y + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.y;
                particle.spin.z = particle.spin.z + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.z;
            }
            if particle.wind.effect != WindEffect::Disabled && particle.wind.parameters.output.factor != 0. {
                // TODO: Verify wind factor
                particle.spin.x += self.half_time_step * particle.wind.parameters.output.factor * particle.spin.x;
                particle.spin.y += self.half_time_step * particle.wind.parameters.output.factor * particle.spin.y;
                particle.spin.z += self.half_time_step * particle.wind.parameters.output.factor * particle.spin.z;
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
                particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.x;
                particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.y;
                particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.z;
            } else {
                particle.spin.x = particle.spin.x + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.x;
                particle.spin.y = particle.spin.y + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.y;
                particle.spin.z = particle.spin.z + self.half_time_step * particle.dangular_momentum_dt_per_moment_of_inertia.z;
            }
            if particle.wind.effect != WindEffect::Disabled && particle.wind.parameters.output.factor != 0. {
                // TODO: Verify wind factor
                particle.spin.x += self.half_time_step * particle.wind.parameters.output.factor * particle.spin.x;
                particle.spin.y += self.half_time_step * particle.wind.parameters.output.factor * particle.spin.y;
                particle.spin.z += self.half_time_step * particle.wind.parameters.output.factor * particle.spin.z;
            }
        }
        self.current_time += self.half_time_step;
    }

}
