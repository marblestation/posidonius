use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use super::Integrator;
use super::super::particles::Universe;
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
    universe: Universe,
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
        let mut s = DefaultHasher::new();
        universe_integrator.hash(&mut s);
        universe_integrator.hash = s.finish();
        universe_integrator
    }

    pub fn restore_snapshot(universe_integrator_snapshot_path: &Path) -> Result<LeapFrog, String> {
        let universe_integrator: LeapFrog;
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
            println!("INFO: Restored previous simulation from '{}'", universe_integrator_snapshot_path.display());
            return Ok(universe_integrator);
        } else {
            return Err(format!("File does not exist"));
        }
    }

}

impl Integrator for LeapFrog {


    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.current_time;
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
            let _ = std::io::stdout().flush();
        }

        // Calculate non-gravity accelerations.
        self.universe.calculate_position_velocity_and_spin_dependent_quantities();
        self.universe.calculate_particles_evolving_quantities(self.current_time);
        self.universe.calculate_torque_and_dspin_dt();


        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        let integrator_is_whfasthelio = false;
        self.universe.gravity_calculate_acceleration(integrator_is_whfasthelio);

        // Calculate non-gravity accelerations.
        self.universe.calculate_position_velocity_and_spin_dependent_quantities();
        self.universe.calculate_particles_evolving_quantities(self.current_time);
        self.universe.calculate_torque_and_dspin_dt();
        self.universe.calculate_additional_accelerations();


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
    }

}

impl LeapFrog {
    // Leapfrog integrator (Drift-Kick-Drift)
    // for non-rotating frame.
    #[allow(dead_code)]
    fn integrator_part1(&mut self) {
        for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            particle.position.x += self.half_time_step * particle.velocity.x;
            particle.position.y += self.half_time_step * particle.velocity.y;
            particle.position.z += self.half_time_step * particle.velocity.z;

            particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.half_time_step * particle.dspin_dt.x + self.half_time_step * particle.wind_factor * particle.spin.x;
            particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.half_time_step * particle.dspin_dt.y + self.half_time_step * particle.wind_factor * particle.spin.y;
            particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.half_time_step * particle.dspin_dt.z + self.half_time_step * particle.wind_factor * particle.spin.z;
        }
        self.current_time += self.half_time_step;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            particle.velocity.x += self.time_step * particle.acceleration.x;
            particle.velocity.y += self.time_step * particle.acceleration.y;
            particle.velocity.z += self.time_step * particle.acceleration.z;
            particle.position.x += self.half_time_step * particle.velocity.x;
            particle.position.y += self.half_time_step * particle.velocity.y;
            particle.position.z += self.half_time_step * particle.velocity.z;

            particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.half_time_step * particle.dspin_dt.x + self.half_time_step * particle.wind_factor * particle.spin.x;
            particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.half_time_step * particle.dspin_dt.y + self.half_time_step * particle.wind_factor * particle.spin.y;
            particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.half_time_step * particle.dspin_dt.z + self.half_time_step * particle.wind_factor * particle.spin.z;
        }
        self.current_time += self.half_time_step;
    }

}
