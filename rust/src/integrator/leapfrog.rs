use std;
use std::io::{Write, BufWriter};
use super::Integrator;
use super::super::constants::PRINT_EVERY_N_DAYS;
use super::super::particle::Particles;
use super::output::{write_bin_snapshot};

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


pub struct LeapFrog {
    time_step: f64,
    half_time_step: f64,
    time_limit: f64,
    universe: Particles,
    current_time: f64,
    current_iteration: u32,
    last_print_time: f64,
}

impl Integrator for LeapFrog {

    fn new(time_step: f64, time_limit: f64, particles: Particles) -> LeapFrog {
        LeapFrog {
                    time_step:time_step,
                    half_time_step:0.5*time_step,
                    time_limit:time_limit,
                    last_print_time:-1.,
                    universe:particles,
                    current_time:0.,
                    current_iteration:0,
                    }
    }

    fn iterate<T: Write>(&mut self, output_bin: &mut BufWriter<T>) -> Result<(), String> {
        // Output
        let add_header = self.last_print_time < 0.;
        let time_triger = self.last_print_time + PRINT_EVERY_N_DAYS <= self.current_time;
        if add_header || time_triger {
            write_bin_snapshot(output_bin, &self.universe, self.current_time, self.time_step);
            let current_time_years = self.current_time/365.25;
            print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
            let _ = std::io::stdout().flush();

            if add_header || time_triger {
                self.last_print_time = self.current_time;
            } 
        }

        // Calculate non-gravity accelerations.
        let only_dspin_dt = true;
        self.universe.calculate_additional_forces(self.current_time, only_dspin_dt);

        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        self.universe.gravity_calculate_acceleration();

        // Calculate non-gravity accelerations.
        let only_dspin_dt = false;
        self.universe.calculate_additional_forces(self.current_time, only_dspin_dt);

        // A 'DKD'-like integrator will do the 'KD' part.
        self.integrator_part2();
        self.current_iteration += 1;

        // Return
        if self.current_time+self.time_step > self.time_limit {
            Err("reached maximum time limit.".to_string())
        } else {
            Ok(())
        }
    }

}

impl LeapFrog {
    // Leapfrog integrator (Drift-Kick-Drift)
    // for non-rotating frame.
    #[allow(dead_code)]
    fn integrator_part1(&mut self) {
        for particle in self.universe.particles.iter_mut() {
            particle.position.x += self.half_time_step * particle.velocity.x;
            particle.position.y += self.half_time_step * particle.velocity.y;
            particle.position.z += self.half_time_step * particle.velocity.z;

            particle.spin.x += self.half_time_step * particle.dspin_dt.x;
            particle.spin.y += self.half_time_step * particle.dspin_dt.y;
            particle.spin.z += self.half_time_step * particle.dspin_dt.z;
        }
        self.current_time += self.half_time_step;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        for particle in self.universe.particles.iter_mut() {
            particle.velocity.x += self.time_step * particle.acceleration.x;
            particle.velocity.y += self.time_step * particle.acceleration.y;
            particle.velocity.z += self.time_step * particle.acceleration.z;
            particle.position.x += self.half_time_step * particle.velocity.x;
            particle.position.y += self.half_time_step * particle.velocity.y;
            particle.position.z += self.half_time_step * particle.velocity.z;

            particle.spin.x += self.half_time_step * particle.dspin_dt.x;
            particle.spin.y += self.half_time_step * particle.dspin_dt.y;
            particle.spin.z += self.half_time_step * particle.dspin_dt.z;
        }
        self.current_time += self.half_time_step;
    }

}
