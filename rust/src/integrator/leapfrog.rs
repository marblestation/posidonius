extern crate rusqlite;
use std;
use std::io::{Write, BufWriter};
use super::Integrator;
use super::super::constants::{PRINT_EVERY_N_DAYS, PRINT_EVERY_N_ITERATIONS};
use super::super::particle::Particles;
use super::output::{write_txt_snapshot, write_bin_snapshot, write_db_snapshot};



pub struct LeapFrog {
    time_step: f64,
    time_limit: f64,
    particles: Particles,
    current_time: f64,
    current_iteration: u32,
    last_print_time: f64,
    last_print_iteration: u32,
}

impl Integrator for LeapFrog {

    fn new(time_step: f64, time_limit: f64, particles: Particles) -> LeapFrog {
        LeapFrog {
                    time_step:time_step,
                    time_limit:time_limit,
                    last_print_time:0.,
                    last_print_iteration:0,
                    particles:particles,
                    current_time:0.,
                    current_iteration:0,
                    }
    }

    fn iterate<T: Write>(&mut self, output_txt: &mut BufWriter<T>, output_bin: &mut BufWriter<T>, output_db: &rusqlite::Connection) -> Result<(), String> {
        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        self.particles.gravity_calculate_acceleration();
        // Calculate non-gravity accelerations.
        self.particles.calculate_additional_forces();

        // A 'DKD'-like integrator will do the 'KD' part.
        self.integrator_part2();
        self.current_iteration += 1;

        let add_header = self.last_print_iteration == 0 && self.last_print_time == 0.;
        let iteration_triger = self.last_print_iteration + PRINT_EVERY_N_ITERATIONS <= self.current_iteration;
        let time_triger = self.last_print_time + PRINT_EVERY_N_DAYS <= self.current_time;
        if add_header || iteration_triger || time_triger {
            write_txt_snapshot(output_txt, &self.particles, self.current_time, self.time_step, add_header);
            write_bin_snapshot(output_bin, &self.particles, self.current_time, self.time_step);
            write_db_snapshot(&output_db, &self.particles, self.current_time, self.time_step, add_header);
            let current_time_years = self.current_time/365.25;
            print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
            let _ = std::io::stdout().flush();

            if add_header || time_triger {
                self.last_print_time = self.current_time;
            } else if iteration_triger {
                self.last_print_iteration = self.current_iteration;
            } 
        }

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
        let half_time_step = 0.5 * self.time_step;
        for particle in self.particles.particles.iter_mut() {
            particle.position.x += half_time_step * particle.velocity.x;
            particle.position.y += half_time_step * particle.velocity.y;
            particle.position.z += half_time_step * particle.velocity.z;

            particle.spin.x += half_time_step * particle.dspin_dt.x;
            particle.spin.y += half_time_step * particle.dspin_dt.y;
            particle.spin.z += half_time_step * particle.dspin_dt.z;
        }
        self.current_time += half_time_step;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        let half_time_step = 0.5 * self.time_step;
        for particle in self.particles.particles.iter_mut() {
            particle.velocity.x += self.time_step * particle.acceleration.x;
            particle.velocity.y += self.time_step * particle.acceleration.y;
            particle.velocity.z += self.time_step * particle.acceleration.z;
            particle.position.x += half_time_step * particle.velocity.x;
            particle.position.y += half_time_step * particle.velocity.y;
            particle.position.z += half_time_step * particle.velocity.z;

            //// SBC
            //if i == 1 {
                ////println!("--- {}", self.current_time);
                //println!("{:e} {:e} {:e}", particle.spin.x, particle.spin.y, particle.spin.z);
            //}
            //particle.spin.x += self.time_step * particle.dspin_dt.x;
            //particle.spin.y += self.time_step * particle.dspin_dt.y;
            //particle.spin.z += self.time_step * particle.dspin_dt.z;
            particle.spin.x += half_time_step * particle.dspin_dt.x;
            particle.spin.y += half_time_step * particle.dspin_dt.y;
            particle.spin.z += half_time_step * particle.dspin_dt.z;

            //// SBC
            //if i == 1 {
                //println!(":: {:e} {:e} {:e}", particle.spin.x, particle.spin.y, particle.spin.z);
            //}
        }
        self.current_time += half_time_step;
    }

}
