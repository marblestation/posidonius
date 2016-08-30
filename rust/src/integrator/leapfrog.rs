use std::io::{Write, BufWriter};
use super::Integrator;
use super::super::constants::{PRINT_EVERY_N_YEARS, PRINT_EVERY_N_ITERATIONS};
use super::super::particle::Particles;
use super::output::print_output;



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

    fn iterate<T: Write>(&mut self, output_writer: &mut BufWriter<T>) -> Result<(), String> {
        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        self.particles.gravity_calculate_acceleration();
        // Calculate non-gravity accelerations.
        self.particles.calculate_additional_forces();

        // A 'DKD'-like integrator will do the 'KD' part.
        self.integrator_part2();
        self.current_iteration += 1;

        let add_header = self.last_print_iteration == 0;
        if add_header || self.last_print_iteration + PRINT_EVERY_N_ITERATIONS <= self.current_iteration {
            print_output(output_writer, &self.particles, self.current_time, self.time_step, add_header);
            self.last_print_iteration = self.current_iteration;
        } else if self.last_print_time + 365.25*PRINT_EVERY_N_YEARS <= self.current_time {
            print_output(output_writer, &self.particles, self.current_time, self.time_step, add_header);
            self.last_print_time = self.current_time;
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
        for particle in self.particles.particles.iter_mut() {
            particle.position.x += 0.5 * self.time_step * particle.velocity.x;
            particle.position.y += 0.5 * self.time_step * particle.velocity.y;
            particle.position.z += 0.5 * self.time_step * particle.velocity.z;

            particle.spin.x += 0.5 * self.time_step * particle.dspin_dt.x;
            particle.spin.y += 0.5 * self.time_step * particle.dspin_dt.y;
            particle.spin.z += 0.5 * self.time_step * particle.dspin_dt.z;
        }
        self.current_time += self.time_step/2.;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        for particle in self.particles.particles.iter_mut() {
            particle.velocity.x += self.time_step * particle.acceleration.x;
            particle.velocity.y += self.time_step * particle.acceleration.y;
            particle.velocity.z += self.time_step * particle.acceleration.z;
            particle.position.x += 0.5 * self.time_step * particle.velocity.x;
            particle.position.y += 0.5 * self.time_step * particle.velocity.y;
            particle.position.z += 0.5 * self.time_step * particle.velocity.z;

            //// SBC
            //if i == 1 {
                ////println!("--- {}", self.current_time);
                //println!("{:e} {:e} {:e}", particle.spin.x, particle.spin.y, particle.spin.z);
            //}
            //particle.spin.x += self.time_step * particle.dspin_dt.x;
            //particle.spin.y += self.time_step * particle.dspin_dt.y;
            //particle.spin.z += self.time_step * particle.dspin_dt.z;
            particle.spin.x += 0.5 * self.time_step * particle.dspin_dt.x;
            particle.spin.y += 0.5 * self.time_step * particle.dspin_dt.y;
            particle.spin.z += 0.5 * self.time_step * particle.dspin_dt.z;

            //// SBC
            //if i == 1 {
                //println!(":: {:e} {:e} {:e}", particle.spin.x, particle.spin.y, particle.spin.z);
            //}
        }
        self.current_time += self.time_step/2.;
    }

}
