use std;
use std::fs::{File};
use std::io::{Write};
use super::super::constants::*;
use super::super::particle::Particles;
use super::super::general::calculate_keplerian_orbital_elements;

pub trait Integrator {
    fn new(time_step: f64, time_limit: f64, particles: Particles) -> Self;
    fn iterate(&mut self, output_file: &mut File) -> Result<(), String>;
}


pub struct LeapFrog {
    time_step: f64,
    time_limit: f64,
    particles: Particles,
    current_time: f64,
    last_print_time: f64,
}

//pub struct Simulation : SimulationA {
    //current_time: f64,
    //last_print_time: f64,
//}

//pub struct Simulation {
    //time_step: f64,
    //time_limit: f64,
    //particles : [Particle; N_PARTICLES],
    //current_time: f64,
    //last_print_time: f64,
    ////// Integrator IAS15 data:
    //integrator_iterations_max_exceeded : i32,  // Count how many times the iteration did not converge
    //time_step_last_success: f64,			// Last accepted timestep (corresponding to br and er)
    //b: [[f64; 3*N_PARTICLES]; 7], // Coefficient b: acceleration dimension
    //br: [[f64; 3*N_PARTICLES]; 7], // Previous b
    //g: [[f64; 3*N_PARTICLES]; 7], // Coefficient g (it can also be expressed in terms of b)
    //e: [[f64; 3*N_PARTICLES]; 7],
    //er: [[f64; 3*N_PARTICLES]; 7], // Previous g
    //at: [f64; 3*N_PARTICLES], // Temporary buffer for acceleration
    //x0: [f64; 3*N_PARTICLES], // Temporary buffer for position (used for initial values at h=0)
    //v0: [f64; 3*N_PARTICLES], // Temporary buffer for velcoity (used for initial values at h=0)
    //a0: [f64; 3*N_PARTICLES], // Temporary buffer for acceleration (used for initial values at h=0)
    //// Compensated summation coefficients
    //csx : [f64; 3*N_PARTICLES],
    //csv : [f64; 3*N_PARTICLES],
    //s: [f64; 9], // Summation coefficients
//}

impl Integrator for LeapFrog {

    fn new(time_step: f64, time_limit: f64, particles: Particles) -> LeapFrog {
        LeapFrog {
                    time_step:time_step,
                    time_limit:time_limit,
                    last_print_time:0.,
                    particles:particles,
                    current_time:0.,
                    //integrator_iterations_max_exceeded:0,
                    //time_step_last_success:0.,
                    //b  :  [[0.; 3*N_PARTICLES]; 7],
                    //g  :  [[0.; 3*N_PARTICLES]; 7],
                    //e  :  [[0.; 3*N_PARTICLES]; 7],
                    //br :  [[0.; 3*N_PARTICLES]; 7],
                    //er :  [[0.; 3*N_PARTICLES]; 7],
                    //at  : [0.; 3*N_PARTICLES],
                    //x0  : [0.; 3*N_PARTICLES],
                    //v0  : [0.; 3*N_PARTICLES],
                    //a0  : [0.; 3*N_PARTICLES],
                    //csx : [0.; 3*N_PARTICLES],
                    //csv : [0.; 3*N_PARTICLES],
                    //s   : [0.; 9],
                    }
    }

    fn iterate(&mut self, output_file: &mut File) -> Result<(), String> {
        // A 'DKD'-like integrator will do the first 'D' part.
        self.integrator_part1();

        // Calculate accelerations.
        self.particles.gravity_calculate_acceleration();
        // Calculate non-gravity accelerations.
        self.particles.calculate_additional_forces();

        // A 'DKD'-like integrator will do the 'KD' part.
        self.integrator_part2();

        self.print_output(output_file);

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
        for i in 0..N_PARTICLES {
            self.particles.particles[i].position.x += 0.5 * self.time_step * self.particles.particles[i].velocity.x;
            self.particles.particles[i].position.y += 0.5 * self.time_step * self.particles.particles[i].velocity.y;
            self.particles.particles[i].position.z += 0.5 * self.time_step * self.particles.particles[i].velocity.z;

            self.particles.particles[i].spin.x += 0.5 * self.time_step * self.particles.particles[i].dspin_dt.x;
            self.particles.particles[i].spin.y += 0.5 * self.time_step * self.particles.particles[i].dspin_dt.y;
            self.particles.particles[i].spin.z += 0.5 * self.time_step * self.particles.particles[i].dspin_dt.z;
        }
        self.current_time += self.time_step/2.;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        for i in 0..N_PARTICLES{
            self.particles.particles[i].velocity.x += self.time_step * self.particles.particles[i].acceleration.x;
            self.particles.particles[i].velocity.y += self.time_step * self.particles.particles[i].acceleration.y;
            self.particles.particles[i].velocity.z += self.time_step * self.particles.particles[i].acceleration.z;
            self.particles.particles[i].position.x += 0.5 * self.time_step * self.particles.particles[i].velocity.x;
            self.particles.particles[i].position.y += 0.5 * self.time_step * self.particles.particles[i].velocity.y;
            self.particles.particles[i].position.z += 0.5 * self.time_step * self.particles.particles[i].velocity.z;

            //// SBC
            //if i == 1 {
                ////println!("--- {}", self.current_time);
                //println!("{:e} {:e} {:e}", self.particles.particles[i].spin.x, self.particles.particles[i].spin.y, self.particles.particles[i].spin.z);
            //}
            //self.particles.particles[i].spin.x += self.time_step * self.particles.particles[i].dspin_dt.x;
            //self.particles.particles[i].spin.y += self.time_step * self.particles.particles[i].dspin_dt.y;
            //self.particles.particles[i].spin.z += self.time_step * self.particles.particles[i].dspin_dt.z;
            self.particles.particles[i].spin.x += 0.5 * self.time_step * self.particles.particles[i].dspin_dt.x;
            self.particles.particles[i].spin.y += 0.5 * self.time_step * self.particles.particles[i].dspin_dt.y;
            self.particles.particles[i].spin.z += 0.5 * self.time_step * self.particles.particles[i].dspin_dt.z;

            //// SBC
            //if i == 1 {
                //println!(":: {:e} {:e} {:e}", self.particles.particles[i].spin.x, self.particles.particles[i].spin.y, self.particles.particles[i].spin.z);
            //}
        }
        self.current_time += self.time_step/2.;
    }


    fn print_output(&mut self, output_file: &mut File) {
        // Serialization: https://blog.safaribooksonline.com/2014/01/28/network-communication-serialization-rust/
        //let n_iterations: f64 = 2000.; // Print every N iterations
        //if (self.current_time == 0.) || (self.current_time >= self.time_limit)
            //||  ((self.current_time/n_iterations).floor() != ((self.current_time - self.time_step)/n_iterations).floor()) {
        // Print ever n_years
        let n_years: f64 = 100.;
        if (self.current_time == 0.) || (self.last_print_time == 0.) || (self.current_time >= self.time_limit)
            ||  (self.last_print_time + 365.25*n_years < self.current_time) {

            if self.last_print_time == 0. {
                // Header
                let _ = write!(output_file, "current_time\ttime_step\tparticle\t");
                let _ = write!(output_file, "position_x\tposition_y\tposition_z\t");
                let _ = write!(output_file, "velocity_x\tvelocity_y\tvelocity_z\t");
                let _ = write!(output_file, "acceleration_x\tacceleration_y\tacceleration_z\t");
                let _ = write!(output_file, "spin_x\tspin_y\tspin_z\t");
                let _ = write!(output_file, "dspin_dt_x\tdspin_dt_y\tdspin_dt_z\t");
                let _ = write!(output_file, "torque_x\ttorque_y\ttorque_z\t");
                let _ = write!(output_file, "orthogonal_component_of_the_tidal_force\t");
                let _ = write!(output_file, "radial_component_of_the_tidal_force\tradial_component_of_the_tidal_force_conservative_part\tradial_component_of_the_tidal_force_dissipative_part\t");
                let _ = write!(output_file, "tidal_acceleration_x\ttidal_acceleration_y\ttidal_acceleration_z\t");
                let _ = write!(output_file, "radial_velocity\tnorm_velocity_vector\tdistance\t");
                let _ = write!(output_file, "semi-major_axis\tperihelion_distance\teccentricity\tinclination\tlongitude_of_perihelion\tlongitude_of_ascending_node\tmean_anomaly\t");
                let _ = write!(output_file, "orbital_angular_momentum_x\torbital_angular_momentum_y\torbital_angular_momentum_z\torbital_angular_momentum\t");
                let _ = write!(output_file, "denergy_dt\t");
                let _ = write!(output_file, "mass\tradius\t");
                let _ = write!(output_file, "dissipation_factor\tradius_of_gyration_2\tlove_number\t");
                let _ = write!(output_file, "\n");
            }

            self.last_print_time = self.current_time;
            let current_time_years = self.current_time/365.25;

            print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
            let _ = std::io::stdout().flush();

            //// SBC
            //for i in 1..N_PARTICLES  {
            for i in 0..N_PARTICLES  {
                let _ = write!(output_file, "{}\t{}\t{}\t", current_time_years, self.time_step, i);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].position.x, self.particles.particles[i].position.y, self.particles.particles[i].position.z);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].velocity.x, self.particles.particles[i].velocity.y, self.particles.particles[i].velocity.z);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].acceleration.x, self.particles.particles[i].acceleration.y, self.particles.particles[i].acceleration.z);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].spin.x, self.particles.particles[i].spin.y, self.particles.particles[i].spin.z);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].dspin_dt.x, self.particles.particles[i].dspin_dt.y, self.particles.particles[i].dspin_dt.z);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].torque.x, self.particles.particles[i].torque.y, self.particles.particles[i].torque.z);
                let _ = write!(output_file, "{:e}\t", self.particles.particles[i].orthogonal_component_of_the_tidal_force);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].radial_component_of_the_tidal_force, self.particles.particles[i].radial_component_of_the_tidal_force_conservative_part, self.particles.particles[i].radial_component_of_the_tidal_force_dissipative_part);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].tidal_acceleration.x, self.particles.particles[i].tidal_acceleration.y, self.particles.particles[i].tidal_acceleration.z);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].radial_velocity, self.particles.particles[i].norm_velocity_vector, self.particles.particles[i].distance);
                if i > 0 {
                    //// Not for the central body (star)
                    let star = 0;
                    let (semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly) = calculate_keplerian_orbital_elements(self.particles.particles[star].mass_g+self.particles.particles[i].mass_g, self.particles.particles[i].position, self.particles.particles[i].velocity);
                    // Calculation of orbital angular momentum (without mass and in AU^2/day)
                    let horb_x = self.particles.particles[i].position.y * self.particles.particles[i].velocity.z - self.particles.particles[i].position.z * self.particles.particles[i].velocity.y;
                    let horb_y = self.particles.particles[i].position.z * self.particles.particles[i].velocity.x - self.particles.particles[i].position.x * self.particles.particles[i].velocity.z;
                    let horb_z = self.particles.particles[i].position.x * self.particles.particles[i].velocity.y - self.particles.particles[i].position.y * self.particles.particles[i].velocity.x;
                    let horbn = (horb_x.powf(2.) + horb_y.powf(2.) + horb_z.powf(2.)).sqrt();
                    let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t", semimajor_axis, perihelion_distance, eccentricity, inclination, longitude_of_perihelion, longitude_of_ascending_node, mean_anomaly);
                    let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t{:e}\t", horb_x/horbn, horb_y/horbn, horb_z/horbn, horbn);
                    let _ = write!(output_file, "{:e}\t", self.particles.particles[i].denergy_dt);
                } else {
                    let _ = write!(output_file, "\t\t\t\t\t\t\t");
                    let _ = write!(output_file, "\t\t\t\t");
                    let _ = write!(output_file, "\t");
                }
                let _ = write!(output_file, "{:e}\t{:e}\t", self.particles.particles[i].mass, self.particles.particles[i].radius);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].dissipation_factor, self.particles.particles[i].radius_of_gyration_2, self.particles.particles[i].love_number);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].dissipation_factor, self.particles.particles[i].radius_of_gyration_2, self.particles.particles[i].love_number);
                let _ = write!(output_file, "{:e}\t{:e}\t{:e}\t", self.particles.particles[i].dissipation_factor, self.particles.particles[i].radius_of_gyration_2, self.particles.particles[i].love_number);
                let _ = write!(output_file, "\n");
            }
        } 
    }

}
