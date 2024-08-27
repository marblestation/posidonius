extern crate time;
use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use serde::{Serialize, Deserialize};
use serde_big_array::BigArray;
use super::Integrator;
use super::super::constants::{INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT, INTEGRATOR_EPSILON, INTEGRATOR_EPSILON_GLOBAL, INTEGRATOR_MIN_DT, INTEGRATOR_MAX_DT, SAFETY_FACTOR, MAX_PARTICLES};
use super::super::particles::Universe;
use super::super::particles::IgnoreGravityTerms;
use super::super::effects::GeneralRelativityImplementation;
use super::super::effects::EvolutionType;
use super::output::{write_recovery_snapshot, write_historic_snapshot};
use time::{OffsetDateTime, format_description};
use std::path::Path;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::any::Any;


///https://arxiv.org/abs/1409.4779
///IAS15: A fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine
///precision over a billion orbits
///
/// https://arxiv.org/pdf/1110.4876v2.pdf
/// variable time-steps also
/// break the symplectic nature of an integrator.

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Ias15 {
    time_step: f64,
    pub universe: Universe,
    pub current_time: f64,
    current_iteration: u32,
    pub recovery_snapshot_period: f64,
    pub historic_snapshot_period: f64,
    last_recovery_snapshot_time: f64,
    last_historic_snapshot_time: f64,
    pub n_historic_snapshots: usize,
    pub hash: u64,
    //// Integrator IAS15 data:
    n_particles: usize,
    integrator_iterations_max_exceeded : i32,  // Count how many times the iteration did not converge
    time_step_last_success: f64,			// Last accepted timestep (corresponding to br and er)
    //#[serde(with = "BigArray")]
    //b: [[f64; 3*MAX_PARTICLES]; 7], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_0: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_1: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_2: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_3: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_4: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_5: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    b_6: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    //#[serde(with = "BigArray")]
    //br: [[f64; 3*MAX_PARTICLES]; 7], // Previous b
    #[serde(with = "BigArray")]
    br_0: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    br_1: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    br_2: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    br_3: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    br_4: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    br_5: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    br_6: [f64; 3*MAX_PARTICLES], // Previous b
    //#[serde(with = "BigArray")]
    //g: [[f64; 3*MAX_PARTICLES]; 7], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_0: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_1: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_2: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_3: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_4: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_5: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    g_6: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    //#[serde(with = "BigArray")]
    //e: [[f64; 3*MAX_PARTICLES]; 7],
    #[serde(with = "BigArray")]
    e_0: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    e_1: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    e_2: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    e_3: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    e_4: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    e_5: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    e_6: [f64; 3*MAX_PARTICLES],
    //#[serde(with = "BigArray")]
    //er: [[f64; 3*MAX_PARTICLES]; 7], // Previous g
    #[serde(with = "BigArray")]
    er_0: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    er_1: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    er_2: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    er_3: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    er_4: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    er_5: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    er_6: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    at: [f64; 3*MAX_PARTICLES], // Temporary buffer for acceleration
    #[serde(with = "BigArray")]
    x0: [f64; 3*MAX_PARTICLES], // Temporary buffer for position (used for initial values at h=0)
    #[serde(with = "BigArray")]
    v0: [f64; 3*MAX_PARTICLES], // Temporary buffer for velocity (used for initial values at h=0)
    #[serde(with = "BigArray")]
    a0: [f64; 3*MAX_PARTICLES], // Temporary buffer for acceleration (used for initial values at h=0)
    // spin coeff
    //#[serde(with = "BigArray")]
    //sb: [[f64; 3*MAX_PARTICLES]; 7], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_0: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_1: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_2: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_3: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_4: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_5: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    #[serde(with = "BigArray")]
    sb_6: [f64; 3*MAX_PARTICLES], // Coefficient b: acceleration dimension
    //#[serde(with = "BigArray")]
    //sbr: [[f64; 3*MAX_PARTICLES]; 7], // Previous b
    #[serde(with = "BigArray")]
    sbr_0: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    sbr_1: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    sbr_2: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    sbr_3: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    sbr_4: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    sbr_5: [f64; 3*MAX_PARTICLES], // Previous b
    #[serde(with = "BigArray")]
    sbr_6: [f64; 3*MAX_PARTICLES], // Previous b
    //#[serde(with = "BigArray")]
    //sg: [[f64; 3*MAX_PARTICLES]; 7], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_0: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_1: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_2: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_3: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_4: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_5: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    #[serde(with = "BigArray")]
    sg_6: [f64; 3*MAX_PARTICLES], // Coefficient g (it can also be expressed in terms of b)
    //#[serde(with = "BigArray")]
    //se: [[f64; 3*MAX_PARTICLES]; 7],
    #[serde(with = "BigArray")]
    se_0: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    se_1: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    se_2: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    se_3: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    se_4: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    se_5: [f64; 3*MAX_PARTICLES],
    #[serde(with = "BigArray")]
    se_6: [f64; 3*MAX_PARTICLES],
    //#[serde(with = "BigArray")]
    //ser: [[f64; 3*MAX_PARTICLES]; 7], // Previous g
    #[serde(with = "BigArray")]
    ser_0: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    ser_1: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    ser_2: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    ser_3: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    ser_4: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    ser_5: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    ser_6: [f64; 3*MAX_PARTICLES], // Previous g
    #[serde(with = "BigArray")]
    dangular_momentum_dtt: [f64; 3*MAX_PARTICLES], // Temporary buffer for dangular_momentum_dt
    #[serde(with = "BigArray")]
    dangular_momentum_dt0: [f64; 3*MAX_PARTICLES], // Temporary buffer for dangular_momentum_dt (used for initial values at h=0)
    #[serde(with = "BigArray")]
    angular_momentum0: [f64; 3*MAX_PARTICLES], // Temporary buffer for angular_momentum (used for initial values at h=0)
    // Compensated summation coefficients
    #[serde(with = "BigArray")]
    csx : [f64; 3*MAX_PARTICLES], // position
    #[serde(with = "BigArray")]
    csv : [f64; 3*MAX_PARTICLES], // velocity
    #[serde(with = "BigArray")]
    css : [f64; 3*MAX_PARTICLES], // spin
    s: [f64; 9], // Summation coefficients
}

impl Hash for Ias15 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash only works with certain types (for instance, it does not work with f64 values)
        // thus we convert the whole integrator to a string thanks to the debug trait
        // and we hash that value
        format!("{:?}", self).hash(state); 
    }
}

impl Ias15 {
    pub fn new(time_step: f64, recovery_snapshot_period: f64, historic_snapshot_period: f64, universe: Universe) -> Ias15 {
        let n_particles = universe.n_particles;
        let universe_integrator = Ias15 {
                    time_step:time_step,
                    recovery_snapshot_period:recovery_snapshot_period,
                    historic_snapshot_period:historic_snapshot_period,
                    last_recovery_snapshot_time: -1.,
                    last_historic_snapshot_time: -1.,
                    n_historic_snapshots:0,
                    hash: 0,
                    universe:universe,
                    current_time:0.,
                    current_iteration:0,
                    n_particles:n_particles,
                    integrator_iterations_max_exceeded:0,
                    time_step_last_success:0.,
                    //b :   [[0.; 3*MAX_PARTICLES]; 7],
                    b_0 :   [0.; 3*MAX_PARTICLES],
                    b_1 :   [0.; 3*MAX_PARTICLES],
                    b_2 :   [0.; 3*MAX_PARTICLES],
                    b_3 :   [0.; 3*MAX_PARTICLES],
                    b_4 :   [0.; 3*MAX_PARTICLES],
                    b_5 :   [0.; 3*MAX_PARTICLES],
                    b_6 :   [0.; 3*MAX_PARTICLES],
                    //g  :  [[0.; 3*MAX_PARTICLES]; 7],
                    g_0  :  [0.; 3*MAX_PARTICLES],
                    g_1  :  [0.; 3*MAX_PARTICLES],
                    g_2  :  [0.; 3*MAX_PARTICLES],
                    g_3  :  [0.; 3*MAX_PARTICLES],
                    g_4  :  [0.; 3*MAX_PARTICLES],
                    g_5  :  [0.; 3*MAX_PARTICLES],
                    g_6  :  [0.; 3*MAX_PARTICLES],
                    //e  :  [[0.; 3*MAX_PARTICLES]; 7],
                    e_0  :  [0.; 3*MAX_PARTICLES],
                    e_1  :  [0.; 3*MAX_PARTICLES],
                    e_2  :  [0.; 3*MAX_PARTICLES],
                    e_3  :  [0.; 3*MAX_PARTICLES],
                    e_4  :  [0.; 3*MAX_PARTICLES],
                    e_5  :  [0.; 3*MAX_PARTICLES],
                    e_6  :  [0.; 3*MAX_PARTICLES],
                    //br :  [[0.; 3*MAX_PARTICLES]; 7],
                    br_0 :  [0.; 3*MAX_PARTICLES],
                    br_1 :  [0.; 3*MAX_PARTICLES],
                    br_2 :  [0.; 3*MAX_PARTICLES],
                    br_3 :  [0.; 3*MAX_PARTICLES],
                    br_4 :  [0.; 3*MAX_PARTICLES],
                    br_5 :  [0.; 3*MAX_PARTICLES],
                    br_6 :  [0.; 3*MAX_PARTICLES],
                    //er :  [[0.; 3*MAX_PARTICLES]; 7],
                    er_0 :  [0.; 3*MAX_PARTICLES],
                    er_1 :  [0.; 3*MAX_PARTICLES],
                    er_2 :  [0.; 3*MAX_PARTICLES],
                    er_3 :  [0.; 3*MAX_PARTICLES],
                    er_4 :  [0.; 3*MAX_PARTICLES],
                    er_5 :  [0.; 3*MAX_PARTICLES],
                    er_6 :  [0.; 3*MAX_PARTICLES],
                    at  : [0.; 3*MAX_PARTICLES],
                    x0  : [0.; 3*MAX_PARTICLES],
                    v0  : [0.; 3*MAX_PARTICLES],
                    a0  : [0.; 3*MAX_PARTICLES],
                    //sb :   [[0.; 3*MAX_PARTICLES]; 7],
                    sb_0 :   [0.; 3*MAX_PARTICLES],
                    sb_1 :   [0.; 3*MAX_PARTICLES],
                    sb_2 :   [0.; 3*MAX_PARTICLES],
                    sb_3 :   [0.; 3*MAX_PARTICLES],
                    sb_4 :   [0.; 3*MAX_PARTICLES],
                    sb_5 :   [0.; 3*MAX_PARTICLES],
                    sb_6 :   [0.; 3*MAX_PARTICLES],
                    //sg  :  [[0.; 3*MAX_PARTICLES]; 7],
                    sg_0  :  [0.; 3*MAX_PARTICLES],
                    sg_1  :  [0.; 3*MAX_PARTICLES],
                    sg_2  :  [0.; 3*MAX_PARTICLES],
                    sg_3  :  [0.; 3*MAX_PARTICLES],
                    sg_4  :  [0.; 3*MAX_PARTICLES],
                    sg_5  :  [0.; 3*MAX_PARTICLES],
                    sg_6  :  [0.; 3*MAX_PARTICLES],
                    //se  :  [[0.; 3*MAX_PARTICLES]; 7],
                    se_0  :  [0.; 3*MAX_PARTICLES],
                    se_1  :  [0.; 3*MAX_PARTICLES],
                    se_2  :  [0.; 3*MAX_PARTICLES],
                    se_3  :  [0.; 3*MAX_PARTICLES],
                    se_4  :  [0.; 3*MAX_PARTICLES],
                    se_5  :  [0.; 3*MAX_PARTICLES],
                    se_6  :  [0.; 3*MAX_PARTICLES],
                    //sbr :  [[0.; 3*MAX_PARTICLES]; 7],
                    sbr_0 :  [0.; 3*MAX_PARTICLES],
                    sbr_1 :  [0.; 3*MAX_PARTICLES],
                    sbr_2 :  [0.; 3*MAX_PARTICLES],
                    sbr_3 :  [0.; 3*MAX_PARTICLES],
                    sbr_4 :  [0.; 3*MAX_PARTICLES],
                    sbr_5 :  [0.; 3*MAX_PARTICLES],
                    sbr_6 :  [0.; 3*MAX_PARTICLES],
                    //ser :  [[0.; 3*MAX_PARTICLES]; 7],
                    ser_0 :  [0.; 3*MAX_PARTICLES],
                    ser_1 :  [0.; 3*MAX_PARTICLES],
                    ser_2 :  [0.; 3*MAX_PARTICLES],
                    ser_3 :  [0.; 3*MAX_PARTICLES],
                    ser_4 :  [0.; 3*MAX_PARTICLES],
                    ser_5 :  [0.; 3*MAX_PARTICLES],
                    ser_6 :  [0.; 3*MAX_PARTICLES],
                    dangular_momentum_dtt  : [0.; 3*MAX_PARTICLES],
                    dangular_momentum_dt0  : [0.; 3*MAX_PARTICLES],
                    angular_momentum0  : [0.; 3*MAX_PARTICLES],
                    csx : [0.; 3*MAX_PARTICLES],
                    csv : [0.; 3*MAX_PARTICLES],
                    css : [0.; 3*MAX_PARTICLES],
                    s   : [0.; 9],
                    };
        universe_integrator
    }

}

impl Integrator for Ias15 {

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
            println!("[INFO {} UTC] The time limit changed from {} to {} days", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), self.universe.time_limit, time_limit);
            self.universe.time_limit = time_limit;
        }
    }

    fn set_snapshot_periods(&mut self, historic_snapshot_period: f64, recovery_snapshot_period: f64) {
        if historic_snapshot_period > 0. && self.historic_snapshot_period != historic_snapshot_period {
            println!("[INFO {} UTC] The historic snapshot period changed from {} to {} days", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), self.historic_snapshot_period, historic_snapshot_period);
            self.historic_snapshot_period = historic_snapshot_period;
        } else {
            println!("[INFO {} UTC] A historic snapshot will be saved every {} days", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), self.historic_snapshot_period);
        }
        
        if recovery_snapshot_period > 0. && self.recovery_snapshot_period != recovery_snapshot_period {
            println!("[INFO {} UTC] The recovery snapshot period changed from {} to {} days", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), self.recovery_snapshot_period, recovery_snapshot_period);
            self.recovery_snapshot_period = recovery_snapshot_period;
        } else {
            println!("[INFO {} UTC] A recovery snapshot will be saved every {} days", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), self.recovery_snapshot_period);
        }
    }

    fn initialize_physical_values(&mut self) {
        if self.current_time != 0. {
            panic!("Physical values cannot be initialized on a resumed simulation");
        }
        let evolution = true;
        self.universe.calculate_spin_and_evolving_quantities(self.current_time, evolution); // Make sure we start with the good initial values
        self.universe.calculate_roche_radiuses(); // Needed for collision detection
    }


    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            self.universe.inertial_to_heliocentric();
            let evolution = true;
            self.universe.calculate_spin_and_evolving_quantities(self.current_time, evolution);
            if self.universe.consider_effects.tides {
                self.universe.calculate_denergy_dt();
            }
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.current_time;
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            if ! silent_mode {
                print!("Year: {:0.0} ({:0.1e}) | Time step: {:0.3} days                    \r", current_time_years, current_time_years, self.time_step);
                let _ = std::io::stdout().flush();
            }
        }

        // Calculate accelerations.
        let ignore_gravity_terms = IgnoreGravityTerms::None;
        let ignored_gravity_terms = ignore_gravity_terms;
        self.universe.gravity_calculate_acceleration(ignore_gravity_terms);
        self.universe.inertial_to_heliocentric();
        // Calculate non-gravity accelerations.
        let evolution = true;
        let dangular_momentum_dt = true;
        let accelerations = true;
        self.universe.calculate_additional_effects(self.current_time, evolution, dangular_momentum_dt, accelerations, ignored_gravity_terms);
        self.universe.apply_acceleration_corrections();

        self.integrator();
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

impl Ias15 {

    #[allow(dead_code)]
    fn integrator(&mut self) {
    
        // Gauss-Radau spacings for substeps within a sequence, for the 15th order 
        // integrator. The sum of the h values should be 3.733333333333333
        let h : [f64; 8]	= [ 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626 ]; 


        ////// Constants to compute g values
        //
        //let mut n = 0;
        //let r : [f64, ..28];
        //for j in 1..8 {
            //for k in 0..j {
                //r[n] = 1. / (h[j] - h[k]);
                //n += 1;
            //}
        //}
        //
        let r : [f64; 28] = [0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147];


        ////// Constants to convert between b and g arrays (c = c21, c31, c32, c41, c42...)
        //
        //let mut c : [f64, ..21];
        //let mut d : [f64, ..21];
        //c[0] = - h[1];
        //d[0] =   h[1];
        //let mut n = 0;
        //for j in 2..7 {
            //n += 1;
            //c[n] = -h[j] * c[n-j+1]
            //d[n] =  h[1] * d[n-j+1]
            //for k in 2..j {
                //n += 1;
                //c[n] = c[n-j]  -  h[j] * c[n-j+1];
                //d[n] = d[n-j]  +  h[k] * d[n-j+1];
            //}
            //n += 1;
            //c[n] = c[n-j] - h[j];
            //d[n] = d[n-j] + h[j];
        //}
        //
        let c : [f64; 21] = [-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588];

        let d : [f64; 21] = [0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588];

        let mut consider_dangular_momentum_dt_from_general_relativity = false;
        if self.universe.consider_effects.general_relativity && self.universe.general_relativity_implementation == GeneralRelativityImplementation::Kidder1995 {
            consider_dangular_momentum_dt_from_general_relativity = true;
        }
        let integrate_spin = self.universe.consider_effects.tides || self.universe.consider_effects.rotational_flattening
                            || self.universe.consider_effects.evolution || consider_dangular_momentum_dt_from_general_relativity;

        loop {

            for (k, particle) in self.universe.particles[..self.universe.n_particles].iter().enumerate() {
                self.x0[3*k]   = particle.inertial_position.x;
                self.x0[3*k+1] = particle.inertial_position.y;
                self.x0[3*k+2] = particle.inertial_position.z;
                self.v0[3*k]   = particle.inertial_velocity.x;
                self.v0[3*k+1] = particle.inertial_velocity.y;
                self.v0[3*k+2] = particle.inertial_velocity.z;
                self.a0[3*k]   = particle.inertial_acceleration.x;
                self.a0[3*k+1] = particle.inertial_acceleration.y;  
                self.a0[3*k+2] = particle.inertial_acceleration.z;
                if integrate_spin {
                    self.dangular_momentum_dt0[3*k]   = particle.dangular_momentum_dt.x;
                    self.dangular_momentum_dt0[3*k+1] = particle.dangular_momentum_dt.y;
                    self.dangular_momentum_dt0[3*k+2] = particle.dangular_momentum_dt.z;
                    self.angular_momentum0[3*k] = particle.angular_momentum.x;
                    self.angular_momentum0[3*k+1] = particle.angular_momentum.y;
                    self.angular_momentum0[3*k+2] = particle.angular_momentum.z;
                }
            }

            let mut csb = [[0.; 3*MAX_PARTICLES]; 7];
            let mut cssb = [[0.; 3*MAX_PARTICLES]; 7];

            // Find g values from b values predicted at the last call (Eqs. 7 of Everhart)
            for k in 0..3*self.n_particles {
                self.g_0[k] = self.b_6[k]*d[15] + self.b_5[k]*d[10] + self.b_4[k]*d[6] + self.b_3[k]*d[3]  + self.b_2[k]*d[1]  + self.b_1[k]*d[0]  + self.b_0[k];
                self.g_1[k] = self.b_6[k]*d[16] + self.b_5[k]*d[11] + self.b_4[k]*d[7] + self.b_3[k]*d[4]  + self.b_2[k]*d[2]  + self.b_1[k];
                self.g_2[k] = self.b_6[k]*d[17] + self.b_5[k]*d[12] + self.b_4[k]*d[8] + self.b_3[k]*d[5]  + self.b_2[k];
                self.g_3[k] = self.b_6[k]*d[18] + self.b_5[k]*d[13] + self.b_4[k]*d[9] + self.b_3[k];
                self.g_4[k] = self.b_6[k]*d[19] + self.b_5[k]*d[14] + self.b_4[k];
                self.g_5[k] = self.b_6[k]*d[20] + self.b_5[k];
                self.g_6[k] = self.b_6[k];

                if integrate_spin {
                    self.sg_0[k] = self.sb_6[k]*d[15] + self.sb_5[k]*d[10] + self.sb_4[k]*d[6] + self.sb_3[k]*d[3]  + self.sb_2[k]*d[1]  + self.sb_1[k]*d[0]  + self.sb_0[k];
                    self.sg_1[k] = self.sb_6[k]*d[16] + self.sb_5[k]*d[11] + self.sb_4[k]*d[7] + self.sb_3[k]*d[4]  + self.sb_2[k]*d[2]  + self.sb_1[k];
                    self.sg_2[k] = self.sb_6[k]*d[17] + self.sb_5[k]*d[12] + self.sb_4[k]*d[8] + self.sb_3[k]*d[5]  + self.sb_2[k];
                    self.sg_3[k] = self.sb_6[k]*d[18] + self.sb_5[k]*d[13] + self.sb_4[k]*d[9] + self.sb_3[k];
                    self.sg_4[k] = self.sb_6[k]*d[19] + self.sb_5[k]*d[14] + self.sb_4[k];
                    self.sg_5[k] = self.sb_6[k]*d[20] + self.sb_5[k];
                    self.sg_6[k] = self.sb_6[k];
                }
            }

            let t_beginning = self.current_time;
            //let mut predictor_corrector_error: f64 = 1e10;
            let mut predictor_corrector_error: f64 = 1e300;
            let mut predictor_corrector_error_last: f64 = 2.;
            let mut iterations : i32 = 0;	

            // Predictor corrector loop
            // Stops if one of the following conditions is satisfied: 
            //   1) predictor_corrector_error better than 1e-16 
            //   2) predictor_corrector_error starts to oscillate
            //   3) more than 12 iterations
            //
            //  NOTE:
            //  In the original Everhart algorithm, there were six iterations 
            //  for first call to the subroutine and two for the rest of calls
            loop {
                if predictor_corrector_error < 1e-16 {
                    break;								// Quit predictor corrector loop
                }
                if iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error {
                    break;								// Quit predictor corrector loop
                }
                if iterations>=12 {
                    self.integrator_iterations_max_exceeded += 1;
                    const INTEGRATOR_ITERATIONS_WARNING: i32 = 10;
                    if self.integrator_iterations_max_exceeded == INTEGRATOR_ITERATIONS_WARNING {
                        println!("[WARNING {} UTC] At least {} predictor corrector loops in integrator IAS15 did not converge. This is typically an indication of the timestep being too large.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), INTEGRATOR_ITERATIONS_WARNING);
                    }
                    break;								// Quit predictor corrector loop
                }
                predictor_corrector_error_last = predictor_corrector_error;
                predictor_corrector_error = 0.;
                iterations += 1;

                for n in 1..8 {							// Loop over interval using Gauss-Radau spacings

                    // Calculate position predictors using Eqn. 9 of Everhart
                    self.s[0] = self.time_step * h[n];
                    self.s[1] = self.s[0] * self.s[0] / 2.;
                    self.s[2] = self.s[1] * h[n] / 3.;
                    self.s[3] = self.s[2] * h[n] / 2.;
                    self.s[4] = 3. * self.s[3] * h[n] / 5.;
                    self.s[5] = 2. * self.s[4] * h[n] / 3.;
                    self.s[6] = 5. * self.s[5] * h[n] / 7.;
                    self.s[7] = 3. * self.s[6] * h[n] / 4.;
                    self.s[8] = 7. * self.s[7] * h[n] / 9.;
                    
                    self.current_time = t_beginning + self.s[0];

                    //// Prepare particles arrays for force calculation
                    // Predict positions at interval n using self.b values
                    for (i, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                        let k0 : usize = 3*i+0;
                        let k1 : usize = 3*i+1;
                        let k2 : usize = 3*i+2;

                        // Equation 7 in paper 2015MNRAS.446.1424R
                        let xk0: f64  = -self.csx[k0] + (self.s[8]*self.b_6[k0] + self.s[7]*self.b_5[k0] + self.s[6]*self.b_4[k0] + self.s[5]*self.b_3[k0] + self.s[4]*self.b_2[k0] + self.s[3]*self.b_1[k0] + self.s[2]*self.b_0[k0] + self.s[1]*self.a0[k0] + self.s[0]*self.v0[k0] );
                        particle.inertial_position.x = xk0 + self.x0[k0];

                        let xk1: f64  = -self.csx[k1] + (self.s[8]*self.b_6[k1] + self.s[7]*self.b_5[k1] + self.s[6]*self.b_4[k1] + self.s[5]*self.b_3[k1] + self.s[4]*self.b_2[k1] + self.s[3]*self.b_1[k1] + self.s[2]*self.b_0[k1] + self.s[1]*self.a0[k1] + self.s[0]*self.v0[k1] );
                        particle.inertial_position.y = xk1 + self.x0[k1];

                        let xk2: f64 = -self.csx[k2] + (self.s[8]*self.b_6[k2] + self.s[7]*self.b_5[k2] + self.s[6]*self.b_4[k2] + self.s[5]*self.b_3[k2] + self.s[4]*self.b_2[k2] + self.s[3]*self.b_1[k2] + self.s[2]*self.b_0[k2] + self.s[1]*self.a0[k2] + self.s[0]*self.v0[k2] );
                        particle.inertial_position.z = xk2 + self.x0[k2];
                    }
                
                    if INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT {
                        // If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
                        self.s[0] = self.time_step * h[n];
                        self.s[1] =      self.s[0] * h[n] / 2.;
                        self.s[2] = 2. * self.s[1] * h[n] / 3.;
                        self.s[3] = 3. * self.s[2] * h[n] / 4.;
                        self.s[4] = 4. * self.s[3] * h[n] / 5.;
                        self.s[5] = 5. * self.s[4] * h[n] / 6.;
                        self.s[6] = 6. * self.s[5] * h[n] / 7.;
                        self.s[7] = 7. * self.s[6] * h[n] / 8.;

                        // Predict velocities at interval n using self.b values
                        for (i, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                            let k0 = 3*i+0;
                            let k1 = 3*i+1;
                            let k2 = 3*i+2;

                            // Equation 6 in paper 2015MNRAS.446.1424R
                            let vk0 =  -self.csv[k0] + self.s[7]*self.b_6[k0] + self.s[6]*self.b_5[k0] + self.s[5]*self.b_4[k0] + self.s[4]*self.b_3[k0] + self.s[3]*self.b_2[k0] + self.s[2]*self.b_1[k0] + self.s[1]*self.b_0[k0] + self.s[0]*self.a0[k0];
                            particle.inertial_velocity.x = vk0 + self.v0[k0];
                            let vk1 =  -self.csv[k1] + self.s[7]*self.b_6[k1] + self.s[6]*self.b_5[k1] + self.s[5]*self.b_4[k1] + self.s[4]*self.b_3[k1] + self.s[3]*self.b_2[k1] + self.s[2]*self.b_1[k1] + self.s[1]*self.b_0[k1] + self.s[0]*self.a0[k1];
                            particle.inertial_velocity.y = vk1 + self.v0[k1];
                            let vk2 =  -self.csv[k2] + self.s[7]*self.b_6[k2] + self.s[6]*self.b_5[k2] + self.s[5]*self.b_4[k2] + self.s[4]*self.b_3[k2] + self.s[3]*self.b_2[k2] + self.s[2]*self.b_1[k2] + self.s[1]*self.b_0[k2] + self.s[0]*self.a0[k2];
                            particle.inertial_velocity.z = vk2 + self.v0[k2];
                        }
                        
                        if integrate_spin {
                            // Spin integration
                            // Predict velocities at interval n using self.b values
                            for (i, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                                let k0 = 3*i+0;
                                let k1 = 3*i+1;
                                let k2 = 3*i+2;

                                // Equation 6 in paper 2015MNRAS.446.1424R
                                let angular_momentum_k0 =  -self.css[k0] + self.s[7]*self.sb_6[k0] + self.s[6]*self.sb_5[k0] + self.s[5]*self.sb_4[k0] + self.s[4]*self.sb_3[k0] + self.s[3]*self.sb_2[k0] + self.s[2]*self.sb_1[k0] + self.s[1]*self.sb_0[k0] + self.s[0]*self.dangular_momentum_dt0[k0];
                                particle.angular_momentum.x = angular_momentum_k0 + self.angular_momentum0[k0];
                                let angular_momentum_k1 =  -self.css[k1] + self.s[7]*self.sb_6[k1] + self.s[6]*self.sb_5[k1] + self.s[5]*self.sb_4[k1] + self.s[4]*self.sb_3[k1] + self.s[3]*self.sb_2[k1] + self.s[2]*self.sb_1[k1] + self.s[1]*self.sb_0[k1] + self.s[0]*self.dangular_momentum_dt0[k1];
                                particle.angular_momentum.y = angular_momentum_k1 + self.angular_momentum0[k1];
                                let angular_momentum_k2 =  -self.css[k2] + self.s[7]*self.sb_6[k2] + self.s[6]*self.sb_5[k2] + self.s[5]*self.sb_4[k2] + self.s[4]*self.sb_3[k2] + self.s[3]*self.sb_2[k2] + self.s[2]*self.sb_1[k2] + self.s[1]*self.sb_0[k2] + self.s[0]*self.dangular_momentum_dt0[k2];
                                particle.angular_momentum.z = angular_momentum_k2 + self.angular_momentum0[k2];
                            }
                        }
                    }

                    // Calculate accelerations.
                    let ignore_gravity_terms = IgnoreGravityTerms::None;
                    let ignored_gravity_terms = ignore_gravity_terms;
                    self.universe.gravity_calculate_acceleration(ignore_gravity_terms);
                    self.universe.inertial_to_heliocentric();
                    // Calculate non-gravity accelerations.
                    let evolution = true;
                    let dangular_momentum_dt = true;
                    let accelerations = true;
                    self.universe.calculate_additional_effects(self.current_time, evolution, dangular_momentum_dt, accelerations, ignored_gravity_terms);
                    self.universe.apply_acceleration_corrections();

                    for (k, particle) in self.universe.particles[..self.universe.n_particles].iter().enumerate() {
                        self.at[3*k]   = particle.inertial_acceleration.x;
                        self.at[3*k+1] = particle.inertial_acceleration.y;  
                        self.at[3*k+2] = particle.inertial_acceleration.z;

                        if integrate_spin {
                            // angular_momentum_dt = dmoment_of_inertia_spin_dt
                            self.dangular_momentum_dtt[3*k]   = particle.dangular_momentum_dt.x;
                            self.dangular_momentum_dtt[3*k+1] = particle.dangular_momentum_dt.y;
                            self.dangular_momentum_dtt[3*k+2] = particle.dangular_momentum_dt.z;
                        }
                    }

                    // Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
                    let always_zero = 0.;
                    match n {
                        1 => {
                                for k in 0..3*self.n_particles {
                                    let tmp = self.g_0[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_0[k]  = gk / r[0];
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], self.g_0[k] - tmp);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;

                                    if integrate_spin {
                                        let stmp = self.sg_0[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_0[k]  = sgk / r[0];
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], self.sg_0[k] - stmp);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                    }
                                }
                            },
                        2 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g_1[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_1[k]  = (gk/r[1] - self.g_0[k])/r[2];
                                    tmp = self.g_1[k] - tmp;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], tmp*c[0]);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_1[k], csb[1][k], tmp);
                                    self.b_1[k] = tmp_b;
                                    csb[1][k] = tmp_csb;

                                    if integrate_spin {
                                        let mut stmp = self.sg_1[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_1[k]  = (sgk/r[1] - self.sg_0[k])/r[2];
                                        stmp = self.sg_1[k] - stmp;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], stmp*c[0]);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_1[k], cssb[1][k], stmp);
                                        self.sb_1[k] = tmp_sb;
                                        cssb[1][k] = tmp_cssb;
                                    }
                                }
                            },
                        3 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g_2[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_2[k] = ((gk/r[3] - self.g_0[k])/r[4] - self.g_1[k])/r[5];
                                    tmp = self.g_2[k] - tmp;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], tmp*c[1]);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_1[k], csb[1][k], tmp*c[2]);
                                    self.b_1[k] = tmp_b;
                                    csb[1][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_2[k], csb[2][k], tmp);
                                    self.b_2[k] = tmp_b;
                                    csb[2][k] = tmp_csb;

                                    if integrate_spin {
                                        let mut stmp = self.sg_2[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_2[k] = ((sgk/r[3] - self.sg_0[k])/r[4] - self.sg_1[k])/r[5];
                                        stmp = self.sg_2[k] - stmp;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], stmp*c[1]);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_1[k], cssb[1][k], stmp*c[2]);
                                        self.sb_1[k] = tmp_sb;
                                        cssb[1][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_2[k], cssb[2][k], stmp);
                                        self.sb_2[k] = tmp_sb;
                                        cssb[2][k] = tmp_cssb;
                                    }
                                }
                            },
                        4 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g_3[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_3[k] = (((gk/r[6] - self.g_0[k])/r[7] - self.g_1[k])/r[8] - self.g_2[k])/r[9];
                                    tmp = self.g_3[k] - tmp;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], tmp*c[3]);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_1[k], csb[1][k], tmp*c[4]);
                                    self.b_1[k] = tmp_b;
                                    csb[1][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_2[k], csb[2][k], tmp*c[5]);
                                    self.b_2[k] = tmp_b;
                                    csb[2][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_3[k], csb[3][k], tmp);
                                    self.b_3[k] = tmp_b;
                                    csb[3][k] = tmp_csb;

                                    if integrate_spin {
                                        let mut stmp = self.sg_3[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_3[k] = (((sgk/r[6] - self.sg_0[k])/r[7] - self.sg_1[k])/r[8] - self.sg_2[k])/r[9];
                                        stmp = self.sg_3[k] - stmp;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], stmp*c[3]);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_1[k], cssb[1][k], stmp*c[4]);
                                        self.sb_1[k] = tmp_sb;
                                        cssb[1][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_2[k], cssb[2][k], stmp*c[5]);
                                        self.sb_2[k] = tmp_sb;
                                        cssb[2][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_3[k], cssb[3][k], stmp);
                                        self.sb_3[k] = tmp_sb;
                                        cssb[3][k] = tmp_cssb;
                                    }
                                }
                            },
                        5 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g_4[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_4[k] = ((((gk/r[10] - self.g_0[k])/r[11] - self.g_1[k])/r[12] - self.g_2[k])/r[13] - self.g_3[k])/r[14];
                                    tmp = self.g_4[k] - tmp;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], tmp*c[6]);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_1[k], csb[1][k], tmp*c[7]);
                                    self.b_1[k] = tmp_b;
                                    csb[1][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_2[k], csb[2][k], tmp*c[8]);
                                    self.b_2[k] = tmp_b;
                                    csb[2][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_3[k], csb[3][k], tmp*c[9]);
                                    self.b_3[k] = tmp_b;
                                    csb[3][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_4[k], csb[4][k], tmp);
                                    self.b_4[k] = tmp_b;
                                    csb[4][k] = tmp_csb;

                                    if integrate_spin {
                                        let mut stmp = self.sg_4[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_4[k] = ((((sgk/r[10] - self.sg_0[k])/r[11] - self.sg_1[k])/r[12] - self.sg_2[k])/r[13] - self.sg_3[k])/r[14];
                                        stmp = self.sg_4[k] - stmp;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], stmp*c[6]);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_1[k], cssb[1][k], stmp*c[7]);
                                        self.sb_1[k] = tmp_sb;
                                        cssb[1][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_2[k], cssb[2][k], stmp*c[8]);
                                        self.sb_2[k] = tmp_sb;
                                        cssb[2][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_3[k], cssb[3][k], stmp*c[9]);
                                        self.sb_3[k] = tmp_sb;
                                        cssb[3][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_4[k], cssb[4][k], stmp);
                                        self.sb_4[k] = tmp_sb;
                                        cssb[4][k] = tmp_cssb;
                                    }
                                }
                            },
                        6 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g_5[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_5[k] = (((((gk/r[15] - self.g_0[k])/r[16] - self.g_1[k])/r[17] - self.g_2[k])/r[18] - self.g_3[k])/r[19] - self.g_4[k])/r[20];
                                    tmp = self.g_5[k] - tmp;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], tmp*c[10]);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_1[k], csb[1][k], tmp*c[11]);
                                    self.b_1[k] = tmp_b;
                                    csb[1][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_2[k], csb[2][k], tmp*c[12]);
                                    self.b_2[k] = tmp_b;
                                    csb[2][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_3[k], csb[3][k], tmp*c[13]);
                                    self.b_3[k] = tmp_b;
                                    csb[3][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_4[k], csb[4][k], tmp*c[14]);
                                    self.b_4[k] = tmp_b;
                                    csb[4][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_5[k], csb[5][k], tmp);
                                    self.b_5[k] = tmp_b;
                                    csb[5][k] = tmp_csb;

                                    if integrate_spin {
                                        let mut stmp = self.sg_5[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_5[k] = (((((sgk/r[15] - self.sg_0[k])/r[16] - self.sg_1[k])/r[17] - self.sg_2[k])/r[18] - self.sg_3[k])/r[19] - self.sg_4[k])/r[20];
                                        stmp = self.sg_5[k] - stmp;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], stmp*c[10]);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_1[k], cssb[1][k], stmp*c[11]);
                                        self.sb_1[k] = tmp_sb;
                                        cssb[1][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_2[k], cssb[2][k], stmp*c[12]);
                                        self.sb_2[k] = tmp_sb;
                                        cssb[2][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_3[k], cssb[3][k], stmp*c[13]);
                                        self.sb_3[k] = tmp_sb;
                                        cssb[3][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_4[k], cssb[4][k], stmp*c[14]);
                                        self.sb_4[k] = tmp_sb;
                                        cssb[4][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_5[k], cssb[5][k], stmp);
                                        self.sb_5[k] = tmp_sb;
                                        cssb[5][k] = tmp_cssb;
                                    }
                                }
                            },
                        7 => {
                                let mut maxak: f64 = 0.0;
                                let mut maxb6ktmp: f64 = 0.0;
                                let mut max_dangular_momentum_dtk: f64 = 0.0;
                                let mut maxb6kstmp: f64 = 0.0;
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g_6[k];
                                    let (gk, gk_cs) = self.add_cs(self.at[k], always_zero, -self.a0[k]);
                                    let (gk, _gk_cs) = self.add_cs(gk, gk_cs, always_zero);
                                    self.g_6[k] = ((((((gk/r[21] - self.g_0[k])/r[22] - self.g_1[k])/r[23] - self.g_2[k])/r[24] - self.g_3[k])/r[25] - self.g_4[k])/r[26] - self.g_5[k])/r[27];
                                    tmp = self.g_6[k] - tmp;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_0[k], csb[0][k], tmp*c[15]);
                                    self.b_0[k] = tmp_b;
                                    csb[0][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_1[k], csb[1][k], tmp*c[16]);
                                    self.b_1[k] = tmp_b;
                                    csb[1][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_2[k], csb[2][k], tmp*c[17]);
                                    self.b_2[k] = tmp_b;
                                    csb[2][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_3[k], csb[3][k], tmp*c[18]);
                                    self.b_3[k] = tmp_b;
                                    csb[3][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_4[k], csb[4][k], tmp*c[19]);
                                    self.b_4[k] = tmp_b;
                                    csb[4][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_5[k], csb[5][k], tmp*c[20]);
                                    self.b_5[k] = tmp_b;
                                    csb[5][k] = tmp_csb;
                                    let (tmp_b, tmp_csb) = self.add_cs(self.b_6[k], csb[6][k], tmp);
                                    self.b_6[k] = tmp_b;
                                    csb[6][k] = tmp_csb;

                                    let mut stmp = 0.;
                                    if integrate_spin {
                                        stmp = self.sg_6[k];
                                        let (sgk, sgk_cs) = self.add_cs(self.dangular_momentum_dtt[k], always_zero, -self.dangular_momentum_dt0[k]);
                                        let (sgk, _sgk_cs) = self.add_cs(sgk, sgk_cs, always_zero);
                                        self.sg_6[k] = ((((((sgk/r[21] - self.sg_0[k])/r[22] - self.sg_1[k])/r[23] - self.sg_2[k])/r[24] - self.sg_3[k])/r[25] - self.sg_4[k])/r[26] - self.sg_5[k])/r[27];
                                        stmp = self.sg_6[k] - stmp;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_0[k], cssb[0][k], stmp*c[15]);
                                        self.sb_0[k] = tmp_sb;
                                        cssb[0][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_1[k], cssb[1][k], stmp*c[16]);
                                        self.sb_1[k] = tmp_sb;
                                        cssb[1][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_2[k], cssb[2][k], stmp*c[17]);
                                        self.sb_2[k] = tmp_sb;
                                        cssb[2][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_3[k], cssb[3][k], stmp*c[18]);
                                        self.sb_3[k] = tmp_sb;
                                        cssb[3][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_4[k], cssb[4][k], stmp*c[19]);
                                        self.sb_4[k] = tmp_sb;
                                        cssb[4][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_5[k], cssb[5][k], stmp*c[20]);
                                        self.sb_5[k] = tmp_sb;
                                        cssb[5][k] = tmp_cssb;
                                        let (tmp_sb, tmp_cssb) = self.add_cs(self.sb_6[k], cssb[6][k], stmp);
                                        self.sb_6[k] = tmp_sb;
                                        cssb[6][k] = tmp_cssb;
                                    }
                                    
                                    // Monitor change in self.b_6[k] relative to self.at[k]. The predictor corrector scheme is converged if it is close to 0.
                                    if INTEGRATOR_EPSILON_GLOBAL {
                                        // global error estimate (equation 9)
                                        let ak  = self.at[k].abs();
                                        if ak.is_normal() && ak>maxak {
                                            maxak = ak;
                                        }
                                        let b6ktmp = tmp.abs();  // change of b6ktmp coefficient
                                        if b6ktmp.is_normal() && b6ktmp>maxb6ktmp {
                                            maxb6ktmp = b6ktmp;
                                        }
                                        if integrate_spin {
                                            let dangular_momentum_dttk  = self.dangular_momentum_dtt[k].abs();
                                            if dangular_momentum_dttk.is_normal() && dangular_momentum_dttk>max_dangular_momentum_dtk {
                                                max_dangular_momentum_dtk = dangular_momentum_dttk;
                                            }
                                            let b6kstmp = stmp.abs();  // change of b6ktmp coefficient
                                            if b6kstmp.is_normal() && b6kstmp>maxb6kstmp {
                                                maxb6kstmp = b6kstmp;
                                            }
                                        }
                                    } else {
                                        // local error estimate (equation 9)
                                        let mut predictor_corrector_error_a = 0.;
                                        let mut predictor_corrector_error_s = 0.;
                                        let ak  = self.at[k];
                                        let b6ktmp = tmp; 
                                        let errork = (b6ktmp/ak).abs();
                                        if errork.is_normal() && errork>predictor_corrector_error {
                                            predictor_corrector_error_a = errork;
                                        }
                                        //
                                        if integrate_spin {
                                            let dangular_momentum_dttk  = self.dangular_momentum_dtt[k];
                                            let b6kstmp = stmp; 
                                            let errorks = (b6kstmp/dangular_momentum_dttk).abs();
                                            if errorks.is_normal() && errorks>predictor_corrector_error {
                                                predictor_corrector_error_s = errorks;
                                            }
                                        }
                                        predictor_corrector_error = self.max(predictor_corrector_error_a, predictor_corrector_error_s);
                                    }
                                } 
                                if INTEGRATOR_EPSILON_GLOBAL {
                                    predictor_corrector_error = self.max(maxb6ktmp, maxb6kstmp)/self.max(maxak, max_dangular_momentum_dtk);
                                }
                            },
                        _ => { println!("[WARNING {} UTC] This should not happen because the loop stops at 7!", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap()); }

                    } // end match
                } // end loop over interval using Gauss-Radau spacings
            } // end predictor corrector loop


            // Set time back to initial value (will be updated below) 
            self.current_time = t_beginning;

            ////////////////////////////////////////////////////////////////////
            //// Find new timestep
            ////////////////////////////////////////////////////////////////////
            let dt_done = self.time_step;
            if INTEGRATOR_EPSILON > 0. {
                // Estimate error (given by last term in series expansion) 
                // There are two options:
                // INTEGRATOR_EPSILON_GLOBAL==true  (default)
                //   First, we determine the maximum acceleration and the maximum of the last term in the series. 
                //   Then, the two are divided.
                // INTEGRATOR_EPSILON_GLOBAL==false
                //   Here, the fractional error is calculated for each particle individually and we use the maximum of the fractional error.
                //   This might fail in cases where a particle does not experience any (physical) acceleration besides roundoff errors. 
                let mut integrator_error: f64 = 0.0;
                if INTEGRATOR_EPSILON_GLOBAL {
                    let mut maxak: f64 = 0.0;
                    let mut maxb6k: f64 = 0.0;
                    let mut max_dangular_momentum_dtk: f64 = 0.0;
                    let mut maxb6ks: f64 = 0.0;
                    // Looping over all particles and all 3 components of the acceleration. 
                    for (i, particle) in self.universe.particles[..self.universe.n_particles].iter().enumerate() {
                        let v2 = particle.inertial_velocity.x*particle.inertial_velocity.x
                                    +particle.inertial_velocity.y*particle.inertial_velocity.y
                                    +particle.inertial_velocity.z*particle.inertial_velocity.z;
                        let x2 = particle.inertial_position.x*particle.inertial_position.x
                                    +particle.inertial_position.y*particle.inertial_position.y
                                    +particle.inertial_position.z*particle.inertial_position.z;

                        if integrate_spin {
                            let dangular_momentum_dt2 = particle.dangular_momentum_dt.x*particle.dangular_momentum_dt.x
                                        +particle.dangular_momentum_dt.y*particle.dangular_momentum_dt.y
                                        +particle.dangular_momentum_dt.z*particle.dangular_momentum_dt.z;
                            let angular_momentum2 = particle.angular_momentum.x*particle.angular_momentum.x
                                        +particle.angular_momentum.y*particle.angular_momentum.y
                                        +particle.angular_momentum.z*particle.angular_momentum.z;
                            
                            // Skip slowly varying accelerations and angular momentums (spin)
                            if (v2*self.time_step*self.time_step/x2).abs() < 1e-16
                                && (dangular_momentum_dt2*self.time_step*self.time_step/angular_momentum2).abs() < 1e-16 {
                                continue;
                            }
                        } else {
                            // Skip slowly varying accelerations
                            if (v2*self.time_step*self.time_step/x2).abs() < 1e-16 {
                                continue;
                            }
                        }


                        for k in 3*i..3*(i+1) {
                            let ak  = self.at[k].abs();
                            if ak.is_normal() && ak>maxak {
                                maxak = ak;
                            }
                            let b6k = self.b_6[k].abs(); 
                            if b6k.is_normal() && b6k>maxb6k {
                                maxb6k = b6k;
                            }
                            //
                            if integrate_spin {
                                let dangular_momentum_dttk  = self.dangular_momentum_dtt[k].abs();
                                if dangular_momentum_dttk.is_normal() && dangular_momentum_dttk>maxak {
                                    max_dangular_momentum_dtk = dangular_momentum_dttk;
                                }
                                let b6ks = self.sb_6[k].abs(); 
                                if b6ks.is_normal() && b6ks>maxb6ks {
                                    maxb6ks = b6ks;
                                }
                            }
                        }
                    }
                    integrator_error = self.max(maxb6k, maxb6ks)/self.max(maxak, max_dangular_momentum_dtk);
                } else {
                    for k in 0..3*self.n_particles {
                        let mut integrator_error_a = 0.;
                        let mut integrator_error_s = 0.;
                        let ak  = self.at[k];
                        let b6k = self.b_6[k]; 
                        let errork = (b6k/ak).abs();
                        if errork.is_normal() && errork>integrator_error {
                            integrator_error_a = errork;
                        }
                        //
                        if integrate_spin {
                            let dangular_momentum_dttk  = self.dangular_momentum_dtt[k];
                            let b6ks = self.sb_6[k]; 
                            let errorks = (b6ks/dangular_momentum_dttk).abs();
                            if errorks.is_normal() && errorks>integrator_error {
                                integrator_error_s = errorks;
                            }
                        }
                        integrator_error = self.max(integrator_error_a, integrator_error_s);
                    }
                }

                let mut dt_new: f64;
                if  integrator_error.is_normal() { 	
                    // if error estimate is available, then increase by more educated guess
                    //dt_new = (INTEGRATOR_EPSILON/integrator_error).powf(1./7.) * dt_done;
                    dt_new = self.sqrt7(INTEGRATOR_EPSILON/integrator_error) * dt_done;
                } else {
                    // In the rare case that the error estimate doesn't give a finite number (e.g. when all forces accidentally cancel up to machine precission).
                    dt_new = dt_done/SAFETY_FACTOR; // by default, increase timestep a little
                }
                
                if dt_new.abs() < INTEGRATOR_MIN_DT {
                    //// copysignf is unstable and it cannot be used with rust 1.0
                    //unsafe {
                        //dt_new = std::intrinsics::copysignf64(INTEGRATOR_MIN_DT, dt_new);
                    //}
                    if dt_new > 0. {
                        dt_new = INTEGRATOR_MIN_DT.abs();
                    } else {
                        dt_new = -1. * INTEGRATOR_MIN_DT.abs();
                    }
                }

                if INTEGRATOR_MAX_DT > 0. && dt_new.abs() > INTEGRATOR_MAX_DT {
                    dt_new = INTEGRATOR_MAX_DT.abs();
                }
                
                if (dt_new/dt_done).abs() < SAFETY_FACTOR {	// New timestep is significantly smaller.
                    // Reset particles
                    for (k, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                        particle.inertial_position.x = self.x0[3*k+0];	// Set inital position
                        particle.inertial_position.y = self.x0[3*k+1];
                        particle.inertial_position.z = self.x0[3*k+2];

                        particle.inertial_velocity.x = self.v0[3*k+0];	// Set inital velocity
                        particle.inertial_velocity.y = self.v0[3*k+1];
                        particle.inertial_velocity.z = self.v0[3*k+2];

                        particle.inertial_acceleration.x = self.a0[3*k+0];	// Set inital acceleration
                        particle.inertial_acceleration.y = self.a0[3*k+1];
                        particle.inertial_acceleration.z = self.a0[3*k+2];

                        if integrate_spin {
                            particle.dangular_momentum_dt.x = self.dangular_momentum_dt0[3*k+0];	// Set inital dangular_momentum_dt
                            particle.dangular_momentum_dt.y = self.dangular_momentum_dt0[3*k+1];
                            particle.dangular_momentum_dt.z = self.dangular_momentum_dt0[3*k+2];

                            particle.angular_momentum.x = self.angular_momentum0[3*k+0];	// Set inital angular_momentum
                            particle.angular_momentum.y = self.angular_momentum0[3*k+1];
                            particle.angular_momentum.z = self.angular_momentum0[3*k+2];
                        }
                    }
                    self.time_step = dt_new;
                    if self.time_step_last_success != 0. {		// Do not predict next self.e/self.b values if this is the first time step.
                        let ratio = self.time_step/self.time_step_last_success;
                        self.predict_next_step(ratio);
                    }
                    
                    //println!("Step rejected, repeating with new timestep {}", self.time_step);
                    continue; // Step rejected. Do again. 
                }
                if (dt_new/dt_done).abs() > 1.0 {	// New timestep is larger.
                    if dt_new/dt_done > 1./SAFETY_FACTOR {
                        //println!("Time step is too large ({:.15}), do not increase so much", dt_new);
                        dt_new = dt_done / SAFETY_FACTOR;	// Don't increase the timestep by too much compared to the last one.
                    }
                }
                self.time_step = dt_new;
            }
            ////////////////////////////////////////////////////////////////////
            //// END: Find new timestep
            ////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////
            // Find new position and velocity values at end of the sequence 
            // (Eqs. 11, 12 of Everhart)
            ////////////////////////////////////////////////////////////////////
            let dt_done2 = dt_done * dt_done;
            for k in 0..3*self.n_particles {
                {
                    let (tmp_x0, tmp_csx) = self.add_cs(self.x0[k], self.csx[k], self.b_6[k]/72.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.b_5[k]/56.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.b_4[k]/42.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.b_3[k]/30.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.b_2[k]/20.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.b_1[k]/12.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.b_0[k]/6.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.a0[k]/2.*dt_done2);
                    let (tmp_x0, tmp_csx) = self.add_cs(tmp_x0, tmp_csx, self.v0[k]*dt_done);
                    self.csx[k] = tmp_csx;
                    self.x0[k] = tmp_x0;
                }
                {
                    let (tmp_v0, tmp_csv) = self.add_cs(self.v0[k], self.csv[k], self.b_6[k]/8.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.b_5[k]/7.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.b_4[k]/6.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.b_3[k]/5.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.b_2[k]/4.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.b_1[k]/3.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.b_0[k]/2.*dt_done);
                    let (tmp_v0, tmp_csv) = self.add_cs(tmp_v0, tmp_csv, self.a0[k]*dt_done);
                    self.csv[k] = tmp_csv;
                    self.v0[k] = tmp_v0;

                    if integrate_spin {
                        // spin
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(self.angular_momentum0[k], self.css[k], self.sb_6[k]/8.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.sb_5[k]/7.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.sb_4[k]/6.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.sb_3[k]/5.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.sb_2[k]/4.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.sb_1[k]/3.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.sb_0[k]/2.*dt_done);
                        let (tmp_angular_momentum0, tmp_css) = self.add_cs(tmp_angular_momentum0, tmp_css, self.dangular_momentum_dt0[k]*dt_done);
                        self.css[k] = tmp_css;
                        self.angular_momentum0[k] = tmp_angular_momentum0;
                    }
                }
            }

            self.current_time += dt_done;
            self.time_step_last_success = dt_done;

            // Swap particle buffers
            for (k, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                particle.inertial_position.x = self.x0[3*k+0];	// Set final position
                particle.inertial_position.y = self.x0[3*k+1];
                particle.inertial_position.z = self.x0[3*k+2];

                particle.inertial_velocity.x = self.v0[3*k+0];	// Set final velocity
                particle.inertial_velocity.y = self.v0[3*k+1];
                particle.inertial_velocity.z = self.v0[3*k+2];

                if integrate_spin {
                    particle.angular_momentum.x = self.angular_momentum0[3*k+0];  // Set final angular momentum
                    particle.angular_momentum.y = self.angular_momentum0[3*k+1];
                    particle.angular_momentum.z = self.angular_momentum0[3*k+2];
                } 
            }

            //self.copybuffers(&mut self.e, &mut self.er);		
            //self.copybuffers(&mut self.b, &mut self.br);		
            //self.er = self.e.clone();
            self.er_0 = self.e_0.clone();
            self.er_1 = self.e_1.clone();
            self.er_2 = self.e_2.clone();
            self.er_3 = self.e_3.clone();
            self.er_4 = self.e_4.clone();
            self.er_5 = self.e_5.clone();
            self.er_6 = self.e_6.clone();
            //self.br = self.b.clone();
            self.br_0 = self.b_0.clone();
            self.br_1 = self.b_1.clone();
            self.br_2 = self.b_2.clone();
            self.br_3 = self.b_3.clone();
            self.br_4 = self.b_4.clone();
            self.br_5 = self.b_5.clone();
            self.br_6 = self.b_6.clone();
            //self.ser = self.se.clone();
            self.ser_0 = self.se_0.clone();
            self.ser_1 = self.se_1.clone();
            self.ser_2 = self.se_2.clone();
            self.ser_3 = self.se_3.clone();
            self.ser_4 = self.se_4.clone();
            self.ser_5 = self.se_5.clone();
            self.ser_6 = self.se_6.clone();
            //self.sbr = self.sb.clone();
            self.sbr_0 = self.sb_0.clone();
            self.sbr_1 = self.sb_1.clone();
            self.sbr_2 = self.sb_2.clone();
            self.sbr_3 = self.sb_3.clone();
            self.sbr_4 = self.sb_4.clone();
            self.sbr_5 = self.sb_5.clone();
            self.sbr_6 = self.sb_6.clone();
            let ratio = self.time_step/dt_done;
            self.predict_next_step(ratio);
            break; // Success.

        } // end main loop
    }


    fn predict_next_step(&mut self, ratio: f64) {
        let mut consider_dangular_momentum_dt_from_general_relativity = false;
        if self.universe.consider_effects.general_relativity && self.universe.general_relativity_implementation == GeneralRelativityImplementation::Kidder1995 {
            consider_dangular_momentum_dt_from_general_relativity = true;
        }
        let integrate_spin = self.universe.consider_effects.tides || self.universe.consider_effects.rotational_flattening
                            || self.universe.consider_effects.evolution || consider_dangular_momentum_dt_from_general_relativity;
        if ratio > 20. {
            // Do not predict if stepsize increase is very large.
            for k in 0..3*self.n_particles {
                self.e_0[k] = 0.;
                self.e_1[k] = 0.;
                self.e_2[k] = 0.;
                self.e_3[k] = 0.;
                self.e_4[k] = 0.;
                self.e_5[k] = 0.;
                self.e_6[k] = 0.;
                self.b_0[k] = 0.;
                self.b_1[k] = 0.;
                self.b_2[k] = 0.;
                self.b_3[k] = 0.;
                self.b_4[k] = 0.;
                self.b_5[k] = 0.;
                self.b_6[k] = 0.;
                if integrate_spin {
                    self.se_0[k] = 0.;
                    self.se_1[k] = 0.;
                    self.se_2[k] = 0.;
                    self.se_3[k] = 0.;
                    self.se_4[k] = 0.;
                    self.se_5[k] = 0.;
                    self.se_6[k] = 0.;
                    self.sb_0[k] = 0.;
                    self.sb_1[k] = 0.;
                    self.sb_2[k] = 0.;
                    self.sb_3[k] = 0.;
                    self.sb_4[k] = 0.;
                    self.sb_5[k] = 0.;
                    self.sb_6[k] = 0.;
                }
            }
        } else {
            // Predict new B values to use at the start of the next sequence. The predicted
            // values from the last call are saved as E. The correction, BD, between the
            // actual and predicted values of B is applied in advance as a correction.
            //
            let q1 = ratio;
            let q2 = q1 * q1;
            let q3 = q1 * q2;
            let q4 = q2 * q2;
            let q5 = q2 * q3;
            let q6 = q3 * q3;
            let q7 = q3 * q4;

            for k in 0..3*self.n_particles {
                let be0 = self.br_0[k] - self.er_0[k];
                let be1 = self.br_1[k] - self.er_1[k];
                let be2 = self.br_2[k] - self.er_2[k];
                let be3 = self.br_3[k] - self.er_3[k];
                let be4 = self.br_4[k] - self.er_4[k];
                let be5 = self.br_5[k] - self.er_5[k];
                let be6 = self.br_6[k] - self.er_6[k];

                // Estimate B values for the next sequence (Eqs. 13 of Everhart).
                self.e_0[k] = q1*(self.br_6[k]* 7.0 + self.br_5[k]* 6.0 + self.br_4[k]* 5.0 + self.br_3[k]* 4.0 + self.br_2[k]* 3.0 + self.br_1[k]*2.0 + self.br_0[k]);
                self.e_1[k] = q2*(self.br_6[k]*21.0 + self.br_5[k]*15.0 + self.br_4[k]*10.0 + self.br_3[k]* 6.0 + self.br_2[k]* 3.0 + self.br_1[k]);
                self.e_2[k] = q3*(self.br_6[k]*35.0 + self.br_5[k]*20.0 + self.br_4[k]*10.0 + self.br_3[k]* 4.0 + self.br_2[k]);
                self.e_3[k] = q4*(self.br_6[k]*35.0 + self.br_5[k]*15.0 + self.br_4[k]* 5.0 + self.br_3[k]);
                self.e_4[k] = q5*(self.br_6[k]*21.0 + self.br_5[k]* 6.0 + self.br_4[k]);
                self.e_5[k] = q6*(self.br_6[k]* 7.0 + self.br_5[k]);
                self.e_6[k] = q7* self.br_6[k];
                
                self.b_0[k] = self.e_0[k] + be0;
                self.b_1[k] = self.e_1[k] + be1;
                self.b_2[k] = self.e_2[k] + be2;
                self.b_3[k] = self.e_3[k] + be3;
                self.b_4[k] = self.e_4[k] + be4;
                self.b_5[k] = self.e_5[k] + be5;
                self.b_6[k] = self.e_6[k] + be6;

                if integrate_spin {
                    // Spin
                    let sbe0 = self.sbr_0[k] - self.ser_0[k];
                    let sbe1 = self.sbr_1[k] - self.ser_1[k];
                    let sbe2 = self.sbr_2[k] - self.ser_2[k];
                    let sbe3 = self.sbr_3[k] - self.ser_3[k];
                    let sbe4 = self.sbr_4[k] - self.ser_4[k];
                    let sbe5 = self.sbr_5[k] - self.ser_5[k];
                    let sbe6 = self.sbr_6[k] - self.ser_6[k];

                    // Estimate B values for the next sequence (Eqs. 13 of Everhart).
                    self.se_0[k] = q1*(self.sbr_6[k]* 7.0 + self.sbr_5[k]* 6.0 + self.sbr_4[k]* 5.0 + self.sbr_3[k]* 4.0 + self.sbr_2[k]* 3.0 + self.sbr_1[k]*2.0 + self.sbr_0[k]);
                    self.se_1[k] = q2*(self.sbr_6[k]*21.0 + self.sbr_5[k]*15.0 + self.sbr_4[k]*10.0 + self.sbr_3[k]* 6.0 + self.sbr_2[k]* 3.0 + self.sbr_1[k]);
                    self.se_2[k] = q3*(self.sbr_6[k]*35.0 + self.sbr_5[k]*20.0 + self.sbr_4[k]*10.0 + self.sbr_3[k]* 4.0 + self.sbr_2[k]);
                    self.se_3[k] = q4*(self.sbr_6[k]*35.0 + self.sbr_5[k]*15.0 + self.sbr_4[k]* 5.0 + self.sbr_3[k]);
                    self.se_4[k] = q5*(self.sbr_6[k]*21.0 + self.sbr_5[k]* 6.0 + self.sbr_4[k]);
                    self.se_5[k] = q6*(self.sbr_6[k]* 7.0 + self.sbr_5[k]);
                    self.se_6[k] = q7* self.sbr_6[k];
                    
                    self.sb_0[k] = self.se_0[k] + sbe0;
                    self.sb_1[k] = self.se_1[k] + sbe1;
                    self.sb_2[k] = self.se_2[k] + sbe2;
                    self.sb_3[k] = self.se_3[k] + sbe3;
                    self.sb_4[k] = self.se_4[k] + sbe4;
                    self.sb_5[k] = self.se_5[k] + sbe5;
                    self.sb_6[k] = self.se_6[k] + sbe6;
                }
            }
        }
    }

    fn sqrt7(&self, a: f64) -> f64 {
        // Machine independent implementation of powf(1./7.)
        let mut x: f64 = 1.;
        for _k in 0..20 { // A smaller number should be ok too.
            let x6 = x*x*x*x*x*x;
            x += (a/x6-x)/7.;
        }
        x
    }

    fn max(&self, a: f64, b: f64) -> f64 {
        if a > b {
            a
        } else {
            b
        }
    }

    fn add_cs(&self, p: f64, csp: f64, inp: f64) -> (f64, f64) {
        let y = inp - csp;
        let new_p = p + y;
        let new_csp = (new_p - p) - y;
        (new_p, new_csp)
    }

}
