extern crate time;
use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use super::Integrator;
use super::super::constants::{INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT, INTEGRATOR_EPSILON, INTEGRATOR_EPSILON_GLOBAL, INTEGRATOR_MIN_DT, SAFETY_FACTOR};
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

///https://arxiv.org/abs/1409.4779
///IAS15: A fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine
///precision over a billion orbits
///
/// https://arxiv.org/pdf/1110.4876v2.pdf
/// variable time-steps also
/// break the symplectic nature of an integrator.

#[derive(Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub struct Ias15 {
    time_step: f64,
    pub universe: Universe,
    current_time: f64,
    current_iteration: u32,
    recovery_snapshot_period: f64,
    historic_snapshot_period: f64,
    last_recovery_snapshot_time: f64,
    last_historic_snapshot_time: f64,
    pub n_historic_snapshots: usize,
    pub hash: u64,
    //// Integrator IAS15 data:
    n_particles: usize,
    integrator_iterations_max_exceeded : i32,  // Count how many times the iteration did not converge
    time_step_last_success: f64,			// Last accepted timestep (corresponding to br and er)
    b: Vec<Vec<f64>>, // Coefficient b: acceleration dimension
    br: Vec<Vec<f64>>, // Previous b
    g: Vec<Vec<f64>>, // Coefficient g (it can also be expressed in terms of b)
    e: Vec<Vec<f64>>,
    er: Vec<Vec<f64>>, // Previous g
    at: Vec<f64>, // Temporary buffer for acceleration
    x0: Vec<f64>, // Temporary buffer for position (used for initial values at h=0)
    v0: Vec<f64>, // Temporary buffer for velcoity (used for initial values at h=0)
    a0: Vec<f64>, // Temporary buffer for acceleration (used for initial values at h=0)
    // Compensated summation coefficients
    csx : Vec<f64>,
    csv : Vec<f64>,
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
        let mut universe_integrator = Ias15 {
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
                    b :   (0..7).map(|_| vec![0.; 3*n_particles]).collect(),
                    g  :  (0..7).map(|_| vec![0.; 3*n_particles]).collect(),
                    e  :  (0..7).map(|_| vec![0.; 3*n_particles]).collect(),
                    br :  (0..7).map(|_| vec![0.; 3*n_particles]).collect(),
                    er :  (0..7).map(|_| vec![0.; 3*n_particles]).collect(),
                    at  : vec![0.; 3*n_particles],
                    x0  : vec![0.; 3*n_particles],
                    v0  : vec![0.; 3*n_particles],
                    a0  : vec![0.; 3*n_particles],
                    csx : vec![0.; 3*n_particles],
                    csv : vec![0.; 3*n_particles],
                    s   : [0.; 9],
                    };
        // Initialize physical values
        let current_time = 0.;
        universe_integrator.universe.calculate_norm_spin(); // Needed for evolution
        universe_integrator.universe.calculate_particles_evolving_quantities(current_time); // Make sure we start with the good initial values
        universe_integrator
    }

    pub fn restore_snapshot(universe_integrator_snapshot_path: &Path, verify_integrity: bool) -> Result<Ias15, String> {
        let mut universe_integrator: Ias15;
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

impl Integrator for Ias15 {


    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.current_time;
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            if ! silent_mode {
                print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
                let _ = std::io::stdout().flush();
            }
        }

        // Calculate accelerations.
        let integrator_is_whfasthelio = false;
        self.universe.gravity_calculate_acceleration(integrator_is_whfasthelio);
        // Calculate non-gravity accelerations.
        let evolution = true;
        let dspin_dt = true;
        let accelerations = true;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations);

        self.integrator();
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

impl Ias15 {

    #[allow(dead_code)]
    fn integrator(&mut self) {
    
        // Gauss-Radau spacings for substeps within a sequence, for the 15th order 
        // integrator. The sum of the h values should be 3.733333333333333
        let h : [f64; 8]	= [ 0.0, 0.0562625605369221464656521910, 0.1802406917368923649875799428, 0.3526247171131696373739077702, 0.5471536263305553830014485577, 0.7342101772154105410531523211, 0.8853209468390957680903597629, 0.9775206135612875018911745004]; 

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
        let r : [f64; 28] = [0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105410531523, 0.6779476166784883945875001, 0.5539694854785181760655724, 0.3815854601022409036792446, 0.1870565508848551580517038, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852270372074, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035946, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769608380222, 0.0921996667221917338008147];

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
        let c : [f64; 21] = [-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915610919, 0.0421585277212687082291130, -0.3600995965020568162530901, 1.2501507118406910366792415, -1.8704917729329500728817408, 0.0012717903090268677658020, -0.0387603579159067708505249, 0.3609622434528459872559689, -1.4668842084004269779203515, 2.9061362593084293206895457, -2.7558127197720458409721005];
        let d : [f64; 21] = [0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500728817408, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691598182153473, 2.7558127197720458409721005];




        loop {

            for (k, particle) in self.universe.particles[..self.universe.n_particles].iter().enumerate() {
                self.x0[3*k]   = particle.position.x;
                self.x0[3*k+1] = particle.position.y;
                self.x0[3*k+2] = particle.position.z;
                self.v0[3*k]   = particle.velocity.x;
                self.v0[3*k+1] = particle.velocity.y;
                self.v0[3*k+2] = particle.velocity.z;
                self.a0[3*k]   = particle.acceleration.x;
                self.a0[3*k+1] = particle.acceleration.y;  
                self.a0[3*k+2] = particle.acceleration.z;
            }

            // Find g values from b values predicted at the last call (Eqs. 7 of Everhart)
            for k in 0..3*self.n_particles {
                self.g[0][k] = self.b[6][k]*d[15] + self.b[5][k]*d[10] + self.b[4][k]*d[6] + self.b[3][k]*d[3]  + self.b[2][k]*d[1]  + self.b[1][k]*d[0]  + self.b[0][k];
                self.g[1][k] = self.b[6][k]*d[16] + self.b[5][k]*d[11] + self.b[4][k]*d[7] + self.b[3][k]*d[4]  + self.b[2][k]*d[2]  + self.b[1][k];
                self.g[2][k] = self.b[6][k]*d[17] + self.b[5][k]*d[12] + self.b[4][k]*d[8] + self.b[3][k]*d[5]  + self.b[2][k];
                self.g[3][k] = self.b[6][k]*d[18] + self.b[5][k]*d[13] + self.b[4][k]*d[9] + self.b[3][k];
                self.g[4][k] = self.b[6][k]*d[19] + self.b[5][k]*d[14] + self.b[4][k];
                self.g[5][k] = self.b[6][k]*d[20] + self.b[5][k];
                self.g[6][k] = self.b[6][k];
            }

            let t_beginning = self.current_time;
            let mut predictor_corrector_error: f64 = 1e10;
            //let mut predictor_corrector_error: f64 = 1e300;
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
                        println!("[WARNING {} UTC] At least {} predictor corrector loops in integrator IAS15 did not converge. This is typically an indication of the timestep being too large.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), INTEGRATOR_ITERATIONS_WARNING);
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
                    //println!("Interval {} with substep {} of {}", n, self.s[0], self.time_step);

                    //// Prepare particles arrays for force calculation
                    // Predict positions at interval n using self.b values
                    for (i, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                        let k0 : usize = 3*i+0;
                        let k1 : usize = 3*i+1;
                        let k2 : usize = 3*i+2;

                        // Equation 7 in paper 2015MNRAS.446.1424R
                        let xk0: f64  = self.csx[k0] + (self.s[8]*self.b[6][k0] + self.s[7]*self.b[5][k0] + self.s[6]*self.b[4][k0] + self.s[5]*self.b[3][k0] + self.s[4]*self.b[2][k0] + self.s[3]*self.b[1][k0] + self.s[2]*self.b[0][k0] + self.s[1]*self.a0[k0] + self.s[0]*self.v0[k0] );
                        particle.position.x = xk0 + self.x0[k0];

                        let xk1: f64  = self.csx[k1] + (self.s[8]*self.b[6][k1] + self.s[7]*self.b[5][k1] + self.s[6]*self.b[4][k1] + self.s[5]*self.b[3][k1] + self.s[4]*self.b[2][k1] + self.s[3]*self.b[1][k1] + self.s[2]*self.b[0][k1] + self.s[1]*self.a0[k1] + self.s[0]*self.v0[k1] );
                        particle.position.y = xk1 + self.x0[k1];

                        let xk2: f64 = self.csx[k2] + (self.s[8]*self.b[6][k2] + self.s[7]*self.b[5][k2] + self.s[6]*self.b[4][k2] + self.s[5]*self.b[3][k2] + self.s[4]*self.b[2][k2] + self.s[3]*self.b[1][k2] + self.s[2]*self.b[0][k2] + self.s[1]*self.a0[k2] + self.s[0]*self.v0[k2] );
                        particle.position.z = xk2 + self.x0[k2];
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
                            let vk0 =  self.csv[k0] + self.s[7]*self.b[6][k0] + self.s[6]*self.b[5][k0] + self.s[5]*self.b[4][k0] + self.s[4]*self.b[3][k0] + self.s[3]*self.b[2][k0] + self.s[2]*self.b[1][k0] + self.s[1]*self.b[0][k0] + self.s[0]*self.a0[k0];
                            particle.velocity.x = vk0 + self.v0[k0];
                            let vk1 =  self.csv[k1] + self.s[7]*self.b[6][k1] + self.s[6]*self.b[5][k1] + self.s[5]*self.b[4][k1] + self.s[4]*self.b[3][k1] + self.s[3]*self.b[2][k1] + self.s[2]*self.b[1][k1] + self.s[1]*self.b[0][k1] + self.s[0]*self.a0[k1];
                            particle.velocity.y = vk1 + self.v0[k1];
                            let vk2 =  self.csv[k2] + self.s[7]*self.b[6][k2] + self.s[6]*self.b[5][k2] + self.s[5]*self.b[4][k2] + self.s[4]*self.b[3][k2] + self.s[3]*self.b[2][k2] + self.s[2]*self.b[1][k2] + self.s[1]*self.b[0][k2] + self.s[0]*self.a0[k2];
                            particle.velocity.z = vk2 + self.v0[k2];
                        }
                    }


                    // Calculate accelerations.
                    let integrator_is_whfasthelio = false;
                    self.universe.gravity_calculate_acceleration(integrator_is_whfasthelio);
                    // Calculate non-gravity accelerations.
                    let evolution = true;
                    let dspin_dt = true;
                    let accelerations = true;
                    self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations);

                    for (k, particle) in self.universe.particles[..self.universe.n_particles].iter().enumerate() {
                        self.at[3*k]   = particle.acceleration.x;
                        self.at[3*k+1] = particle.acceleration.y;  
                        self.at[3*k+2] = particle.acceleration.z;
                    }

                    // Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
                    match n {
                        1 => {
                                for k in 0..3*self.n_particles {
                                    let tmp = self.g[0][k];
                                    self.g[0][k]  = (self.at[k] - self.a0[k]) / r[0];
                                    self.b[0][k] += self.g[0][k] - tmp;
                                }
                            },
                        2 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g[1][k];
                                    let gk = self.at[k] - self.a0[k];
                                    self.g[1][k] = (gk/r[1] - self.g[0][k])/r[2];
                                    tmp = self.g[1][k] - tmp;
                                    self.b[0][k] += tmp * c[0];
                                    self.b[1][k] += tmp;
                                }
                            },
                        3 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g[2][k];
                                    let gk = self.at[k] - self.a0[k];
                                    self.g[2][k] = ((gk/r[3] - self.g[0][k])/r[4] - self.g[1][k])/r[5];
                                    tmp = self.g[2][k] - tmp;
                                    self.b[0][k] += tmp * c[1];
                                    self.b[1][k] += tmp * c[2];
                                    self.b[2][k] += tmp;
                                }
                            },
                        4 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g[3][k];
                                    let gk = self.at[k] - self.a0[k];
                                    self.g[3][k] = (((gk/r[6] - self.g[0][k])/r[7] - self.g[1][k])/r[8] - self.g[2][k])/r[9];
                                    tmp = self.g[3][k] - tmp;
                                    self.b[0][k] += tmp * c[3];
                                    self.b[1][k] += tmp * c[4];
                                    self.b[2][k] += tmp * c[5];
                                    self.b[3][k] += tmp;
                                }
                            },
                        5 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g[4][k];
                                    let gk = self.at[k] - self.a0[k];
                                    self.g[4][k] = ((((gk/r[10] - self.g[0][k])/r[11] - self.g[1][k])/r[12] - self.g[2][k])/r[13] - self.g[3][k])/r[14];
                                    tmp = self.g[4][k] - tmp;
                                    self.b[0][k] += tmp * c[6];
                                    self.b[1][k] += tmp * c[7];
                                    self.b[2][k] += tmp * c[8];
                                    self.b[3][k] += tmp * c[9];
                                    self.b[4][k] += tmp;
                                }
                            },
                        6 => {
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g[5][k];
                                    let gk = self.at[k] - self.a0[k];
                                    self.g[5][k] = (((((gk/r[15] - self.g[0][k])/r[16] - self.g[1][k])/r[17] - self.g[2][k])/r[18] - self.g[3][k])/r[19] - self.g[4][k])/r[20];
                                    tmp = self.g[5][k] - tmp;
                                    self.b[0][k] += tmp * c[10];
                                    self.b[1][k] += tmp * c[11];
                                    self.b[2][k] += tmp * c[12];
                                    self.b[3][k] += tmp * c[13];
                                    self.b[4][k] += tmp * c[14];
                                    self.b[5][k] += tmp;
                                }
                            },
                        7 => {
                                let mut maxak: f64 = 0.0;
                                let mut maxb6ktmp: f64 = 0.0;
                                for k in 0..3*self.n_particles {
                                    let mut tmp = self.g[6][k];
                                    let gk = self.at[k] - self.a0[k];
                                    self.g[6][k] = ((((((gk/r[21] - self.g[0][k])/r[22] - self.g[1][k])/r[23] - self.g[2][k])/r[24] - self.g[3][k])/r[25] - self.g[4][k])/r[26] - self.g[5][k])/r[27];
                                    tmp = self.g[6][k] - tmp;	
                                    self.b[0][k] += tmp * c[15];
                                    self.b[1][k] += tmp * c[16];
                                    self.b[2][k] += tmp * c[17];
                                    self.b[3][k] += tmp * c[18];
                                    self.b[4][k] += tmp * c[19];
                                    self.b[5][k] += tmp * c[20];
                                    self.b[6][k] += tmp;
                                    
                                    // Monitor change in self.b[6][k] relative to self.at[k]. The predictor corrector scheme is converged if it is close to 0.
                                    if INTEGRATOR_EPSILON_GLOBAL {
                                        let ak  = self.at[k].abs();
                                        if ak.is_normal() && ak>maxak {
                                            maxak = ak;
                                        }
                                        let b6ktmp = tmp.abs();  // change of b6ktmp coefficient
                                        if b6ktmp.is_normal() && b6ktmp>maxb6ktmp {
                                            maxb6ktmp = b6ktmp;
                                        }
                                    } else {
                                        let ak  = self.at[k];
                                        let b6ktmp = tmp; 
                                        let errork = (b6ktmp/ak).abs();
                                        if errork.is_normal() && errork>predictor_corrector_error {
                                            predictor_corrector_error = errork;
                                        }
                                    }
                                } 
                                if INTEGRATOR_EPSILON_GLOBAL {
                                    predictor_corrector_error = maxb6ktmp/maxak;
                                }
                            },
                        _ => { println!("[WARNING {} UTC] This should not happen because the loop stops at 7!", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap()); }

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
                // INTEGRATOR_EPSILON_GLOBAL==true
                //   Here, the fractional error is calculated for each particle individually and we use the maximum of the fractional error.
                //   This might fail in cases where a particle does not experience any (physical) acceleration besides roundoff errors. 
                let mut integrator_error: f64 = 0.0;
                if INTEGRATOR_EPSILON_GLOBAL {
                    let mut maxak: f64 = 0.0;
                    let mut maxb6k: f64 = 0.0;
                    // Looping over all particles and all 3 components of the acceleration. 
                    for (i, particle) in self.universe.particles[..self.universe.n_particles].iter().enumerate() {
                        let v2 = particle.velocity.x*particle.velocity.x
                                    +particle.velocity.y*particle.velocity.y
                                    +particle.velocity.z*particle.velocity.z;
                        let x2 = particle.position.x*particle.position.x
                                    +particle.position.y*particle.position.y
                                    +particle.position.z*particle.position.z;

                        // Skip slowly varying accelerations
                        if (v2*self.time_step*self.time_step/x2).abs() < 1e-16 {
                            continue;
                        }

                        for k in 3*i..3*(i+1) {
                            let ak  = self.at[k].abs();
                            if ak.is_normal() && ak>maxak {
                                maxak = ak;
                            }
                            let b6k = self.b[6][k].abs(); 
                            if b6k.is_normal() && b6k>maxb6k {
                                maxb6k = b6k;
                            }
                        }
                    }
                    integrator_error = maxb6k/maxak;
                } else {
                    for k in 0..3*self.n_particles {
                        let ak  = self.at[k];
                        let b6k = self.b[6][k]; 
                        let errork = (b6k/ak).abs();
                        if errork.is_normal() && errork>integrator_error {
                            integrator_error = errork;
                        }
                    }
                }

                let mut dt_new: f64;
                if  integrator_error.is_normal() { 	
                    // if error estimate is available, then increase by more educated guess
                    dt_new = (INTEGRATOR_EPSILON/integrator_error).powf(1./7.) * dt_done;
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
                
                if (dt_new/dt_done).abs() < SAFETY_FACTOR {	// New timestep is significantly smaller.
                    // Reset particles
                    for (k, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                        particle.position.x = self.x0[3*k+0];	// Set inital position
                        particle.position.y = self.x0[3*k+1];
                        particle.position.z = self.x0[3*k+2];

                        particle.velocity.x = self.v0[3*k+0];	// Set inital velocity
                        particle.velocity.y = self.v0[3*k+1];
                        particle.velocity.z = self.v0[3*k+2];
                    }
                    self.time_step = dt_new;
                    if self.time_step_last_success != 0. {		// Do not predict next self.e/self.b values if this is the first time step.
                        let ratio = self.time_step/self.time_step_last_success;
                        //self.predict_next_step(&mut self.e, &mut self.b, self.er, self.br, ratio);
                        self.predict_next_step(ratio);
                    }
                    
                    continue; // Step rejected. Do again. 
                }		
                if (dt_new/dt_done).abs() > 1.0 {	// New timestep is larger.
                    if dt_new/dt_done > 1./SAFETY_FACTOR {
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
                    let a = self.x0[k];
                    self.csx[k]  +=  (self.b[6][k]/72. + self.b[5][k]/56. + self.b[4][k]/42. 
                                    + self.b[3][k]/30. + self.b[2][k]/20. + self.b[1][k]/12. + self.b[0][k]/6. + self.a0[k]/2.) 
                                    * dt_done2 + self.v0[k] * dt_done;
                    self.x0[k]    = a + self.csx[k];
                    self.csx[k]  += a - self.x0[k]; 
                }
                {
                    let a = self.v0[k]; 
                    self.csv[k]  += (self.b[6][k]/8. + self.b[5][k]/7. + self.b[4][k]/6. 
                                    + self.b[3][k]/5. + self.b[2][k]/4. + self.b[1][k]/3. + self.b[0][k]/2. + self.a0[k])
                                    * dt_done;
                    self.v0[k]    = a + self.csv[k];
                    self.csv[k]  += a - self.v0[k];
                }
            }

            self.current_time += dt_done;

            // Swap particle buffers
            for (k, particle) in self.universe.particles[..self.universe.n_particles].iter_mut().enumerate() {
                particle.position.x = self.x0[3*k+0];	// Set final position
                particle.position.y = self.x0[3*k+1];
                particle.position.z = self.x0[3*k+2];

                particle.velocity.x = self.v0[3*k+0];	// Set final velocity
                particle.velocity.y = self.v0[3*k+1];
                particle.velocity.z = self.v0[3*k+2];

                if particle.moment_of_inertia_ratio != 1. {
                    particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + self.time_step * particle.dspin_dt.x;
                    particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + self.time_step * particle.dspin_dt.y;
                    particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + self.time_step * particle.dspin_dt.z;
                } else {
                    particle.spin.x = particle.spin.x + self.time_step * particle.dspin_dt.x;
                    particle.spin.y = particle.spin.y + self.time_step * particle.dspin_dt.y;
                    particle.spin.z = particle.spin.z + self.time_step * particle.dspin_dt.z;
                }
                if particle.wind_factor != 0. {
                    // TODO: Verify wind factor
                    particle.spin.x += self.time_step * particle.wind_factor * particle.spin.x;
                    particle.spin.y += self.time_step * particle.wind_factor * particle.spin.y;
                    particle.spin.z += self.time_step * particle.wind_factor * particle.spin.z;
                }
            }

            self.time_step_last_success = dt_done;

            //self.copybuffers(&mut self.e, &mut self.er);		
            //self.copybuffers(&mut self.b, &mut self.br);		
            self.er = self.e.clone();
            self.br = self.b.clone();
            let ratio = self.time_step/dt_done;
            //self.predict_next_step(&mut self.e, &mut self.b, self.er, self.br, ratio);
            self.predict_next_step(ratio);
            break; // Success.

        } // end main loop
    }


    fn predict_next_step(&mut self, ratio: f64) {
    //fn predict_next_step(self, e : &mut [[f64; 3*N_PARTICLES]; 7], b : &mut [[f64; 3*N_PARTICLES]; 7], 
                            //_e : [[f64; 3*N_PARTICLES]; 7], _b: [[f64; 3*N_PARTICLES]; 7],
                            //ratio: f64) {
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
            let be0 = self.br[0][k] - self.er[0][k];
            let be1 = self.br[1][k] - self.er[1][k];
            let be2 = self.br[2][k] - self.er[2][k];
            let be3 = self.br[3][k] - self.er[3][k];
            let be4 = self.br[4][k] - self.er[4][k];
            let be5 = self.br[5][k] - self.er[5][k];
            let be6 = self.br[6][k] - self.er[6][k];

            // Estimate B values for the next sequence (Eqs. 13 of Everhart).
            self.e[0][k] = q1*(self.br[6][k]* 7.0 + self.br[5][k]* 6.0 + self.br[4][k]* 5.0 + self.br[3][k]* 4.0 + self.br[2][k]* 3.0 + self.br[1][k]*2.0 + self.br[0][k]);
            self.e[1][k] = q2*(self.br[6][k]*21.0 + self.br[5][k]*15.0 + self.br[4][k]*10.0 + self.br[3][k]* 6.0 + self.br[2][k]* 3.0 + self.br[1][k]);
            self.e[2][k] = q3*(self.br[6][k]*35.0 + self.br[5][k]*20.0 + self.br[4][k]*10.0 + self.br[3][k]* 4.0 + self.br[2][k]);
            self.e[3][k] = q4*(self.br[6][k]*35.0 + self.br[5][k]*15.0 + self.br[4][k]* 5.0 + self.br[3][k]);
            self.e[4][k] = q5*(self.br[6][k]*21.0 + self.br[5][k]* 6.0 + self.br[4][k]);
            self.e[5][k] = q6*(self.br[6][k]* 7.0 + self.br[5][k]);
            self.e[6][k] = q7* self.br[6][k];
            
            self.b[0][k] = self.e[0][k] + be0;
            self.b[1][k] = self.e[1][k] + be1;
            self.b[2][k] = self.e[2][k] + be2;
            self.b[3][k] = self.e[3][k] + be3;
            self.b[4][k] = self.e[4][k] + be4;
            self.b[5][k] = self.e[5][k] + be5;
            self.b[6][k] = self.e[6][k] + be6;
        }
    }

}
