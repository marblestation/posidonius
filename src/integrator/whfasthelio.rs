extern crate time;
use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use super::Integrator;
use super::super::constants::{PI, WHFAST_NMAX_QUART, WHFAST_NMAX_NEWT, MAX_PARTICLES};
use super::super::particles::Universe;
use super::super::particles::Particle;
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

/// WHFastHelio (symplectic integrator) to be used always in safe mode (always sync because it is required by tides)
/// and without correction (i.e. 2nd order integrator, comparable to mercury symplectic part of the hybrid integrator)
///
/// "The original WHFast algorithm described in Rein & Tamayo
/// (2015) was implemented in Jacobi coordinates. Jacobi coordinates
/// lead to a better precision compared to heliocentric
/// coordinates if orbits are well separated and do not cross each
/// other. If close encounter occur, then heliocentric coordinates
/// can help improve the integrator’s accuracy. For that reason we
/// implemented a heliocentric version of WHFast in REBOUND. We
/// call it WHFastHelio"
/// Source: HERMES: a hybrid integrator for simulating close encounters and planetesimal migration
///         Ari Silburt, Hanno Rein & Dan Tamayo
///
///  https://en.wikipedia.org/wiki/N-body_problem#Few_bodies
///  for N > 2, the N-body problem is chaotic,[37] which means that even small errors in
///  integration may grow exponentially in time. Third, a simulation may be over large stretches of
///  model time (e.g. millions of years) and numerical errors accumulate as integration time
///  increases.
///  There are a number of techniques to reduce errors in numerical integration.[18] Local
///  coordinate systems are used to deal with widely differing scales in some problems, for example
///  an Earth-Moon coordinate system in the context of a solar system simulation. Variational
///  methods and perturbation theory can yield approximate analytic trajectories upon which the
///  numerical integration can be a correction. The use of a symplectic integrator ensures that the
///  simulation obeys Hamilton's equations to a high degree of accuracy and in particular that
///  energy is conserved.
///
///  https://arxiv.org/pdf/1110.4876v2.pdf
///  These integrators are second order accurate
///  and symplectic, their symplectic nature is formally lost
///  as soon as self-gravity or collisions are approximated or when
///  velocity dependent forces are included.
///  NOTE: Which is the case for tidal forces
///
///  https://arxiv.org/pdf/1110.4876v2.pdf
///  A symplectic Wisdom-Holman mapping (WH, Wisdom
///  & Holman 1991) is implemented as a module in
///  integrator_wh.c. The implementation follows closely
///  that by the SWIFT code4. The WH mapping is a mixed variable
///  integrator that calculates the Keplerian motion of two bodies
///  orbiting each other exactly up to machine precision during the
///  drift sub-step. Thus, it is very accurate for problems in which
///  the particle motion is dominated by a central 1/r potential and
///  perturbations added in the kick sub-step are small. However,
///  the WH integrator is substantially slower than the leap-frog
///  integrator because Kepler’s equation is solved iteratively every
///  time-step for every particle.
///  The integrator assumes that the central object has the index 0
///  in the particle array, that it is located at the origin and that it does
///  not move. The coordinates of all particles are assumed to be the
///  heliocentric frame. During the sub-time-steps the coordinates are
///  converted to Jacobi coordinates (and back) according to their
///  index. The particle with index 1 has the first Jacobi index, and
///  so on. This works best if the particles are sorted according to
///  their semi-major axis. Note that this is not done automatically
///
/// https://arxiv.org/pdf/1506.01084v1.pdf
/// a complete reimplementation
/// of the Wisdom-Holman integrator. We show how to
/// speed up the algorithm in several ways and dramatically increase
/// its accuracy. Many of the improvements are related to finite double
/// floating-point precision on modern computers (IEEE754, ISO
/// 2011). The fact that almost all real numbers cannot be represented
/// exactly in floating-point precision leads to important consequences
/// for the numerical stability of any algorithm and the growth of numerical
/// round-off error.
/// To our knowledge, we present the first publicly available
/// Wisdom-Holman integrator that is unbiased, i.e. the errors are random
/// and uncorrelated. This leads to a very slow error growth. For
/// sufficiently small timesteps, we achieve Brouwer’s law, i.e., the energy
/// error grows as time to the power of one half.
/// We have also sped up the integrator through various improvements
/// to the integrator’s Kepler-solver. Our implementation allows
/// for the evolution of variational equations (to determine whether orbits
/// are chaotic) at almost no additional cost. Additionally, we implement
/// so-called symplectic correctors up to order eleven to increase
/// the accurary (Wisdom et al. 1996), allow for arbitrary unit
/// choices, and do not tie the integration to a particular frame of reference.

#[derive(Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub struct WHFastHelio {
    time_step: f64,
    half_time_step: f64,
    pub universe: Universe,
    last_spin: [Axes; MAX_PARTICLES], // For spin integration with the midpoint method
    current_time: f64,
    current_iteration: usize,
    recovery_snapshot_period: f64,
    historic_snapshot_period: f64,
    last_recovery_snapshot_time: f64,
    last_historic_snapshot_time: f64,
    pub n_historic_snapshots: usize,
    pub hash: u64,
    /**
     * @brief Heliocentric coordinates
     * @details This array contains the heliocentric coordinates of all particles.
     * It is automatically filled and updated by WHfastDemocratic.
     * Access this array with caution.
     */
    universe_heliocentric: Universe,

    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    ///< Flag to determine if current particle structure is synchronized
    is_synchronized: bool, 
    timestep_warning: usize ,
    set_to_center_of_mass: bool
}

impl Hash for WHFastHelio {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash only works with certain types (for instance, it does not work with f64 values)
        // thus we convert the whole integrator to a string thanks to the debug trait
        // and we hash that value
        format!("{:?}", self).hash(state); 
    }
}

impl WHFastHelio {
    pub fn new(time_step: f64, recovery_snapshot_period: f64, historic_snapshot_period: f64, universe: Universe) -> WHFastHelio {
        let universe_heliocentric = universe.clone(); // Clone will not clone evolving types
        let mut universe_integrator = WHFastHelio {
                    time_step:time_step,
                    half_time_step:0.5*time_step,
                    recovery_snapshot_period:recovery_snapshot_period,
                    historic_snapshot_period:historic_snapshot_period,
                    last_recovery_snapshot_time:-1.,
                    last_historic_snapshot_time:-1.,
                    n_historic_snapshots:0,
                    hash: 0,
                    universe:universe,
                    last_spin:[Axes{x:0., y:0., z:0. }; MAX_PARTICLES],
                    current_time:0.,
                    current_iteration:0,
                    // WHFastHelio specifics:
                    universe_heliocentric: universe_heliocentric,
                    is_synchronized: true,
                    timestep_warning: 0,
                    set_to_center_of_mass: false,
                    };
        // Initialize physical values
        let current_time = 0.;
        universe_integrator.universe.calculate_norm_spin(); // Needed for evolution
        universe_integrator.universe.calculate_particles_evolving_quantities(current_time); // Make sure we start with the good initial values
        universe_integrator
    }
    
    pub fn restore_snapshot(universe_integrator_snapshot_path: &Path, verify_integrity: bool) -> Result<WHFastHelio, String> {
        let mut universe_integrator: WHFastHelio;
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
                    panic!("[PANIC {} UTC] File '{}' has a zeroed hash (i.e., new simulation) but a current time different from zero ({})", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display(), universe_integrator.current_time)
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

impl Integrator for WHFastHelio {

    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            if self.set_to_center_of_mass {
                self.move_to_star_center();
            }
            self.universe.calculate_denergy_dt();
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.n_historic_snapshots as f64*self.historic_snapshot_period; // Instead of self.current_time to avoid small deviations
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            if ! silent_mode {
                print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
                let _ = std::io::stdout().flush();
            }
        }

        self.iterate_position_and_velocity_with_whfasthelio();
        if self.universe.evolving_particles_exist || self.universe.consider_tides || self.universe.consider_rotational_flattening {
            self.iterate_spin_with_midpoint_method();
        }

        // ---------------------------------------------------------------------
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
        self.hash = 0;
        let mut s = DefaultHasher::new();
        self.hash(&mut s);
        self.hash = s.finish();
    }

}


impl WHFastHelio {
    // WHFastHelio integrator
    fn iterate_position_and_velocity_with_whfasthelio(&mut self) {
        let time_step = self.time_step;
        let half_time_step = self.half_time_step;

        self.move_to_center_of_mass();
        // ---------------------------------------------------------------------
        // A 'DKD'-like integrator will do the first 'D' part.
        self.to_helio_posvel();
        self.helio_kepler_steps(half_time_step);
        self.to_inertial_posvel();
        self.current_time += self.half_time_step;
        // ---------------------------------------------------------------------

        // Calculate accelerations.
        let integrator_is_whfasthelio = true;
        self.universe.gravity_calculate_acceleration(integrator_is_whfasthelio);
        
        // Calculate spin variation and non-gravity accelerations.
        self.move_to_star_center();
        let evolution = true;
        let dspin_dt = false;
        let accelerations = true;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations);
        self.move_to_center_of_mass();

        // ---------------------------------------------------------------------
        // A 'DKD'-like integrator will do the 'KD' part.
        self.to_helio_posvel();
        self.helio_interaction_step(time_step);
        self.helio_jump_step(time_step);
        self.helio_kepler_steps(half_time_step);
        self.to_inertial_posvel();
        self.current_time += self.half_time_step;
        // ---------------------------------------------------------------------
        self.move_to_star_center();
    }

    fn iterate_spin_with_midpoint_method(&mut self) {
        let time_step = self.time_step;
        let half_time_step = self.half_time_step;

        // Midpoint method: https://en.wikipedia.org/wiki/Midpoint_method
        let evolution = false; // Don't evolve particles at this point or it messes up the conservation of angular momentum
        let dspin_dt = true;
        let accelerations = false;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations);

        self.save_current_spin();
        self.spin_step(half_time_step);
        let evolution = false; // Don't evolve particles at this point or it messes up the conservation of angular momentum
        let dspin_dt = true;
        let accelerations = false;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations);
        self.restore_last_spin();
        self.spin_step(time_step);
    }

    fn save_current_spin(&mut self) {
        for (particle, spin) in self.universe.particles[..self.universe.n_particles].iter().zip(self.last_spin[..self.universe.n_particles].iter_mut()) {
            spin.x = particle.spin.x;
            spin.y = particle.spin.y;
            spin.z = particle.spin.z;
        }
    }

    fn restore_last_spin(&mut self) {
        for (particle, spin) in self.universe.particles[..self.universe.n_particles].iter_mut().zip(self.last_spin[..self.universe.n_particles].iter()) {
            particle.spin.x = spin.x;
            particle.spin.y = spin.y;
            particle.spin.z = spin.z;
        }
    }
 
    fn spin_step(&mut self, _dt: f64) {
        for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            if particle.moment_of_inertia_ratio != 1. {
                particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + _dt * particle.dspin_dt.x;
                particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + _dt * particle.dspin_dt.y;
                particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + _dt * particle.dspin_dt.z;
            } else {
                particle.spin.x = particle.spin.x + _dt * particle.dspin_dt.x;
                particle.spin.y = particle.spin.y + _dt * particle.dspin_dt.y;
                particle.spin.z = particle.spin.z + _dt * particle.dspin_dt.z;
            }
            if particle.wind_factor != 0. {
                // TODO: Verify wind factor
                particle.spin.x += _dt * particle.wind_factor * particle.spin.x;
                particle.spin.y += _dt * particle.wind_factor * particle.spin.y;
                particle.spin.z += _dt * particle.wind_factor * particle.spin.z;
            }
        }
    }

    /***************************** 
     * Operators                 */
    fn helio_jump_step(&mut self, _dt: f64){
        let m0 = self.universe.particles[0].mass;
        let mut px = 0.;
        let mut py = 0.;
        let mut pz = 0.;
        for particle_heliocentric in self.universe_heliocentric.particles[1..self.universe_heliocentric.n_particles].iter() {
            px += particle_heliocentric.mass* particle_heliocentric.velocity.x;
            py += particle_heliocentric.mass* particle_heliocentric.velocity.y;
            pz += particle_heliocentric.mass* particle_heliocentric.velocity.z;
        }
        for particle_heliocentric in self.universe_heliocentric.particles[1..self.universe_heliocentric.n_particles].iter_mut() {
            particle_heliocentric.position.x += _dt * px/m0;
            particle_heliocentric.position.y += _dt * py/m0;
            particle_heliocentric.position.z += _dt * pz/m0;
        }
        self.is_synchronized = false;
    }

    fn helio_interaction_step(&mut self, _dt: f64){
        for (particle_heliocentric, particle) in self.universe_heliocentric.particles[1..self.universe_heliocentric.n_particles].iter_mut().zip(self.universe.particles[1..self.universe.n_particles].iter()) {
            particle_heliocentric.velocity.x += _dt*particle.acceleration.x;
            particle_heliocentric.velocity.y += _dt*particle.acceleration.y;
            particle_heliocentric.velocity.z += _dt*particle.acceleration.z;
        }
        self.is_synchronized = false;
    }

    fn helio_kepler_steps(&mut self, time_step: f64){
        let star_mass_g = self.universe.particles[0].mass_g;

        for i in 1..self.universe.n_particles {
            self.kepler_individual_step(i, star_mass_g, time_step);
        }
        let star_heliocentric = &mut self.universe_heliocentric.particles[0];
        star_heliocentric.position.x += time_step*star_heliocentric.velocity.x;
        star_heliocentric.position.y += time_step*star_heliocentric.velocity.y;
        star_heliocentric.position.z += time_step*star_heliocentric.velocity.z;
        self.is_synchronized = false;
    }

    //***************************** 
    // Keplerian motion           
    fn kepler_individual_step(&mut self, i: usize, mass_g: f64, _dt: f64){
        // Save a copy of the original position and velocities
        let p1_position;
        let p1_velocity;
        {
            let p1 = &self.universe_heliocentric.particles[i];
            p1_position = p1.position.clone();
            p1_velocity = p1.velocity.clone();
        }

        let r0 = (p1_position.x.powi(2) + p1_position.y.powi(2) + p1_position.z.powi(2)).sqrt();
        let r0i = 1./r0;
        let v2 = p1_velocity.x.powi(2) + p1_velocity.y.powi(2) + p1_velocity.z.powi(2);
        let beta = 2.*mass_g*r0i - v2;
        let eta0 = p1_position.x*p1_velocity.x + p1_position.y*p1_velocity.y + p1_position.z*p1_velocity.z;
        let zeta0 = mass_g - beta*r0;
        let mut x;
        let mut gs;
        let mut invperiod = 0.; // only used for beta>0.
        let x_per_period;
        
        if beta > 0. {
            //// Elliptic orbit
            let sqrt_beta = beta.sqrt();
            invperiod = sqrt_beta*beta / (2.*PI*mass_g);
            x_per_period = 2.*PI / sqrt_beta;
            if _dt.abs()*invperiod > 1. && self.timestep_warning == 0 {
                // Ignoring const qualifiers. This warning should not have any effect on
                // other parts of the code, nor is it vital to show it.
                println!("[WARNING {} UTC] WHFast convergence issue. Timestep is larger than at least one orbital period.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
                self.timestep_warning += 1;
            }
            //x = _dt*invperiod*x_per_period; // first order guess 
            let dtr0i = _dt*r0i;
            //x = dtr0i; // first order guess
            x = dtr0i * (1. - dtr0i*eta0*0.5*r0i); // second order guess
            //x = dtr0i *(1.- 0.5*dtr0i*r0i*(eta0-dtr0i*(eta0*eta0*r0i-1./3.*zeta0))); // third order guess
            //x = _dt*beta/mass_g + eta0/mass_g*(0.85*sqrt(1.+zeta0*zeta0/beta/eta0/eta0) - 1.);  // Dan's version 
        } else {
            //// Hyperbolic orbit
            x = 0.; // Initial guess 
            x_per_period = std::f64::NAN; // only used for beta>0. nan triggers Newton's method for beta<0
        }


        let mut converged = 0;
        let mut old_x = x;

        //// Do one Newton step
        gs = WHFastHelio::stiefel_gs3(beta, x);
        let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
        let mut ri = 1./(r0 + eta0_gs1_zeta0_gs2);
        x  = ri*(x*eta0_gs1_zeta0_gs2-eta0*gs[2]-zeta0*gs[3]+_dt);

        // Choose solver depending on estimated step size
        // Note, for hyperbolic orbits this uses Newton's method.
        if !x_per_period.is_nan() && (x-old_x).abs() > 0.01*x_per_period {
            // Quartic solver
            // Linear initial guess
            x = beta*_dt/mass_g;
            let mut prev_x = [0.; WHFAST_NMAX_QUART+1];
            'outer: for n_lag in 1..WHFAST_NMAX_QUART {
                gs = WHFastHelio::stiefel_gs3(beta, x);
                let f = r0*x + eta0*gs[2] + zeta0*gs[3] - _dt;
                let fp = r0 + eta0*gs[1] + zeta0*gs[2];
                let fpp = eta0*gs[0] + zeta0*gs[1];
                let denom = fp + (16.*fp*fp - 20.*f*fpp).abs().sqrt();
                x = (x*denom - 5.*f)/denom;
                'inner: for i in 1..n_lag {
                    if x == prev_x[i] {
                        // Converged. Exit.
                        //n_lag = WHFAST_NMAX_QUART;
                        converged = 1;
                        break 'outer;
                    }
                }
                prev_x[n_lag] = x;
            }
            let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
            ri = 1./(r0 + eta0_gs1_zeta0_gs2);
        } else {
            // Newton's method
            let mut old_x2;
            for _ in 1..WHFAST_NMAX_NEWT {
                old_x2 = old_x;
                old_x = x;
                gs = WHFastHelio::stiefel_gs3(beta, x);
                let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
                ri = 1./(r0 + eta0_gs1_zeta0_gs2);
                x  = ri*(x*eta0_gs1_zeta0_gs2-eta0*gs[2]-zeta0*gs[3]+_dt);
                
                if x==old_x || x==old_x2 {
                    // Converged. Exit.
                    converged = 1;
                    break; 
                }
            }
        }
            
        // If solver did not work, fallback to bisection 
        if converged == 0 { 
            let mut x_min;
            let mut x_max;
            if beta > 0. {
                //Elliptic
                x_min = x_per_period * (_dt*invperiod).floor();
                x_max = x_min + x_per_period;
            } else {
                //Hyperbolic
                let h2 = r0*r0*v2-eta0*eta0;
                let q = h2/mass_g/(1.+(1.-h2*beta/(mass_g*mass_g)).sqrt());
                let vq = (h2).sqrt()/q;
                x_min = 1./(vq+r0/_dt);
                x_max = _dt/q;
            }
            x = (x_max + x_min)/2.;
            loop {
                gs = WHFastHelio::stiefel_gs3(beta, x);
                let s   = r0*x + eta0*gs[2] + zeta0*gs[3]-_dt;
                if s >= 0. {
                    x_max = x;
                } else {
                    x_min = x;
                }
                x = (x_max + x_min)/2.;

                if (x_max-x_min).abs()/x_max <= 1e-15 {
                    break;
                }
            }
            let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
            ri = 1./(r0 + eta0_gs1_zeta0_gs2);
        }
        if ri.is_nan() {
            // Exception for (almost) straight line motion in hyperbolic case
            ri = 0.;
            gs[1] = 0.;
            gs[2] = 0.;
            gs[3] = 0.;
        }

        // Note: These are not the traditional f and g functions.
        let f = -mass_g*gs[2]*r0i;
        let g = _dt - mass_g*gs[3];
        let fd = -mass_g*gs[1]*r0i*ri; 
        let gd = -mass_g*gs[2]*ri; 
            
        let p_j = &mut self.universe_heliocentric.particles[i];
        p_j.position.x += f*p1_position.x + g*p1_velocity.x;
        p_j.position.y += f*p1_position.y + g*p1_velocity.y;
        p_j.position.z += f*p1_position.z + g*p1_velocity.z;
        
        // WARNING: p1_position used below should not be modified by the previous block (keep p_j
        // independent)
        p_j.velocity.x += fd*p1_position.x + gd*p1_velocity.x;
        p_j.velocity.y += fd*p1_position.y + gd*p1_velocity.y;
        p_j.velocity.z += fd*p1_position.z + gd*p1_velocity.z;
    }

    fn stiefel_gs3(beta: f64, x: f64) -> [f64; 6] {
        let x2 = x.powi(2);
        let mut gs = WHFastHelio::stumpff_cs3(beta*x2);
        gs[1] *= x; 
        gs[2] *= x2; 
        gs[3] *= x2*x;
        gs
    }

    fn stumpff_cs3(z: f64) -> [f64; 6] {
        // Fast inverse factorial lookup table
        let invfactorial: [f64; 35] = [1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.];
        let mut z = z;
        let mut n = 0;
        while z.abs() > 0.1 {
            z = z/4.;
            n += 1;
        }
        let mut cs = [0.; 6];
        let nmax = 13;
        let mut c_odd = invfactorial[nmax];
        let mut c_even = invfactorial[nmax-1];

        let mut np = nmax-2;
        while np >= 3 {
            c_odd  = invfactorial[np]    - z *c_odd;
            c_even = invfactorial[np-1]  - z *c_even;
            np -= 2;
        }
        cs[3] = c_odd;
        cs[2] = c_even;
        cs[1] = invfactorial[1]  - z *c_odd;
        cs[0] = invfactorial[0]  - z *c_even;
        while n > 0 {
            cs[3] = (cs[2]+cs[0]*cs[3])*0.25;
            cs[2] = cs[1]*cs[1]*0.5;
            cs[1] = cs[0]*cs[1];
            cs[0] = 2.*cs[0]*cs[0]-1.;
            n -= 1;
        }
        cs
    }

    //***************************** 
    // Coordinate transformations 
    //***************************** 
    fn to_helio_posvel(&mut self){
        {
            let star_heliocentric = &mut self.universe_heliocentric.particles[0];
            star_heliocentric.position.x  = 0.;
            star_heliocentric.position.y  = 0.;
            star_heliocentric.position.z  = 0.;
            star_heliocentric.velocity.x = 0.;
            star_heliocentric.velocity.y = 0.;
            star_heliocentric.velocity.z = 0.;
            star_heliocentric.mass  = 0.;
            for particle in self.universe.particles[0..self.universe.n_particles].iter() {
                star_heliocentric.position.x  += particle.position.x *particle.mass;
                star_heliocentric.position.y  += particle.position.y *particle.mass;
                star_heliocentric.position.z  += particle.position.z *particle.mass;
                star_heliocentric.velocity.x += particle.velocity.x*particle.mass;
                star_heliocentric.velocity.y += particle.velocity.y*particle.mass;
                star_heliocentric.velocity.z += particle.velocity.z*particle.mass;
                star_heliocentric.mass  += particle.mass;
            }
            star_heliocentric.position.x  /= star_heliocentric.mass;
            star_heliocentric.position.y  /= star_heliocentric.mass;
            star_heliocentric.position.z  /= star_heliocentric.mass;
            star_heliocentric.velocity.x /= star_heliocentric.mass;
            star_heliocentric.velocity.y /= star_heliocentric.mass;
            star_heliocentric.velocity.z /= star_heliocentric.mass;
        }

        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            if let Some((star_heliocentric, particles_heliocentric)) = self.universe_heliocentric.particles[..self.universe_heliocentric.n_particles].split_first_mut() {
                for (particle_heliocentric, particle) in particles_heliocentric.iter_mut().zip(particles.iter()) {
                    particle_heliocentric.position.x  = particle.position.x  - star.position.x ; 
                    particle_heliocentric.position.y  = particle.position.y  - star.position.y ;
                    particle_heliocentric.position.z  = particle.position.z  - star.position.z ;
                    particle_heliocentric.velocity.x = particle.velocity.x - star_heliocentric.velocity.x;
                    particle_heliocentric.velocity.y = particle.velocity.y - star_heliocentric.velocity.y;
                    particle_heliocentric.velocity.z = particle.velocity.z - star_heliocentric.velocity.z;
                    particle_heliocentric.mass  = particle.mass;
                }
            }
        }
        self.is_synchronized = true;
    }

    fn to_inertial_pos(&mut self) {
        let mtot = self.universe_heliocentric.particles[0].mass;
        {
            let mut new_star_position = self.universe_heliocentric.particles[0].position; // Copy
            for (particle_heliocentric, particle) in self.universe_heliocentric.particles[1..self.universe_heliocentric.n_particles].iter().zip(self.universe.particles[1..self.universe.n_particles].iter()) {
                new_star_position.x  -= particle_heliocentric.position.x*particle.mass/mtot;
                new_star_position.y  -= particle_heliocentric.position.y*particle.mass/mtot;
                new_star_position.z  -= particle_heliocentric.position.z*particle.mass/mtot;
            }
            let star = &mut self.universe.particles[0];
            star.position.x  = new_star_position.x;
            star.position.y  = new_star_position.y;
            star.position.z  = new_star_position.z;
        }
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for (particle_heliocentric, particle) in self.universe_heliocentric.particles[1..self.universe_heliocentric.n_particles].iter().zip(particles.iter_mut()) {
                particle.position.x = particle_heliocentric.position.x+star.position.x;
                particle.position.y = particle_heliocentric.position.y+star.position.y;
                particle.position.z = particle_heliocentric.position.z+star.position.z;
            }
        }

    }
    fn to_inertial_posvel(&mut self) {
        self.to_inertial_pos();
        let mtot = self.universe_heliocentric.particles[0].mass;
        let m0 = self.universe.particles[0].mass;
        let factor = mtot/m0;
        if let Some((star_heliocentric, particles_heliocentric)) = self.universe_heliocentric.particles[..self.universe_heliocentric.n_particles].split_first_mut() {
            for (particle_heliocentric, particle) in particles_heliocentric.iter().zip(self.universe.particles[1..self.universe.n_particles].iter_mut()) {
                particle.velocity.x = particle_heliocentric.velocity.x+star_heliocentric.velocity.x;
                particle.velocity.y = particle_heliocentric.velocity.y+star_heliocentric.velocity.y;
                particle.velocity.z = particle_heliocentric.velocity.z+star_heliocentric.velocity.z;
            }
        }

        let mut new_star_velocity = self.universe_heliocentric.particles[0].velocity.clone();
        new_star_velocity.x = new_star_velocity.x*factor;
        new_star_velocity.y = new_star_velocity.y*factor;
        new_star_velocity.z = new_star_velocity.z*factor;

        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                let factor = particle.mass/m0;
                new_star_velocity.x -= particle.velocity.x*factor;
                new_star_velocity.y -= particle.velocity.y*factor;
                new_star_velocity.z -= particle.velocity.z*factor;
            }

            star.velocity.x = new_star_velocity.x;
            star.velocity.y = new_star_velocity.y;
            star.velocity.z = new_star_velocity.z;
        }
        self.is_synchronized = true;
    }

    fn get_center_of_mass_of_pair(&self, center_of_mass_position: &mut Axes, center_of_mass_velocity: &mut Axes, center_of_mass_acceleration: &mut Axes, center_of_mass_mass: f64, particle: &Particle) -> f64 {
        center_of_mass_position.x   = center_of_mass_position.x*center_of_mass_mass + particle.position.x*particle.mass;
        center_of_mass_position.y   = center_of_mass_position.y*center_of_mass_mass + particle.position.y*particle.mass;
        center_of_mass_position.z   = center_of_mass_position.z*center_of_mass_mass + particle.position.z*particle.mass;
        center_of_mass_velocity.x  = center_of_mass_velocity.x*center_of_mass_mass + particle.velocity.x*particle.mass;
        center_of_mass_velocity.y  = center_of_mass_velocity.y*center_of_mass_mass + particle.velocity.y*particle.mass;
        center_of_mass_velocity.z  = center_of_mass_velocity.z*center_of_mass_mass + particle.velocity.z*particle.mass;
        center_of_mass_acceleration.x  = center_of_mass_acceleration.x*center_of_mass_mass + particle.acceleration.x*particle.mass;
        center_of_mass_acceleration.y  = center_of_mass_acceleration.y*center_of_mass_mass + particle.acceleration.y*particle.mass;
        center_of_mass_acceleration.z  = center_of_mass_acceleration.z*center_of_mass_mass + particle.acceleration.z*particle.mass;
        
        let new_center_of_mass_mass = center_of_mass_mass + particle.mass;
        if new_center_of_mass_mass > 0. {
            center_of_mass_position.x  /= new_center_of_mass_mass;
            center_of_mass_position.y  /= new_center_of_mass_mass;
            center_of_mass_position.z  /= new_center_of_mass_mass;
            center_of_mass_velocity.x /= new_center_of_mass_mass;
            center_of_mass_velocity.y /= new_center_of_mass_mass;
            center_of_mass_velocity.z /= new_center_of_mass_mass;
            center_of_mass_acceleration.x /= new_center_of_mass_mass;
            center_of_mass_acceleration.y /= new_center_of_mass_mass;
            center_of_mass_acceleration.z /= new_center_of_mass_mass;
        }
        new_center_of_mass_mass
    }


    fn move_to_center_of_mass(&mut self) {
        if self.set_to_center_of_mass {
            return;
        }
        if !self.is_synchronized {
            panic!("[PANIC {} UTC] Non synchronized particles cannot be moved to center of mass.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap())
        }

        // Compute center of mass
        let mut center_of_mass_position = Axes{x:0., y:0., z:0.};
        let mut center_of_mass_velocity = Axes{x:0., y:0., z:0.};
        let mut center_of_mass_acceleration = Axes{x:0., y:0., z:0.};
        let mut center_of_mass_mass = 0.;

        for particle in self.universe.particles[..self.universe.n_particles].iter() {
            center_of_mass_mass = self.get_center_of_mass_of_pair(&mut center_of_mass_position, 
                                                                    &mut center_of_mass_velocity, 
                                                                    &mut center_of_mass_acceleration,
                                                                    center_of_mass_mass,
                                                                    &particle);
        }

        // Move
        for particle in self.universe.particles[0..self.universe.n_particles].iter_mut() {
            particle.position.x  -= center_of_mass_position.x;
            particle.position.y  -= center_of_mass_position.y;
            particle.position.z  -= center_of_mass_position.z;
            particle.velocity.x -= center_of_mass_velocity.x;
            particle.velocity.y -= center_of_mass_velocity.y;
            particle.velocity.z -= center_of_mass_velocity.z;
            particle.acceleration.x -= center_of_mass_acceleration.x;
            particle.acceleration.y -= center_of_mass_acceleration.y;
            particle.acceleration.z -= center_of_mass_acceleration.z;
        }
        self.set_to_center_of_mass = true;
    }

    fn move_to_star_center(&mut self) {
        if !self.set_to_center_of_mass {
            return;
        }
        if !self.is_synchronized {
            panic!("[PANIC {} UTC] Non synchronized particles cannot be moved to star center.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap())
        }

        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.position.x  -= star.position.x;
                particle.position.y  -= star.position.y;
                particle.position.z  -= star.position.z;
                particle.velocity.x -= star.velocity.x;
                particle.velocity.y -= star.velocity.y;
                particle.velocity.z -= star.velocity.z;
                particle.acceleration.x -= star.acceleration.x;
                particle.acceleration.y -= star.acceleration.y;
                particle.acceleration.z -= star.acceleration.z;
            }
            star.position.x = 0.;
            star.position.y = 0.;
            star.position.z = 0.;
            star.velocity.x = 0.;
            star.velocity.y = 0.;
            star.velocity.z = 0.;
            star.acceleration.x = 0.;
            star.acceleration.y = 0.;
            star.acceleration.z = 0.;
        }


        self.set_to_center_of_mass = false;
    }

}
