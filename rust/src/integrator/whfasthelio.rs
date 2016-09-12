use std;
use std::io::{Write, BufWriter};
use super::Integrator;
use super::super::constants::{N_PARTICLES, PRINT_EVERY_N_DAYS, INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT, PI, WHFAST_NMAX_QUART, WHFAST_NMAX_NEWT};
use super::super::particle::Particles;
use super::super::particle::Particle;
use super::super::particle::Axes;
use super::output::{write_bin_snapshot};

/// WHFastHelio (symplectic integrator) to be used always in safe mode (always sync because it is required by tides)
/// and without correction (i.e. 2nd order integrator, comparable to mercury symplectic part of the hybrid integrator)
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

pub struct WHFastHelio {
    time_step: f64,
    half_time_step: f64,
    time_limit: f64,
    particles: Particles,
    current_time: f64,
    current_iteration: usize,
    last_print_time: f64,
    /**
     * @brief This variable turns on/off different symplectic correctors for WHFastHelio. Same as for WHFast.
     * @details 
     * - 0 (default): turns off all correctors
     * - 3: uses third order (two-stage) corrector 
     * - 5: uses fifth order (four-stage) corrector 
     * - 7: uses seventh order (six-stage) corrector 
     * - 11: uses eleventh order (ten-stage) corrector 
     */
    corrector: usize, 
    /** 
     * @brief Setting this flag to one will recalculate heliocentric coordinates from the particle structure in the next timestep. 
     * @details After the timestep, the flag gets set back to 0. 
     * If you want to change particles after every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    recalculate_heliocentric_this_timestep: bool,
    /**
     * @brief If this flag is set (the default), WHFastHelio will recalculate heliocentric
     * coordinates and synchronize every timestep, to avoid problems with outputs or 
     * particle modifications between timesteps. 
     * @details Setting it to 0 will result in a speedup, but care
     * must be taken to synchronize and recalculate heliocentric coordinates when needed.
     */
    safe_mode: bool,
    /**
     * @brief Heliocentric coordinates
     * @details This array contains the heliocentric coordinates of all particles.
     * It is automatically filled and updated by WHfastDemocratic.
     * Access this array with caution.
     */
    particles_heliocentric: Particles,

    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    ///< Flag to determine if current particle structure is synchronized
    is_synchronized: bool, 
    ///< Counter of heliocentric synchronization errors
    recalculate_heliocentric_but_not_synchronized_warning: usize,
    timestep_warning: usize ,
    set_to_center_of_mass: bool
}

impl Integrator for WHFastHelio {

    fn new(time_step: f64, time_limit: f64, particles: Particles) -> WHFastHelio {
        WHFastHelio {
                    time_step:time_step,
                    half_time_step:0.5*time_step,
                    time_limit:time_limit,
                    last_print_time:-1.,
                    particles:particles,
                    current_time:0.,
                    current_iteration:0,
                    // WHFastHelio specifics:
                    corrector: 0, 
                    recalculate_heliocentric_this_timestep: true,
                    safe_mode: true,
                    particles_heliocentric: particles,
                    is_synchronized: true,
                    recalculate_heliocentric_but_not_synchronized_warning: 0,
                    timestep_warning: 0,
                    set_to_center_of_mass: false,
                    }
    }

    fn iterate<T: Write>(&mut self, output_bin: &mut BufWriter<T>) -> Result<(), String> {
        // Output
        let add_header = self.last_print_time < 0.;
        let time_triger = self.last_print_time + PRINT_EVERY_N_DAYS <= self.current_time;
        if add_header || time_triger {
            if self.set_to_center_of_mass {
                self.move_to_star_center();
            }
            write_bin_snapshot(output_bin, &self.particles, self.current_time, self.time_step);
            let current_time_years = self.current_time/365.25;
            print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
            let _ = std::io::stdout().flush();

            if add_header || time_triger {
                self.last_print_time = self.current_time;
            } 
        }

        //// Calculate non-gravity accelerations.
        if self.set_to_center_of_mass {
            self.move_to_star_center();
        }
        let only_dspin_dt = true;
        self.particles.calculate_additional_forces(only_dspin_dt);
        
        // A 'DKD'-like integrator will do the first 'D' part.
        if !self.set_to_center_of_mass {
            self.move_to_center_of_mass();
        }
        self.integrator_part1();

        // Calculate accelerations.
        self.particles.gravity_calculate_acceleration();
        
        //// Calculate non-gravity accelerations.
        if self.set_to_center_of_mass {
            self.move_to_star_center();
        }
        let only_dspin_dt = false;
        self.particles.calculate_additional_forces(only_dspin_dt);

        // A 'DKD'-like integrator will do the 'KD' part.
        if !self.set_to_center_of_mass {
            self.move_to_center_of_mass();
        }
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


impl WHFastHelio {
    // WHFastHelio integrator
    #[allow(dead_code)]
    fn integrator_part1(&mut self) {
        if self.safe_mode || self.recalculate_heliocentric_this_timestep {
            if !self.is_synchronized {
                self.synchronize();
                if self.recalculate_heliocentric_but_not_synchronized_warning == 0 {
                    println!("Recalculating heliocentric coordinates but pos/vel were not synchronized before.");
                    self.recalculate_heliocentric_but_not_synchronized_warning += 1;
                }

            }
            self.recalculate_heliocentric_this_timestep = false;
            self.to_helio_posvel();
        }

        if self.is_synchronized {
            // First half DRIFT step
            if self.corrector != 0 {
                //reb_whfast_apply_corrector(1., reb_whfasthelio_corrector_Z);
                panic!("Correctors not implemented yet!");
            }
            let half_time_step = true;
            self.kepler_steps(half_time_step);
        }else{
            // Combined DRIFT step
            let half_time_step = false;
            self.kepler_steps(half_time_step);
        }

        // For force calculation:
        if INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT {
            self.to_inertial_posvel();
        } else {
            self.to_inertial_pos();
        }

        for particle in self.particles.particles.iter_mut() {
            particle.spin.x += self.half_time_step * particle.dspin_dt.x;
            particle.spin.y += self.half_time_step * particle.dspin_dt.y;
            particle.spin.z += self.half_time_step * particle.dspin_dt.z;
        }
    
        self.current_time += self.half_time_step;
    }

    #[allow(dead_code)]
    fn integrator_part2(&mut self) {
        let time_step = self.time_step;
        self.interaction_step(time_step);
        self.jump_step(time_step);
        
        self.is_synchronized = false;
        if self.safe_mode {
            self.synchronize();
        }

        for particle in self.particles.particles.iter_mut() {
            particle.spin.x += self.half_time_step * particle.dspin_dt.x;
            particle.spin.y += self.half_time_step * particle.dspin_dt.y;
            particle.spin.z += self.half_time_step * particle.dspin_dt.z;
        }

        self.current_time += self.half_time_step;
    }

    fn synchronize(&mut self) {
        if !self.is_synchronized {
            let half_time_step = true;
            self.kepler_steps(half_time_step);
            if self.corrector != 0 {
                //reb_whfast_apply_corrector(-1., reb_whfasthelio_corrector_Z);
                panic!("Correctors not implemented yet!");
            }
            self.to_inertial_posvel();
        }
        self.is_synchronized = true;
    }

    /***************************** 
     * Operators                 */
    fn jump_step(&mut self, _dt: f64){
        let m0 = self.particles.particles[0].mass;
        let mut px = 0.;
        let mut py = 0.;
        let mut pz = 0.;
        for particle_heliocentric in self.particles_heliocentric.particles[1..].iter() {
            px += particle_heliocentric.mass* particle_heliocentric.velocity.x;
            py += particle_heliocentric.mass* particle_heliocentric.velocity.y;
            pz += particle_heliocentric.mass* particle_heliocentric.velocity.z;
        }
        for particle_heliocentric in self.particles_heliocentric.particles[1..].iter_mut() {
            particle_heliocentric.position.x += _dt * px/m0;
            particle_heliocentric.position.y += _dt * py/m0;
            particle_heliocentric.position.z += _dt * pz/m0;
        }
    }

    fn interaction_step(&mut self, _dt: f64){
        for (particle_heliocentric, particle) in self.particles_heliocentric.particles[1..].iter_mut().zip(self.particles.particles[1..].iter()) {
            particle_heliocentric.velocity.x += _dt*particle.acceleration.x;
            particle_heliocentric.velocity.y += _dt*particle.acceleration.y;
            particle_heliocentric.velocity.z += _dt*particle.acceleration.z;
        }
    }

    fn kepler_steps(&mut self, half_time_step: bool){
        let star_mass_g = self.particles.particles[0].mass_g;

        for i in 1..N_PARTICLES {
            self.kepler_individual_step(i, star_mass_g, half_time_step);
        }
        let time_step = match half_time_step {
            true => self.half_time_step,
            false => self.time_step
        };
        let star_heliocentric = &mut self.particles_heliocentric.particles[0];
        star_heliocentric.position.x += time_step*star_heliocentric.velocity.x;
        star_heliocentric.position.y += time_step*star_heliocentric.velocity.y;
        star_heliocentric.position.z += time_step*star_heliocentric.velocity.z;
    }

    //***************************** 
    // Keplerian motion           
    fn kepler_individual_step(&mut self, i: usize, mass_g: f64, half_time_step: bool){
        let _dt = match half_time_step {
            true => self.half_time_step,
            false => self.time_step
        };
        let p1 = self.particles_heliocentric.particles[i];

        let r0 = (p1.position.x.powi(2) + p1.position.y.powi(2) + p1.position.z.powi(2)).sqrt();
        let r0i = 1./r0;
        let v2 = p1.velocity.x.powi(2) + p1.velocity.y.powi(2) + p1.velocity.z.powi(2);
        let beta = 2.*mass_g*r0i - v2;
        let eta0 = p1.position.x*p1.velocity.x + p1.position.y*p1.velocity.y + p1.position.z*p1.velocity.z;
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
                println!("WHFast convergence issue. Timestep is larger than at least one orbital period.");
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
        gs = self.stiefel_gs3(beta, x);
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
                gs = self.stiefel_gs3(beta, x);
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
                gs = self.stiefel_gs3(beta, x);
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
                gs = self.stiefel_gs3(beta, x);
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
            
        let p_j = &mut self.particles_heliocentric.particles[i];
        p_j.position.x += f*p1.position.x + g*p1.velocity.x;
        p_j.position.y += f*p1.position.y + g*p1.velocity.y;
        p_j.position.z += f*p1.position.z + g*p1.velocity.z;
        
        // WARNING: p1.position used below should not be modified by the previous block (keep p_j
        // independent)
        p_j.velocity.x += fd*p1.position.x + gd*p1.velocity.x;
        p_j.velocity.y += fd*p1.position.y + gd*p1.velocity.y;
        p_j.velocity.z += fd*p1.position.z + gd*p1.velocity.z;
    }

    fn stiefel_gs3(&mut self, beta: f64, x: f64) -> [f64; 6] {
        let x2 = x.powi(2);
        let mut gs = self.stumpff_cs3(beta*x2);
        gs[1] *= x; 
        gs[2] *= x2; 
        gs[3] *= x2*x;
        gs
    }

    fn stumpff_cs3(&mut self, z: f64) -> [f64; 6] {
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
            let star_heliocentric = &mut self.particles_heliocentric.particles[0];
            star_heliocentric.position.x  = 0.;
            star_heliocentric.position.y  = 0.;
            star_heliocentric.position.z  = 0.;
            star_heliocentric.velocity.x = 0.;
            star_heliocentric.velocity.y = 0.;
            star_heliocentric.velocity.z = 0.;
            star_heliocentric.mass  = 0.;
            for particle in self.particles.particles[0..].iter() {
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
        let local_copy_star = self.particles.particles[0];
        let local_copy_star_heliocentric = self.particles_heliocentric.particles[0];
        for (particle_heliocentric, particle) in self.particles_heliocentric.particles[1..].iter_mut().zip(self.particles.particles[1..].iter()) {
            particle_heliocentric.position.x  = particle.position.x  - local_copy_star.position.x ; 
            particle_heliocentric.position.y  = particle.position.y  - local_copy_star.position.y ;
            particle_heliocentric.position.z  = particle.position.z  - local_copy_star.position.z ;
            particle_heliocentric.velocity.x = particle.velocity.x - local_copy_star_heliocentric.velocity.x;
            particle_heliocentric.velocity.y = particle.velocity.y - local_copy_star_heliocentric.velocity.y;
            particle_heliocentric.velocity.z = particle.velocity.z - local_copy_star_heliocentric.velocity.z;
            particle_heliocentric.mass  = particle.mass;
        }
    }

    fn to_inertial_pos(&mut self) {
        let mtot = self.particles_heliocentric.particles[0].mass;
        {
            let mut new_star_position = self.particles_heliocentric.particles[0].position; // Copy
            for (particle_heliocentric, particle) in self.particles_heliocentric.particles[1..].iter().zip(self.particles.particles[1..].iter()) {
                new_star_position.x  -= particle_heliocentric.position.x*particle.mass/mtot;
                new_star_position.y  -= particle_heliocentric.position.y*particle.mass/mtot;
                new_star_position.z  -= particle_heliocentric.position.z*particle.mass/mtot;
            }
            let star = &mut self.particles.particles[0];
            star.position.x  = new_star_position.x;
            star.position.y  = new_star_position.y;
            star.position.z  = new_star_position.z;
        }
        let local_copy_star = self.particles.particles[0];
        for (particle_heliocentric, particle) in self.particles_heliocentric.particles[1..].iter().zip(self.particles.particles[1..].iter_mut()) {
            particle.position.x = particle_heliocentric.position.x+local_copy_star.position.x;
            particle.position.y = particle_heliocentric.position.y+local_copy_star.position.y;
            particle.position.z = particle_heliocentric.position.z+local_copy_star.position.z;
        }
    }
    fn to_inertial_posvel(&mut self) {
        self.to_inertial_pos();
        let mtot = self.particles_heliocentric.particles[0].mass;
        let m0 = self.particles.particles[0].mass;
        let factor = mtot/m0;
        let local_copy_star_heliocentric = self.particles_heliocentric.particles[0];
        for (particle_heliocentric, particle) in self.particles_heliocentric.particles[1..].iter().zip(self.particles.particles[1..].iter_mut()) {
            particle.velocity.x = particle_heliocentric.velocity.x+local_copy_star_heliocentric.velocity.x;
            particle.velocity.y = particle_heliocentric.velocity.y+local_copy_star_heliocentric.velocity.y;
            particle.velocity.z = particle_heliocentric.velocity.z+local_copy_star_heliocentric.velocity.z;
        }

        let mut new_star_velocity = self.particles_heliocentric.particles[0].velocity; // Copy
        new_star_velocity.x = new_star_velocity.x*factor;
        new_star_velocity.y = new_star_velocity.y*factor;
        new_star_velocity.z = new_star_velocity.z*factor;

        for particle in self.particles.particles[1..].iter_mut() {
            let factor = particle.mass/m0;
            new_star_velocity.x -= particle.velocity.x*factor;
            new_star_velocity.y -= particle.velocity.y*factor;
            new_star_velocity.z -= particle.velocity.z*factor;
        }

        let star = &mut self.particles.particles[0];
        star.velocity.x = new_star_velocity.x;
        star.velocity.y = new_star_velocity.y;
        star.velocity.z = new_star_velocity.z;
    }

    fn get_center_of_mass_of_pair(&self, center_of_mass: &mut Particle, particle: &Particle) {
        center_of_mass.position.x   = center_of_mass.position.x*center_of_mass.mass + particle.position.x*particle.mass;
        center_of_mass.position.y   = center_of_mass.position.y*center_of_mass.mass + particle.position.y*particle.mass;
        center_of_mass.position.z   = center_of_mass.position.z*center_of_mass.mass + particle.position.z*particle.mass;
        center_of_mass.velocity.x  = center_of_mass.velocity.x*center_of_mass.mass + particle.velocity.x*particle.mass;
        center_of_mass.velocity.y  = center_of_mass.velocity.y*center_of_mass.mass + particle.velocity.y*particle.mass;
        center_of_mass.velocity.z  = center_of_mass.velocity.z*center_of_mass.mass + particle.velocity.z*particle.mass;
        center_of_mass.acceleration.x  = center_of_mass.acceleration.x*center_of_mass.mass + particle.acceleration.x*particle.mass;
        center_of_mass.acceleration.y  = center_of_mass.acceleration.y*center_of_mass.mass + particle.acceleration.y*particle.mass;
        center_of_mass.acceleration.z  = center_of_mass.acceleration.z*center_of_mass.mass + particle.acceleration.z*particle.mass;
        
        center_of_mass.mass  += particle.mass;
        if center_of_mass.mass > 0. {
            center_of_mass.position.x  /= center_of_mass.mass;
            center_of_mass.position.y  /= center_of_mass.mass;
            center_of_mass.position.z  /= center_of_mass.mass;
            center_of_mass.velocity.x /= center_of_mass.mass;
            center_of_mass.velocity.y /= center_of_mass.mass;
            center_of_mass.velocity.z /= center_of_mass.mass;
            center_of_mass.acceleration.x /= center_of_mass.mass;
            center_of_mass.acceleration.y /= center_of_mass.mass;
            center_of_mass.acceleration.z /= center_of_mass.mass;
        }
    }


    fn move_to_center_of_mass(&mut self) {
        if self.set_to_center_of_mass {
            return;
        }
        if !self.is_synchronized {
            panic!("Non synchronized particles cannot be moved to center of mass.")
        }

        // Compute center of mass
        let mass = 0.;
        let radius = 0.;
        let dissipation_factor = 0.;
        let radius_of_gyration_2 = 0.;
        let love_number = 0.;
        let position = Axes{x:0., y:0., z:0. };
        let velocity = Axes{x:0., y:0., z:0. };
        let acceleration = Axes{x:0., y:0., z:0. };
        let spin = Axes{x:0., y:0., z:0. };
        let mut center_of_mass = Particle::new(mass, radius, dissipation_factor, radius_of_gyration_2, love_number,
                                                position, velocity, acceleration, spin);

        for particle in self.particles.particles.iter() {
            self.get_center_of_mass_of_pair(&mut center_of_mass, &particle);
        }

        // Move
        for particle in self.particles.particles[0..].iter_mut() {
            particle.position.x  -= center_of_mass.position.x;
            particle.position.y  -= center_of_mass.position.y;
            particle.position.z  -= center_of_mass.position.z;
            particle.velocity.x -= center_of_mass.velocity.x;
            particle.velocity.y -= center_of_mass.velocity.y;
            particle.velocity.z -= center_of_mass.velocity.z;
            //particle.acceleration.x -= center_of_mass.acceleration.x; // Always zero
            //particle.acceleration.y -= center_of_mass.acceleration.y;
            //particle.acceleration.z -= center_of_mass.acceleration.z;
        }
        self.to_helio_posvel();
        self.set_to_center_of_mass = true;
    }

    fn move_to_star_center(&mut self) {
        if !self.set_to_center_of_mass {
            return;
        }
        if !self.is_synchronized {
            panic!("Non synchronized particles cannot be moved to star center.")
        }
        let center = self.particles.particles[0];
        for particle in self.particles.particles[0..].iter_mut() {
            particle.position.x  -= center.position.x;
            particle.position.y  -= center.position.y;
            particle.position.z  -= center.position.z;
            particle.velocity.x -= center.velocity.x;
            particle.velocity.y -= center.velocity.y;
            particle.velocity.z -= center.velocity.z;
            //particle.acceleration.x -= center.acceleration.x; // Should be always zero
            //particle.acceleration.y -= center.acceleration.y;
            //particle.acceleration.z -= center.acceleration.z;
        }
        self.to_helio_posvel();
        self.set_to_center_of_mass = false;
    }

}
