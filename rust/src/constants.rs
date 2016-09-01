use std;
use super::integrator::IntegratorType;

//// Simulation
pub const N_PARTICLES: usize = 2;
pub const TIME_STEP: f64 = 0.08; // in days
pub const TIME_LIMIT: f64 = TIME_STEP * 365.25 * 10.0e8;
//pub const TIME_LIMIT: f64 = TIME_STEP * 365.25 * 10.0e6;
//pub const TIME_LIMIT: f64 = TIME_STEP*2.;
pub const PRINT_EVERY_N_DAYS: f64 = 100.*365.25;
pub const PRINT_EVERY_N_ITERATIONS: u32 = 100_000_000;


//// Integrator
pub const INTEGRATOR: IntegratorType = IntegratorType::LeapFrog;
//pub const INTEGRATOR: IntegratorType = IntegratorType::Ias15;

//// Constants for IAS15 integrator (to be ignored for others)
pub const INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT : bool = true;	// Turn this off to safe some time if the force is not velocity dependent (i.e. radiation forces, tides depend on vel.).
pub const INTEGRATOR_EPSILON_GLOBAL : bool = true;  // if true: estimate the fractional error by max(acceleration_error)/max(acceleration), where max is take over all particles.
                                                    // if false: estimate the fractional error by max(acceleration_error/acceleration).
pub const INTEGRATOR_EPSILON : f64 = 1e-9;          // Precision parameter 
    							                    // If it is zero, then a constant timestep is used. 
pub const INTEGRATOR_MIN_DT : f64 = 0.;             // Minimum timestep used as a floor when adaptive timestepping is enabled.
pub const SAFETY_FACTOR : f64 = 0.25;               // Maximum increase/deacrease of consecutve timesteps.


//// Physical constants
// These are the same units as used by the mercury6 code.
const K  : f64 = 0.01720209895;    // Gaussian constant 
pub const K2 : f64 = K*K; 
pub const G  : f64 = K2;  // Gravitational constant
pub const PI : f64 = std::f64::consts::PI;
pub const TWO_PI : f64 = std::f64::consts::PI * 2.;
pub const DEG2RAD : f64 = std::f64::consts::PI / 180.; // conversion factor from degrees to radians

// Solar system
pub const R_SUN : f64 = 4.67920694e-3;
pub const R_EARTH : f64 = 4.25874677e-5; // AU
pub const M2EARTH : f64 = (1.9891e6/5.9794); // Factor for mass-radius relation (valid only for earth type planets)
