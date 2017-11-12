use std;

pub const MAX_PARTICLES : usize = 10; // The optimal value matches the real number of bodies (it will generate smaller snapshots), but a greater number will work too.
pub const MAX_DISTANCE : f64 = 100.; // AU
pub const MAX_DISTANCE_2 : f64 = MAX_DISTANCE*MAX_DISTANCE; // AU
//
pub const MIN_ORBITAL_PERIOD_TIME_STEP_RATIO : f64 = 10.; // The orbital period should be N times greater than the time step to correctly integrate an orbit

//// Constants for IAS15 integrator (to be ignored for others)
pub const INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT : bool = true;	// Turn this off to safe some time if the force is not velocity dependent (i.e. radiation forces, tides depend on vel.).
pub const INTEGRATOR_EPSILON_GLOBAL : bool = true;  // if true: estimate the fractional error by max(acceleration_error)/max(acceleration), where max is take over all particles.
                                                    // if false: estimate the fractional error by max(acceleration_error/acceleration).
pub const INTEGRATOR_EPSILON : f64 = 1e-6;          // Precision parameter (default: 1e-9)
    							                    // If it is zero, then a constant timestep is used. 
pub const INTEGRATOR_MIN_DT : f64 = 0.;             // Minimum timestep used as a floor when adaptive timestepping is enabled (default: 0.).
pub const SAFETY_FACTOR : f64 = 0.25;               // Maximum increase/deacrease of consecutve timesteps (default: 0.25).

///// Constants for WHFastHelio
pub const WHFAST_NMAX_QUART : usize = 64;               // Maximum number of iterations for quartic solver
pub const WHFAST_NMAX_NEWT : usize = 32;               // Maximum number of iterations for Newton's method


//// Physical constants
// These are the same units as used by the mercury6 code.
const K  : f64 = 0.01720209895;    // Gaussian constant 
pub const K2 : f64 = K*K; 
pub const G  : f64 = K2;  // Gravitational constant in Mercury units
pub const PI : f64 = std::f64::consts::PI;
pub const TWO_PI : f64 = std::f64::consts::PI * 2.;
pub const DEG2RAD : f64 = std::f64::consts::PI / 180.; // conversion factor from degrees to radians
pub const HOUR: f64 = 3600.; // s
pub const DAY: f64 = 24.*HOUR; // s

////////////////////////////////////////////////////////////////////////////////
// The IAU 2009 system of astronomical constants: the report of the IAU working group on numerical standards for Fundamental Astronomy
// https://en.wikipedia.org/wiki/Astronomical_constant
pub const M_SUN : f64 =  1.9818e30; // Kg
pub const M2EARTH : f64 = 332946.050895; // ratio sun/earth mass
//pub const M_EARTH : f64 = M_SUN/M2EARTH; // kg
pub const M_EARTH : f64 = 1./M2EARTH; // M_SUN
pub const G_SI : f64 = 6.67428e-11;  // m^3.kg^-1.s^-2 (S.I. units)
pub const AU : f64 = 1.49597870700e11; // m
pub const M2AU  : f64 = 1./AU; // 6.684587122268445e-12 AU (1 meter in AU)
pub const SPEED_OF_LIGHT : f64 = (2.99792458e8/AU)*DAY; // 173.14463267424034 AU/day
pub const SPEED_OF_LIGHT_2 : f64 = SPEED_OF_LIGHT*SPEED_OF_LIGHT;

// Resolution B3 on recommended nominal conversion constants for selected solar and planetary properties
// http://adsabs.harvard.edu/abs/2015arXiv151007674M
pub const R_SUN : f64 = 6.957e8/AU;  // AU
pub const R_EARTH : f64 = 6.3781e6/AU; // AU
////////////////////////////////////////////////////////////////////////////////

// Solar system
pub const SUN_DYN_FREQ : f64 = K2/(R_SUN*R_SUN*R_SUN); // Needed for MathisSolarLike

