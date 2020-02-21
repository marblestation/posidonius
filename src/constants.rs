use std;

pub const MAX_PARTICLES : usize = 10; // The optimal value matches the real number of bodies (it will generate smaller snapshots), but a greater number will work too.
pub const MAX_DISTANCE : f64 = 100.; // AU
pub const MAX_DISTANCE_2 : f64 = MAX_DISTANCE*MAX_DISTANCE; // AU (use a negative value to disable)
//
pub const MIN_ORBITAL_PERIOD_TIME_STEP_RATIO : f64 = 5.; // The orbital period should be N times greater than the time step to correctly integrate an orbit (use a negative value to disable)
pub const COMPENSATED_GRAVITY : bool = false;

//// Constants for IAS15 integrator (to be ignored for others)
pub const INTEGRATOR_FORCE_IS_VELOCITYDEPENDENT : bool = true;	// Turn this off to safe some time if the force is not velocity dependent (i.e. radiation forces, tides depend on vel.).
pub const INTEGRATOR_EPSILON_GLOBAL : bool = true;  // if true: estimate the fractional error by max(acceleration_error)/max(acceleration), where max is take over all particles.
                                                    // if false: estimate the fractional error by max(acceleration_error/acceleration).
pub const INTEGRATOR_EPSILON : f64 = 1e-6;          // Precision parameter (default: 1e-9)
    							                    // If it is zero, then a constant timestep is used. 
pub const INTEGRATOR_MIN_DT : f64 = 0.;             // Minimum timestep used as a floor when adaptive timestepping is enabled (default: 0.).
pub const INTEGRATOR_MAX_DT : f64 = 0.;             // Maximum timestep used as a ceiling when adaptive timestepping is enabled (default: 0. => disabled).
pub const SAFETY_FACTOR : f64 = 0.25;               // Maximum increase/deacrease of consecutve timesteps (default: 0.25).

///// Constants for WHFast
pub const WHFAST_NMAX_QUART : usize = 64;               // Maximum number of iterations for quartic solver
pub const WHFAST_NMAX_NEWT : usize = 32;               // Maximum number of iterations for Newton's method

pub const DBL_EPSILON: f64 = 2.2204460492503131e-16; // https://en.wikipedia.org/wiki/Machine_epsilon
pub const DBL_EPSILON_2 : f64 = DBL_EPSILON*DBL_EPSILON;


//// Physical constants
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

////////////////////////////////////////////////////////////////////////////////
// Constants in the units used in the Mercury code
////////////////////////////////////////////////////////////////////////////////
const G_MERCURY : f64 = G_SI/(AU*AU*AU) * M_SUN * (DAY*DAY); // Gravitational constant with Mercury units
pub const G : f64 = G_MERCURY;
pub const K2 : f64 = G_MERCURY;
//const K  : f64 = G_MERCURY.sqrt(); // Gaussian constant (in Mercury is 0.01720209895 but it does not follow from the recommended IAU constants)


// Solar system
pub const SUN_DYN_FREQ : f64 = K2/(R_SUN*R_SUN*R_SUN); // Needed for MathisSolarLike

// Boltzmann constant and mass of an hydrogen atom from CODATA 2017
// Table 3 of Newell et al. 2018 (https://iopscience.iop.org/article/10.1088/1681-7575/aa950a/pdf) 
pub const BOLTZMANN_CONSTANT_SI : f64 = 1.380649e-23; // J.K-1
pub const MASS_HYDROGEN_ATOM_SI : f64 = 1.672162841e-27; // kg
pub const BOLTZMANN_CONSTANT : f64 = BOLTZMANN_CONSTANT_SI/(M_SUN*AU*AU) * DAY*DAY; // M_SUN.AU^2.DAY^-2.K-1
pub const MASS_HYDROGEN_ATOM : f64 = MASS_HYDROGEN_ATOM_SI/M_SUN; // kg
