mod leapfrog;
mod ias15;
mod whfasthelio;
mod output;

pub use self::leapfrog::*;
pub use self::ias15::*;
pub use self::whfasthelio::*;

use super::particles::Universe;
use std::io::{Write, BufWriter};

#[derive(Debug,Copy, Clone)]
pub enum IntegratorType {
    LeapFrog,
    Ias15,
    WHFastHelio,
}

pub trait Integrator {
    fn new(time_step: f64, time_limit: f64, universe: Universe) -> Self;
    fn iterate<T: Write>(&mut self, output_bin: &mut BufWriter<T>) -> Result<(), String>;
}
