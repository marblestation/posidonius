extern crate rusqlite;

mod leapfrog;
mod ias15;
mod output;

pub use self::leapfrog::*;
pub use self::ias15::*;

use super::particle::Particles;
use std::io::{Write, BufWriter};

pub enum IntegratorType {
    LeapFrog,
    Ias15,
}

pub trait Integrator {
    fn new(time_step: f64, time_limit: f64, particles: Particles) -> Self;
    fn iterate<T: Write>(&mut self, output_text: &mut BufWriter<T>, output_bin: &mut BufWriter<T>, output_db: &rusqlite::Connection) -> Result<(), String>;
}
