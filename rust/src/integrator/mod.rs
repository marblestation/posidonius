extern crate rusqlite;

mod leapfrog;
mod ias15;
mod whfasthelio;
mod output;

pub use self::leapfrog::*;
pub use self::ias15::*;
pub use self::whfasthelio::*;

use super::particle::Particles;
use std::io::{Write, BufWriter};

#[derive(Debug,Copy, Clone)]
pub enum IntegratorType {
    LeapFrog,
    Ias15,
    WHFastHelio,
}

pub trait Integrator {
    fn new(time_step: f64, time_limit: f64, particles: Particles) -> Self;
    fn iterate<T: Write>(&mut self, output_text: &mut BufWriter<T>, output_bin: &mut BufWriter<T>, output_db: &rusqlite::Connection) -> Result<(), String>;
}
