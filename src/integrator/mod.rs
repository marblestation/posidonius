mod leapfrog;
mod ias15;
mod whfast;
pub mod output;

pub use self::leapfrog::*;
pub use self::ias15::*;
pub use self::whfast::*;

use std::io::{BufWriter};
use std::fs::File;

pub trait Integrator {
    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String>;
    fn prepare_for_recovery_snapshot(&mut self, universe_history_writer: &mut BufWriter<File>);
}
