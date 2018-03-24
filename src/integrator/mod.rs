mod leapfrog;
mod ias15;
mod whfast;
pub mod output;

pub use self::leapfrog::*;
pub use self::ias15::*;
pub use self::whfast::*;

use std::io::{BufWriter};
use std::fs::File;
use std::path::Path;


pub trait Integrator {
    fn get_n_historic_snapshots(&self) -> usize;
    fn get_n_particles(&self) -> usize;
    fn get_current_time(&self) -> f64;
    fn set_snapshot_periods(&mut self, historic_snapshot_period: f64, recovery_snapshot_period: f64);
    fn initialize_physical_values(&mut self);
    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String>;
    fn write_recovery_snapshot(&mut self, snapshot_path: &Path, universe_history_writer: &mut BufWriter<File>);
}


