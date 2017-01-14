extern crate posidonius;
extern crate time;
use posidonius::Integrator;
use std::path::Path;
use std::fs::{OpenOptions};
use std::io::{BufWriter};
use std::fs;
use std::fs::File;

fn get_universe_history_writer(universe_history_path: &Path) -> BufWriter<File> {
    if universe_history_path.exists() {
        // Backup previous file
        let new_extension = format!("{0}.bin", time::now().strftime("%Y%m%dT%H%M%S").unwrap());
        fs::rename(universe_history_path, universe_history_path.with_extension(new_extension)).unwrap();
    }
    // We create file options to write
    let mut options_bin = OpenOptions::new();
    options_bin.create(true).truncate(true).write(true);

    let universe_history_file = match options_bin.open(&universe_history_path) {
        Ok(f) => f,
        Err(e) => panic!("file error: {}", e),
    };
    ////////////////////////////////////////////////////////////////////////////
    let universe_history_writer = BufWriter::new(universe_history_file);
    ////////////////////////////////////////////////////////////////////////////
    universe_history_writer
}

fn main() {
    let t1 = time::precise_time_s();
    let universe_history_filename = String::from("target/universe_history.bin");
    let universe_integrator_snapshot_filename = String::from("target/universe_integrator_snapshot.bin");

    let universe_history_path = Path::new(&universe_history_filename);
    let mut universe_history_writer = get_universe_history_writer(universe_history_path);

    // Resume from snapshot if it exists
    let universe_integrator_snapshot_path = Path::new(&universe_integrator_snapshot_filename);
    let mut universe_integrator = match posidonius::WHFastHelio::restore_snapshot(&universe_integrator_snapshot_path) {
        Ok(v) => v,
        Err(_) => { 
            //// If there is no snapshot, create a case to init the universe to be simulated:
            //posidonius::cases::bolmont_et_al_2015::case4()
            //posidonius::cases::bolmont_et_al_2015::case3()
            //posidonius::cases::bolmont_et_al_2015::case7()
            posidonius::cases::main_example() 
            //posidonius::cases::example_with_helpers()
        },
    };


    loop {
        match universe_integrator.iterate(&mut universe_history_writer) {
            Ok(recovery_snapshot_time_trigger) => {
                if recovery_snapshot_time_trigger {
                    // Save a universe snapshot so that we can resume in case of failure
                    universe_integrator.prepare_for_recovery_snapshot(&mut universe_history_writer);
                    posidonius::output::write_recovery_snapshot(&universe_integrator_snapshot_path, &universe_integrator);
                }
            },
            Err(e) => { println!("{}", e); break; }
        };
    }

    let t2 = time::precise_time_s();
    let d = t2 - t1;
    println!("Execution time: {} seconds", d);
}

