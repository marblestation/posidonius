extern crate posidonius;
extern crate time;
use posidonius::Integrator;
use std::path::Path;


fn main() {
    let t1 = time::precise_time_s();

    ////////////////////////////////////////////////////////////////////////////
    // 1.- Indicate output filenames
    ////////////////////////////////////////////////////////////////////////////
    let universe_history_filename = String::from("target/universe_history.bin");
    let universe_integrator_snapshot_filename = String::from("target/universe_integrator_snapshot.bin");
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // 2.- Select case to be run:
    ////////////////////////////////////////////////////////////////////////////
    //let run_case = posidonius::cases::bolmont_et_al_2015::case3();
    //let run_case = posidonius::cases::bolmont_et_al_2015::case4();
    //let run_case = posidonius::cases::bolmont_et_al_2015::case7();
    let run_case = posidonius::cases::main_example();
    //let run_case = posidonius::cases::example_with_helpers();
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // >>> Nothing else to do from here
    ////////////////////////////////////////////////////////////////////////////

    // Resume from snapshot if it exists
    let universe_integrator_snapshot_path = Path::new(&universe_integrator_snapshot_filename);
    let mut universe_integrator = match posidonius::WHFastHelio::restore_snapshot(&universe_integrator_snapshot_path) {
        Ok(restored_case) => {
            if run_case.hash != restored_case.hash {
                panic!("Initial conditions from the restored snapshot file does not match the initial conditions of the selected case to be run.");
            }
            restored_case
        },
        Err(_) => { 
            //// If there is no snapshot, create a case to init the universe to be simulated:
            run_case
        },
    };

    // Create/recover historic snapshot
    let expected_n_bytes = (universe_integrator.n_historic_snapshots as u64) 
                                * posidonius::output::n_bytes_per_particle_in_historic_snapshot()
                                * (universe_integrator.universe.n_particles as u64);
    let universe_history_path = Path::new(&universe_history_filename);
    let mut universe_history_writer = posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);

    // Simulate
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

