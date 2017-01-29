extern crate posidonius;
extern crate time;
extern crate clap;
use clap::{Arg, App, SubCommand, AppSettings};
use posidonius::Integrator;
use std::path::Path;


fn main() {
    let t1 = time::precise_time_s();

    let matches = App::new("Posidonius")
                            .version("2017.01.27")
                            .author("Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/")
                            .about("N-Body simulator for planetary system affected by tidal effects. Based on Mercury-T from Emeline Bolmont.")
                            .subcommand(SubCommand::with_name("start")
                                      .about("Start a simulation")
                                      .arg(Arg::with_name("start_case_filename")
                                          .required(true)
                                          .index(1)
                                          .help("JSON case description"))
                                      .arg(Arg::with_name("snapshot_filename")
                                          .required(true)
                                          .index(2)
                                          .help("Recovery snapshot filename")) 
                                      .arg(Arg::with_name("historic_snapshot_filename")
                                          .required(true)
                                          .index(3)
                                          .help("Historic snapshot filename")) 
                                     )
                            .subcommand(SubCommand::with_name("resume")
                                      .about("Resume a simulation")
                                      .arg(Arg::with_name("resume_case_filename")
                                          .required(true)
                                          .help("Recovery snapshot filename"))
                                      .arg(Arg::with_name("historic_snapshot_filename")
                                          .required(true)
                                          .help("Historic snapshot filename")) 
                                      )
                            .subcommand(SubCommand::with_name("default")
                                      .about("Start the default simulation (hard-coded in Posidonius)")
                                      .arg(Arg::with_name("snapshot_filename")
                                          .required(true)
                                          .help("Recovery snapshot filename"))
                                      .arg(Arg::with_name("historic_snapshot_filename")
                                          .required(true)
                                          .help("Historic snapshot filename")) 
                                      )
                            .setting(AppSettings::SubcommandRequiredElseHelp)
                          .get_matches();

    let first_universe_integrator_snapshot_filename;
    let universe_integrator_snapshot_filename;
    let universe_history_filename;
    let resume;
    let default;

    match matches.subcommand() {
        ("start", Some(start_matches)) =>{
            first_universe_integrator_snapshot_filename = start_matches.value_of("start_case_filename").unwrap();
            universe_integrator_snapshot_filename = start_matches.value_of("snapshot_filename").unwrap();
            universe_history_filename = start_matches.value_of("historic_snapshot_filename").unwrap();
            resume = false;
            default = false;
        },
        ("resume", Some(resume_matches)) =>{
            universe_integrator_snapshot_filename = resume_matches.value_of("resume_case_filename").unwrap();
            first_universe_integrator_snapshot_filename = &universe_integrator_snapshot_filename;
            universe_history_filename = resume_matches.value_of("historic_snapshot_filename").unwrap();
            resume = true;
            default = false;
        },
        ("default", Some(default_matches)) =>{
            universe_integrator_snapshot_filename = default_matches.value_of("snapshot_filename").unwrap();
            first_universe_integrator_snapshot_filename = &universe_integrator_snapshot_filename;
            universe_history_filename = default_matches.value_of("historic_snapshot_filename").unwrap();
            resume = false;
            default = true;
        },
        ("", None)   => unreachable!(),
        _            => unreachable!(),
    }

    let universe_integrator_snapshot_path = Path::new(&universe_integrator_snapshot_filename);

    //let default_case = posidonius::cases::bolmont_et_al_2015::case3();
    //let default_case = posidonius::cases::bolmont_et_al_2015::case4();
    //let default_case = posidonius::cases::bolmont_et_al_2015::case7();
    let default_case = posidonius::cases::main_example();
    //let default_case = posidonius::cases::example_with_helpers();

    let first_universe_integrator_snapshot_path = Path::new(&first_universe_integrator_snapshot_filename);
    let mut universe_integrator = match default {
        true => {
            // Start default case
            default_case
        },
        false => {
            // Start/Resume from snapshot
            match posidonius::WHFastHelio::restore_snapshot(&first_universe_integrator_snapshot_path) {
                Ok(restored_case) => { restored_case },
                Err(_) => { 
                    if resume {
                        panic!("It was not possible to resume the simulation");
                    } else {
                        panic!("It was not possible to start the simulation");
                    }
                },
            }
        }
    };

    // Create/recover historic snapshot
    let expected_n_bytes = (universe_integrator.n_historic_snapshots as u64) 
                                * posidonius::output::n_bytes_per_particle_in_historic_snapshot()
                                * (universe_integrator.universe.n_particles as u64);
    let universe_history_path = Path::new(&universe_history_filename);

    if !resume && universe_integrator_snapshot_path.exists() {
        panic!("File '{}' already exists.", universe_integrator_snapshot_filename);
    } else if !resume && universe_history_path.exists() {
        panic!("File '{}' already exists.", universe_history_filename);
    }

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

