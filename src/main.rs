extern crate posidonius;
extern crate time;
#[macro_use]
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
                                    .arg(Arg::with_name("silent")
                                        .short("s")
                                        .long("silent")
                                        .multiple(false)
                                        .help("Only print INFO/WARNING/ERROR messages"))
                                     )
                            .subcommand(SubCommand::with_name("resume")
                                    .about("Resume a simulation")
                                    .arg(Arg::with_name("resume_case_filename")
                                        .required(true)
                                        .help("Recovery snapshot filename"))
                                    .arg(Arg::with_name("historic_snapshot_filename")
                                        .required(true)
                                        .help("Historic snapshot filename")) 
                                    .arg(Arg::with_name("silent")
                                        .short("s")
                                        .long("silent")
                                        .multiple(false)
                                        .help("Only print INFO/WARNING/ERROR messages"))
                                    .arg(Arg::with_name("change_historic_snapshot_period")
                                        .long("historic-snapshot-period")
                                        .value_name("days")
                                        .help("Set new historic snapshots period in days."))
                                    .arg(Arg::with_name("change_recovery_snapshot_period")
                                        .long("recovery-snapshot-period")
                                        .value_name("days")
                                        .help("Set new recovery snapshots period in days."))
                                    )
                            .setting(AppSettings::SubcommandRequiredElseHelp)
                          .get_matches();

    let first_universe_integrator_snapshot_filename;
    let universe_integrator_snapshot_filename;
    let universe_history_filename;
    let silent_mode;
    let resume;
    let new_historic_snapshot_period;
    let new_recovery_snapshot_period;

    match matches.subcommand() {
        ("start", Some(start_matches)) =>{
            first_universe_integrator_snapshot_filename = start_matches.value_of("start_case_filename").unwrap();
            universe_integrator_snapshot_filename = start_matches.value_of("snapshot_filename").unwrap();
            universe_history_filename = start_matches.value_of("historic_snapshot_filename").unwrap();
            silent_mode = start_matches.is_present("silent");
            resume = false;
            new_historic_snapshot_period = -1.0;
            new_recovery_snapshot_period = -1.0;
        },
        ("resume", Some(resume_matches)) =>{
            universe_integrator_snapshot_filename = resume_matches.value_of("resume_case_filename").unwrap();
            first_universe_integrator_snapshot_filename = &universe_integrator_snapshot_filename;
            universe_history_filename = resume_matches.value_of("historic_snapshot_filename").unwrap();
            silent_mode = resume_matches.is_present("silent");
            resume = true;
            new_historic_snapshot_period = value_t!(resume_matches.value_of("change_historic_snapshot_period"), f64).unwrap_or(-1.);
            new_recovery_snapshot_period = value_t!(resume_matches.value_of("change_recovery_snapshot_period"), f64).unwrap_or(-1.);
        },
        ("", None)   => unreachable!(),
        _            => unreachable!(),
    }

    let universe_integrator_snapshot_path = Path::new(&universe_integrator_snapshot_filename);

    let first_universe_integrator_snapshot_path = Path::new(&first_universe_integrator_snapshot_filename);
    
    // Start/Resume from snapshot
    let mut boxed_universe_integrator : Box<Integrator> = match posidonius::output::restore_snapshot(&first_universe_integrator_snapshot_path) {
        Ok(restored_case) => { restored_case },
        Err(_) => { 
            if resume {
                panic!("[PANIC {} UTC] It was not possible to resume the simulation", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
            } else {
                panic!("[PANIC {} UTC] It was not possible to start the simulation", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
            }
        },
    };

    boxed_universe_integrator.set_snapshot_periods(new_historic_snapshot_period, new_recovery_snapshot_period);

    // Create/recover historic snapshot
    let expected_n_bytes = (boxed_universe_integrator.get_n_historic_snapshots() as u64) 
                                * posidonius::output::n_bytes_per_particle_in_historic_snapshot()
                                * (boxed_universe_integrator.get_n_particles() as u64);
    let universe_history_path = Path::new(&universe_history_filename);

    if !resume && universe_integrator_snapshot_path.exists() {
        panic!("[PANIC {} UTC] File '{}' already exists.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_filename);
    } else if !resume && universe_history_path.exists() {
        panic!("[PANIC {} UTC] File '{}' already exists.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_history_filename);
    }

    let mut universe_history_writer = posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);

    // Simulate
    loop {
        match boxed_universe_integrator.iterate(&mut universe_history_writer, silent_mode) {
            Ok(recovery_snapshot_time_trigger) => {
                if recovery_snapshot_time_trigger {
                    // Save a universe snapshot so that we can resume in case of failure
                    boxed_universe_integrator.write_recovery_snapshot(&universe_integrator_snapshot_path, &mut universe_history_writer);
                }
            },
            Err(e) => { println!("[WARNING {} UTC] {} '{}'.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), e, first_universe_integrator_snapshot_filename); break; }
        };
    }

    let t2 = time::precise_time_s();
    let d = t2 - t1;
    if !resume {
        println!("[INFO {} UTC] Execution time: {} seconds", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), d);
    } else {
        println!("[INFO {} UTC] Execution time since last resume: {} seconds", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), d);
    }
}

