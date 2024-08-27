extern crate posidonius;
extern crate time;
extern crate clap;
use clap::{Arg, ArgAction, Command};
use clap::value_parser;
use std::path::Path;
use std::time::{Duration, Instant};
use time::{OffsetDateTime, format_description};

fn main() {
    let timer = Instant::now();

    let matches = Command::new("Posidonius")
                            .version("2017.01.27")
                            .author("Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/")
                            .about("N-Body simulator for planetary system affected by tidal effects. Based on Mercury-T from Emeline Bolmont.")
                            .subcommand(Command::new("start")
                                    .about("Start a simulation")
                                    .arg(Arg::new("start_case_filename")
                                        .required(true)
                                        .index(1)
                                        .help("JSON case description"))
                                    .arg(Arg::new("snapshot_filename")
                                        .required(true)
                                        .index(2)
                                        .help("Recovery snapshot filename")) 
                                    .arg(Arg::new("historic_snapshot_filename")
                                        .required(true)
                                        .index(3)
                                        .help("Historic snapshot filename")) 
                                    .arg(Arg::new("limit")
                                        .short('l')
                                        .long("limit")
                                        .value_name("seconds")
                                        .value_parser(value_parser!(u64))
                                        .help("Maximum execution time limit (default: no limit)"))
                                    .arg(Arg::new("silent")
                                        .short('s')
                                        .long("silent")
                                        .action(ArgAction::SetTrue)
                                        .help("Only print INFO/WARNING/ERROR messages"))
                                     )
                            .subcommand(Command::new("resume")
                                    .about("Resume a simulation")
                                    .arg(Arg::new("resume_case_filename")
                                        .required(true)
                                        .help("Recovery snapshot filename"))
                                    .arg(Arg::new("historic_snapshot_filename")
                                        .required(true)
                                        .help("Historic snapshot filename")) 
                                    .arg(Arg::new("limit")
                                        .short('l')
                                        .long("limit")
                                        .value_name("seconds")
                                        .value_parser(value_parser!(u64))
                                        .help("Maximum execution time limit (default: no limit)"))
                                    .arg(Arg::new("silent")
                                        .short('s')
                                        .long("silent")
                                        .action(ArgAction::SetTrue)
                                        .help("Only print INFO/WARNING/ERROR messages"))
                                    .arg(Arg::new("change_historic_snapshot_period")
                                        .long("historic-snapshot-period")
                                        .value_name("days")
                                        .value_parser(value_parser!(f64))
                                        .help("Set new historic snapshots period in days."))
                                    .arg(Arg::new("change_recovery_snapshot_period")
                                        .long("recovery-snapshot-period")
                                        .value_name("days")
                                        .value_parser(value_parser!(f64))
                                        .help("Set new recovery snapshots period in days."))
                                    .arg(Arg::new("change_time_limit")
                                        .long("time-limit")
                                        .value_name("days")
                                        .value_parser(value_parser!(u64))
                                        .help("Set new simulation time limit in days."))
                                    )
                            .subcommand_required(true)
                            .arg_required_else_help(true)
                          .get_matches();

    let first_universe_integrator_snapshot_filename;
    let universe_integrator_snapshot_filename;
    let universe_history_filename;
    let silent_mode;
    let resume;
    let new_historic_snapshot_period;
    let new_recovery_snapshot_period;
    let new_time_limit;
    let execution_time_limit;

    match matches.subcommand() {
        Some(("start", start_matches)) => {
            first_universe_integrator_snapshot_filename = start_matches.get_one::<String>("start_case_filename").unwrap();
            universe_integrator_snapshot_filename = start_matches.get_one::<String>("snapshot_filename").unwrap();
            universe_history_filename = start_matches.get_one::<String>("historic_snapshot_filename").unwrap();
            silent_mode = start_matches.get_flag("silent");
            resume = false;
            new_historic_snapshot_period = -1.0;
            new_recovery_snapshot_period = -1.0;
            new_time_limit = -1.0;
            execution_time_limit = Duration::from_secs(start_matches.get_one::<u64>("limit").copied().unwrap_or(0));
        },
        Some(("resume", resume_matches)) => {
            universe_integrator_snapshot_filename = resume_matches.get_one::<String>("resume_case_filename").unwrap();
            first_universe_integrator_snapshot_filename = &universe_integrator_snapshot_filename;
            universe_history_filename = resume_matches.get_one::<String>("historic_snapshot_filename").unwrap();
            silent_mode = resume_matches.get_flag("silent");
            resume = true;

            new_historic_snapshot_period = resume_matches.get_one::<f64>("change_historic_snapshot_period").copied().unwrap_or(-1.);
            new_recovery_snapshot_period = resume_matches.get_one::<f64>("change_recovery_snapshot_period").copied().unwrap_or(-1.);
            new_time_limit = resume_matches.get_one::<f64>("change_time_limit").copied().unwrap_or(-1.);
            execution_time_limit = Duration::from_secs(resume_matches.get_one::<u64>("limit").copied().unwrap_or(0));
        },
        _ => unreachable!(),
    }

    let universe_integrator_snapshot_path = Path::new(&universe_integrator_snapshot_filename);

    let first_universe_integrator_snapshot_path = Path::new(&first_universe_integrator_snapshot_filename);
    
    // Start/Resume from snapshot
    let mut boxed_universe_integrator : Box<dyn posidonius::Integrator> = match posidonius::output::restore_snapshot(&first_universe_integrator_snapshot_path) {
        Ok(restored_case) => { restored_case },
        Err(_) => { 
            if resume {
                panic!("[PANIC {} UTC] It was not possible to resume the simulation", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            } else {
                panic!("[PANIC {} UTC] It was not possible to start the simulation", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            }
        },
    };

    boxed_universe_integrator.set_snapshot_periods(new_historic_snapshot_period, new_recovery_snapshot_period);
    boxed_universe_integrator.set_time_limit(new_time_limit);

    // Create/recover historic snapshot
    let expected_n_bytes = (boxed_universe_integrator.get_n_historic_snapshots() as u64) 
                                * posidonius::output::n_bytes_per_particle_in_historic_snapshot()
                                * (boxed_universe_integrator.get_n_particles() as u64);
    let universe_history_path = Path::new(&universe_history_filename);

    if !resume && universe_integrator_snapshot_path.exists() {
        panic!("[PANIC {} UTC] File '{}' already exists.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), universe_integrator_snapshot_filename);
    } else if !resume && universe_history_path.exists() {
        panic!("[PANIC {} UTC] File '{}' already exists.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), universe_history_filename);
    }

    let mut universe_history_writer = posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);

    // Simulate
    let instant = Instant::now();
    let enabled_execution_time_limit = match execution_time_limit.as_secs() {
        0 => false,
        _ => true,
    };
    loop {
        match boxed_universe_integrator.iterate(&mut universe_history_writer, silent_mode) {
            Ok(recovery_snapshot_time_trigger) => {
                if enabled_execution_time_limit {
                    let elapsed = instant.elapsed();
                    if elapsed >= execution_time_limit {
                        // Save a universe snapshot so that we can resume later on
                        boxed_universe_integrator.write_recovery_snapshot(&universe_integrator_snapshot_path, &mut universe_history_writer);
                        println!("[WARNING {} UTC] Reached execution time limit before simulation completion", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap()); 
                        break;
                    }
                } else if recovery_snapshot_time_trigger {
                    // Save a universe snapshot so that we can resume in case of failure
                    boxed_universe_integrator.write_recovery_snapshot(&universe_integrator_snapshot_path, &mut universe_history_writer);
                }
            },
            Err(e) => { println!("[INFO {} UTC] {} '{}'.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), e, first_universe_integrator_snapshot_filename); break; }
        };
    }

    let duration = timer.elapsed(); // duration
    if !resume {
        println!("[INFO {} UTC] Execution time: {} seconds", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), duration.as_secs_f64());
    } else {
        println!("[INFO {} UTC] Execution time since last resume: {} seconds", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), duration.as_secs_f64());
    }
}

