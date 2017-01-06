extern crate posidonius;
extern crate time;

//use std::num;
//use std::io::timer;
//use std::time::Duration;
use std::path::Path;
use std::fs::{OpenOptions};
use std::io::{BufWriter};

use posidonius::Integrator;



fn main() {
    let t1 = time::precise_time_s();
    //timer::sleep(Duration::milliseconds(1000));

    let universe = posidonius::cases::case3();


    ////////////////////////////////////////////////////////////////////////////
    let path_bin = Path::new("target/output.bin");
    // We create file options to write
    let mut options_bin = OpenOptions::new();
    options_bin.create(true).truncate(true).write(true);

    let output_bin_file = match options_bin.open(&path_bin) {
        Ok(f) => f,
        Err(e) => panic!("file error: {}", e),
    };
    ////////////////////////////////////////////////////////////////////////////
    let mut output_bin = BufWriter::new(&output_bin_file);
    ////////////////////////////////////////////////////////////////////////////

    // TODO: Improve (dynamic dispatching?)
    let mut leapfrog = posidonius::LeapFrog::new(posidonius::constants::TIME_STEP, posidonius::constants::TIME_LIMIT, universe.clone());
    let mut ias15 = posidonius::Ias15::new(posidonius::constants::TIME_STEP, posidonius::constants::TIME_LIMIT, universe.clone());
    let mut whfasthelio = posidonius::WHFastHelio::new(posidonius::constants::TIME_STEP, posidonius::constants::TIME_LIMIT, universe.clone());

    loop {
        match posidonius::constants::INTEGRATOR {
            posidonius::IntegratorType::LeapFrog => {
                match leapfrog.iterate(&mut output_bin) {
                    Ok(_) => {},
                    Err(e) => { println!("{}", e); break; }
                };
            },
            posidonius::IntegratorType::Ias15 => {
                match ias15.iterate(&mut output_bin) {
                    Ok(_) => {},
                    Err(e) => { println!("{}", e); break; }
                };
            },
            posidonius::IntegratorType::WHFastHelio => {
                match whfasthelio.iterate(&mut output_bin) {
                    Ok(_) => {},
                    Err(e) => { println!("{}", e); break; }
                };
            }
        };
    }

    let t2 = time::precise_time_s();
    let d = t2 - t1;
    println!("Execution time: {} seconds", d);
}

