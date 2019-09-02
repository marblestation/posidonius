pub mod stars;
pub mod planets;
pub mod universe;

#[allow(dead_code)]
pub fn simulation_properties() -> (f64, f64, f64, f64, f64) {
    let time_step: f64 = 0.08; // in days
    let time_limit: f64 = time_step*200.; // days
    let initial_time: f64 = 1.2e8*365.25; // time [days] where simulation starts
    let historic_snapshot_period: f64 = 100.*365.25; // days
    let recovery_snapshot_period: f64 = 10.*historic_snapshot_period; // days
    (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period)
}

#[allow(dead_code)]
pub fn get_data_dirname(test_name: &String) -> (String, String) {
    let rust_data_dirname = format!("tests/data/{0}/", test_name);
    let python_data_dirname = format!("posidonius/tests/data/{0}/", test_name);
    ((rust_data_dirname, python_data_dirname))
}
