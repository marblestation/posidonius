use super::super::constants::{R_SUN, M2AU};
use super::super::tools::{interpolate_b_spline};
use super::super::{csv};

#[derive(Debug, Copy, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub enum SolarEvolutionType {
    EvolvingDissipation(f64),
    ConstantDissipation,
}

#[derive(Debug, Copy, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub enum EvolutionType {
    BrownDwarf(f64),
    MDwarf,
    SolarLike(SolarEvolutionType),
    Jupiter,
    NonEvolving,
}

#[derive(Debug, RustcEncodable, RustcDecodable, PartialEq)]
pub struct Evolver {
    pub evolution_type: EvolutionType,
    pub time: Vec<f64>,
    pub radius: Vec<f64>,
    pub radius_of_gyration_2: Vec<f64>,
    pub love_number: Vec<f64>,
    pub inverse_tidal_q_factor: Vec<f64>, // Bolmont & Mathis 2016
}

// NOTE: This is a big optimization to reduce the cost of cloning
//       Cloned particles are downgraded to a NonEvolving status
//       to avoid cloning model's data which can be big and expensive.
impl Clone for Evolver {
    fn clone(&self) -> Self {
        Evolver {
            evolution_type:EvolutionType::NonEvolving,
            time:vec![],
            radius:vec![],
            radius_of_gyration_2:vec![],
            love_number:vec![],
            inverse_tidal_q_factor:vec![],
        }

    }
}

impl Evolver {
    pub fn new(evolution_type: EvolutionType, initial_time: f64, time_limit: f64) -> Evolver {
        let mut time: Vec<f64> = Vec::new();
        let mut radius: Vec<f64> = Vec::new();
        let mut radius_of_gyration_2: Vec<f64> = Vec::new();
        let mut love_number: Vec<f64> = Vec::new();
        let mut inverse_tidal_q_factor: Vec<f64> = Vec::new();

        let filename = match evolution_type {
            EvolutionType::BrownDwarf(mass) => {
                if mass <= 0.0101 && mass >= 0.0099 {
                    String::from("input/data_host_body/mass_10.0000.dat")
                } else if mass <= 0.0121 && mass >= 0.0119 {
                    String::from("input/data_host_body/mass_12.0000.dat")
                } else if mass <= 0.0151 && mass >= 0.0149 {
                    String::from("input/data_host_body/mass_15.0000.dat")
                } else if mass <= 0.0201 && mass >= 0.0199 {
                    String::from("input/data_host_body/mass_20.0000.dat")
                } else if mass <= 0.0301 && mass >= 0.0299 {
                    String::from("input/data_host_body/mass_30.0000.dat")
                } else if mass <= 0.0401 && mass >= 0.0399 {
                    String::from("input/data_host_body/mass_40.0000.dat")
                } else if mass <= 0.0501 && mass >= 0.0499 {
                    String::from("input/data_host_body/mass_50.0000.dat")
                } else if mass <= 0.0601 && mass >= 0.0599 {
                    String::from("input/data_host_body/mass_60.0000.dat")
                } else if mass <= 0.0701 && mass >= 0.0699 {
                    String::from("input/data_host_body/mass_70.0000.dat")
                } else if mass <= 0.0721 && mass >= 0.0719 {
                    String::from("input/data_host_body/mass_72.0000.dat")
                } else if mass <= 0.0751 && mass >= 0.0749 {
                    String::from("input/data_host_body/mass_75.0000.dat")
                } else if mass <= 0.0801 && mass >= 0.0799 {
                    String::from("input/data_host_body/mass_80.0000.dat")
                } else {
                    panic!("The evolution type BrownDwarf does not support a mass of {} Msun!", mass);
                }
            },
            EvolutionType::MDwarf => {
                String::from("input/data_host_body/01Msun.dat")
            },
            EvolutionType::SolarLike(model) => {
                match model {
                    SolarEvolutionType::ConstantDissipation => String::from("input/data_host_body/SRad_Spli_M-1_0000.dat"),
                    SolarEvolutionType::EvolvingDissipation(mass) => {
                        if mass <= 0.401 && mass >= 0.399 {
                            String::from("input/data_host_body/L04Z02r.dat")
                        } else if mass <= 0.501 && mass >= 0.499 {
                            String::from("input/data_host_body/L05Z02r.dat")
                        } else if mass <= 0.601 && mass >= 0.599 {
                            String::from("input/data_host_body/L06Z02r.dat")
                        } else if mass <= 0.701 && mass >= 0.699 {
                            String::from("input/data_host_body/L07Z02r.dat")
                        } else if mass <= 0.801 && mass >= 0.799 {
                            String::from("input/data_host_body/L08Z02r.dat")
                        } else if mass <= 0.901 && mass >= 0.899 {
                            String::from("input/data_host_body/L09Z02r.dat")
                        } else if mass <= 1.001 && mass >= 0.999 {
                            String::from("input/data_host_body/L10Z02r.dat")
                        } else if mass <= 1.101 && mass >= 1.099 {
                            String::from("input/data_host_body/L11Z02r.dat")
                        } else if mass <= 1.201 && mass >= 1.199 {
                            String::from("input/data_host_body/L12Z02r.dat")
                        } else if mass <= 1.201 && mass >= 1.499 {
                            String::from("input/data_host_body/L13Z02r.dat")
                        } else if mass <= 1.301 && mass >= 1.299 {
                            String::from("input/data_host_body/L14Z02r.dat")
                        } else if mass <= 1.401 && mass >= 1.399 {
                            String::from("input/data_host_body/L15Z02r.dat")
                        } else {
                            panic!("The evolution type MathisSolarLike does not support a mass of {} Msun!", mass);
                        }
                    },
                }
            },
            EvolutionType::Jupiter => {
                String::from("input/data_host_body/Jupiter.dat")
            },
            EvolutionType::NonEvolving => {
                String::from("input/data_host_body/Empty.dat")
            }
        };
        //println!("Filename {}", filename);
        
        let mut rdr = csv::Reader::from_file(filename).unwrap().has_headers(false).delimiter(b' ').flexible(true);
        for row in rdr.records().map(|r| r.unwrap()) {
            let raw_time = row[0].parse::<f64>().unwrap();
            let raw_radius = row[1].parse::<f64>().unwrap();
            // All types have time and radius
            let (current_time, current_radius) = match evolution_type {
                EvolutionType::BrownDwarf(_) => {
                    (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                },
                EvolutionType::MDwarf => {
                    (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                },
                EvolutionType::SolarLike(model) => {
                    match model {
                        SolarEvolutionType::ConstantDissipation => (raw_time * 365.25 - initial_time, raw_radius * M2AU),
                        SolarEvolutionType::EvolvingDissipation(_) => (raw_time * 365.25 - initial_time, raw_radius * R_SUN),
                    }
                },
                EvolutionType::Jupiter => {
                    (raw_time * 365.25 - initial_time, raw_radius * M2AU)
                },
                _ => {
                    (0., 0.)
                }
            };
            // Fields that only some types have...
            let current_radius_of_gyration_2 = match evolution_type {
                EvolutionType::Jupiter => {
                    let raw_radius_of_gyration_2 = row[3].parse::<f64>().unwrap();
                    raw_radius_of_gyration_2
                },
                _ => {
                    0.
                }
            };
            let current_love_number = match evolution_type {
                EvolutionType::Jupiter => {
                    let raw_love_number = row[2].parse::<f64>().unwrap();
                    raw_love_number
                },
                _ => {
                    0.
                }
            };
            let current_inverse_tidal_q_factor = match evolution_type {
                EvolutionType::SolarLike(model) => {
                    match model {
                        SolarEvolutionType::EvolvingDissipation(_) => {
                            let raw_inverse_tidal_q_factor = row[2].parse::<f64>().unwrap();
                            raw_inverse_tidal_q_factor
                        },
                        _ => 0.,
                    }
                },
                _ => {
                    0.
                }
            };

            // Only save values needed for the evolution model to save memory
            // and computation time when cloning
            match evolution_type {
                EvolutionType::BrownDwarf(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    radius_of_gyration_2.push(current_radius_of_gyration_2);
                },
                EvolutionType::MDwarf => {
                    time.push(current_time);
                    radius.push(current_radius);
                },
                EvolutionType::SolarLike(model) => {
                    match model {
                        SolarEvolutionType::ConstantDissipation => {
                            time.push(current_time);
                            radius.push(current_radius);
                        },
                        SolarEvolutionType::EvolvingDissipation(_) => {
                            time.push(current_time);
                            radius.push(current_radius);
                            inverse_tidal_q_factor.push(current_inverse_tidal_q_factor);
                        },
                    }
                },
                EvolutionType::Jupiter => {
                    time.push(current_time);
                    radius.push(current_radius);
                    radius_of_gyration_2.push(current_radius_of_gyration_2);
                    love_number.push(current_love_number);
                },
                EvolutionType::NonEvolving => {
                }
            };
        }

        if time.len() > 0 && time[0] > 0. {
            panic!("Your initial time ({} days) is smaller than the minimum allowed age of the star ({} days)", initial_time, time[0]+initial_time);
        }
        if time.len() > 0 && time[time.len()-1] < time_limit {
            panic!("Your time limit ({} days) is greater than the maximum allowed age of the star ({} days)", time_limit, time[time.len()-1]);
        }

        // BrownDwarf have a separate file for radius of gyration with a different time sampling
        // that should be homogenized:
        let radius_of_gyration_2 = match evolution_type {
            EvolutionType::BrownDwarf(_) => {
                // Read separate file
                let mut tmp_time: Vec<f64> = Vec::new();
                let mut tmp_radius_of_gyration_2: Vec<f64> = Vec::new();
                let mut rdr = csv::Reader::from_file("input/data_host_body/rg2BD.dat").unwrap().has_headers(false).delimiter(b' ').flexible(true);
                for row in rdr.records().map(|r| r.unwrap()) {
                    let raw_time = row[0].parse::<f64>().unwrap();
                    let raw_radius_of_gyration_2 = row[2].parse::<f64>().unwrap();
                    let current_time = raw_time * 365.25 - initial_time;
                    let current_radius_of_gyration_2 = raw_radius_of_gyration_2;
                    tmp_time.push(current_time);
                    tmp_radius_of_gyration_2.push(current_radius_of_gyration_2);
                    //println!("Raw Rg2 {} {}", current_time, raw_radius_of_gyration_2);
                }

                // Interpolate to match the same time sampling shared by other parameters (e.g., radius)
                let mut resampled_radius_of_gyration_2: Vec<f64> = Vec::new();
                for current_time in time.iter() {
                    let current_radius_of_gyration_2 = interpolate_b_spline(&tmp_time, &tmp_radius_of_gyration_2, *current_time);
                    resampled_radius_of_gyration_2.push(current_radius_of_gyration_2);
                    //println!("Rg2 {} {}", current_time, current_radius_of_gyration_2);
                }
                resampled_radius_of_gyration_2
            },
            _ => {
                radius_of_gyration_2
            }
        };

        Evolver { evolution_type:evolution_type, 
                time:time,
                radius:radius,
                radius_of_gyration_2:radius_of_gyration_2,
                love_number:love_number,
                inverse_tidal_q_factor:inverse_tidal_q_factor,
        }
    }

    pub fn radius(&self, current_time: f64, current_radius: f64) -> f64 {
        let new_radius = match self.evolution_type {
            EvolutionType::NonEvolving => { current_radius },
            _ => { interpolate_b_spline(&self.time, &self.radius, current_time) }
        };
        return new_radius;
    }

    pub fn radius_of_gyration_2(&self, current_time: f64, current_radius_of_gyration_2: f64) -> f64 {
        let new_radius_of_gyration_2 = match self.evolution_type {
            EvolutionType::BrownDwarf(_) => { interpolate_b_spline(&self.time, &self.radius_of_gyration_2, current_time) },
            EvolutionType::Jupiter => { interpolate_b_spline(&self.time, &self.radius_of_gyration_2, current_time) },
            _ => { current_radius_of_gyration_2 }
        };
        return new_radius_of_gyration_2;
    }

    pub fn love_number(&self, current_time: f64, current_love_number: f64) -> f64 {
        let new_love_number = match self.evolution_type {
            EvolutionType::Jupiter => { interpolate_b_spline(&self.time, &self.love_number, current_time) },
            _ => { current_love_number }
        };
        return new_love_number;
    }

    pub fn inverse_tidal_q_factor(&self, current_time: f64, current_inverse_tidal_q_factor: f64) -> f64 {
        // The planet excite/induce waves in the convective region of the star
        // and the Q factor is a quality factor that describes that loss of energy.
        //
        // The inverse_tidal_q_factor was derived from analytical formulas using 
        // parameters from stellar models. This factor is used to derive the 
        // dissipation factor of the star, which depends on the excitation 
        // frequency of the planet (the frequency at what we would see the planet
        // if we were sitting at the surface of the star).
        //
        //  - If the star does not rotate, the excitation frequency is the orbital period (i.e.,
        //  frequency)
        //  - If the planet is very far away, the excitation frequency is the spin
        //  - Planets are treated as source points and their rotation is not important here
        //
        // The tidal theory (and this code) assumes circular orbits because they are easier.
        // Excentric orbits needs more than one frequency to be described and it will be
        // included in future versions of this code.
        let new_inverse_tidal_q_factor = match self.evolution_type {
                EvolutionType::SolarLike(model) => {
                    match model {
                        SolarEvolutionType::EvolvingDissipation(_) => {
                            interpolate_b_spline(&self.time, &self.inverse_tidal_q_factor, current_time)
                        },
                        _ => current_inverse_tidal_q_factor,
                    }
                },
                _ => current_inverse_tidal_q_factor,
        };
        return new_inverse_tidal_q_factor;
    }
}


