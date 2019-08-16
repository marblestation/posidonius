use super::{Particle};
use super::super::constants::{SUN_DYN_FREQ};
use super::super::tools::{interpolate_b_spline};
use super::super::{csv};
use super::super::constants::{R_SUN, M2AU};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum EvolutionType {
    GalletBolmont2017(f64), // SolarLike EvolvingDissipation Evolving dissipation
    BolmontMathis2016(f64), // SolarLike EvolvingDissipation Evolving dissipation
    Baraffe2015(f64), // NEW
    Leconte2011(f64), // BrownDwarf
    Baraffe1998(f64), // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    LeconteChabrier2013(bool), // Jupiter with/without dissipation of dynamical tides
    NonEvolving,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct Evolver {
    pub evolution_type: EvolutionType,
    pub time: Vec<f64>,
    pub radius: Vec<f64>,
    pub radius_of_gyration_2: Vec<f64>,
    pub love_number: Vec<f64>,
    pub inverse_tidal_q_factor: Vec<f64>, // Bolmont & Mathis 2016
    left_index: usize,
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
            left_index: 0,
        }

    }
}

impl Evolver {
    //pub fn new_dummy() -> Evolver {
        //Evolver{
            //evolution_type: EvolutionType::NonEvolving,
            //time: Vec::new(),
            //radius: Vec::new(),
            //radius_of_gyration_2: Vec::new(),
            //love_number: Vec::new(),
            //inverse_tidal_q_factor: Vec::new(),
            //left_index: 0,
        //}
    //}
    pub fn new(evolution_type: EvolutionType, initial_time: f64, time_limit: f64) -> Evolver {
        let mut time: Vec<f64> = Vec::new();
        let mut radius: Vec<f64> = Vec::new();
        let mut radius_of_gyration_2: Vec<f64> = Vec::new();
        let mut love_number: Vec<f64> = Vec::new();
        let mut inverse_tidal_q_factor: Vec<f64> = Vec::new();

        let filename = match evolution_type {
            EvolutionType::GalletBolmont2017(mass) => {
                if mass <= 0.301 && mass >= 0.299 {
                    String::from("input/Gallet_Bolmont_2017/M_03_Z_0134.dat")
                } else if mass <= 0.401 && mass >= 0.399 {
                    String::from("input/Gallet_Bolmont_2017/M_04_Z_0134.dat")
                } else if mass <= 0.601 && mass >= 0.599 {
                    String::from("input/Gallet_Bolmont_2017/M_06_Z_0134.dat")
                } else if mass <= 0.701 && mass >= 0.699 {
                    String::from("input/Gallet_Bolmont_2017/M_07_Z_0134.dat")
                } else if mass <= 0.801 && mass >= 0.799 {
                    String::from("input/Gallet_Bolmont_2017/M_08_Z_0134.dat")
                } else if mass <= 0.901 && mass >= 0.899 {
                    String::from("input/Gallet_Bolmont_2017/M_09_Z_0134.dat")
                } else if mass <= 1.001 && mass >= 0.999 {
                    String::from("input/Gallet_Bolmont_2017/M_10_Z_0134.dat")
                } else if mass <= 1.101 && mass >= 1.099 {
                    String::from("input/Gallet_Bolmont_2017/M_11_Z_0134.dat")
                } else if mass <= 1.201 && mass >= 1.199 {
                    String::from("input/Gallet_Bolmont_2017/M_12_Z_0134.dat")
                } else if mass <= 1.401 && mass >= 1.399 {
                    String::from("input/Gallet_Bolmont_2017/M_14_Z_0134.dat")
                } else {
                    panic!("The evolution type Gallet_Bolmont_2017 does not support a mass of {} Msun!", mass);
                }
            },
            EvolutionType::Baraffe2015(mass) => {
                if mass <= 0.0101 && mass >= 0.0099 {
                    String::from("input/Baraffe_2015/0010_Msun.dat")
                } else if mass <= 0.0151 && mass >= 0.0149 {
                    String::from("input/Baraffe_2015/0015_Msun.dat")
                } else if mass <= 0.0201 && mass >= 0.0199 {
                    String::from("input/Baraffe_2015/0020_Msun.dat")
                } else if mass <= 0.0301 && mass >= 0.0299 {
                    String::from("input/Baraffe_2015/0030_Msun.dat")
                } else if mass <= 0.0401 && mass >= 0.0399 {
                    String::from("input/Baraffe_2015/0040_Msun.dat")
                } else if mass <= 0.0501 && mass >= 0.0499 {
                    String::from("input/Baraffe_2015/0050_Msun.dat")
                } else if mass <= 0.0601 && mass >= 0.0599 {
                    String::from("input/Baraffe_2015/0060_Msun.dat")
                } else if mass <= 0.0701 && mass >= 0.0699 {
                    String::from("input/Baraffe_2015/0070_Msun.dat")
                } else if mass <= 0.0721 && mass >= 0.0719 {
                    String::from("input/Baraffe_2015/0072_Msun.dat")
                } else if mass <= 0.0751 && mass >= 0.0749 {
                    String::from("input/Baraffe_2015/0075_Msun.dat")
                } else if mass <= 0.0801 && mass >= 0.0799 {
                    String::from("input/Baraffe_2015/0080_Msun.dat")
                } else if (mass - 0.09).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0090_Msun.dat")
                } else if (mass - 0.11).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0110_Msun.dat")
                } else if (mass - 0.13).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0130_Msun.dat")
                } else if (mass - 0.15).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0150_Msun.dat")
                } else if (mass - 0.17).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0170_Msun.dat")
                } else if (mass - 0.20).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0200_Msun.dat")
                } else if (mass - 0.30).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0300_Msun.dat")
                } else if (mass - 0.40).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0400_Msun.dat")
                } else if (mass - 0.50).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0500_Msun.dat")
                } else if (mass - 0.60).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0600_Msun.dat")
                } else if (mass - 0.70).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0700_Msun.dat")
                } else if (mass - 0.80).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0800_Msun.dat")
                } else if (mass - 0.90).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/0900_Msun.dat")
                } else if (mass - 1.00).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/1000_Msun.dat")
                } else if (mass - 1.10).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/1100_Msun.dat")
                } else if (mass - 1.20).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/1200_Msun.dat")
                } else if (mass - 1.30).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/1300_Msun.dat")
                } else if (mass - 1.40).abs() < 1e-7 {
                    String::from("input/Baraffe_2015/1400_Msun.dat")
                } else {
                    panic!("The evolution type Baraffe2015 does not support a mass of {} Msun!", mass);
                }
            },
            EvolutionType::Leconte2011(mass) => {
                if mass <= 0.0101 && mass >= 0.0099 {
                    String::from("input/Leconte_2011/mass_10.0000.dat")
                } else if mass <= 0.0121 && mass >= 0.0119 {
                    String::from("input/Leconte_2011/mass_12.0000.dat")
                } else if mass <= 0.0151 && mass >= 0.0149 {
                    String::from("input/Leconte_2011/mass_15.0000.dat")
                } else if mass <= 0.0201 && mass >= 0.0199 {
                    String::from("input/Leconte_2011/mass_20.0000.dat")
                } else if mass <= 0.0301 && mass >= 0.0299 {
                    String::from("input/Leconte_2011/mass_30.0000.dat")
                } else if mass <= 0.0401 && mass >= 0.0399 {
                    String::from("input/Leconte_2011/mass_40.0000.dat")
                } else if mass <= 0.0501 && mass >= 0.0499 {
                    String::from("input/Leconte_2011/mass_50.0000.dat")
                } else if mass <= 0.0601 && mass >= 0.0599 {
                    String::from("input/Leconte_2011/mass_60.0000.dat")
                } else if mass <= 0.0701 && mass >= 0.0699 {
                    String::from("input/Leconte_2011/mass_70.0000.dat")
                } else if mass <= 0.0721 && mass >= 0.0719 {
                    String::from("input/Leconte_2011/mass_72.0000.dat")
                } else if mass <= 0.0751 && mass >= 0.0749 {
                    String::from("input/Leconte_2011/mass_75.0000.dat")
                } else if mass <= 0.0801 && mass >= 0.0799 {
                    String::from("input/Leconte_2011/mass_80.0000.dat")
                } else {
                    panic!("The evolution type Leconte2011 does not support a mass of {} Msun!", mass);
                }
            },
            EvolutionType::Baraffe1998(mass) => {
                if (mass - 0.10).abs() <= 1.0e-7 {
                    String::from("input/Baraffe_1998/01Msun.dat")
                } else if (mass - 1.0).abs() <= 1.0e-7 {
                    String::from("input/Baraffe_1998/SRad_Spli_M-1_0000.dat")
                } else {
                    panic!("The evolution type Baraffe1998 does not support a mass of {} Msun!", mass);
                }
            },
            EvolutionType::BolmontMathis2016(mass) => {
                if mass <= 0.401 && mass >= 0.399 {
                    String::from("input/Bolmont_Mathis_2016/L04Z02r.dat")
                } else if mass <= 0.501 && mass >= 0.499 {
                    String::from("input/Bolmont_Mathis_2016/L05Z02r.dat")
                } else if mass <= 0.601 && mass >= 0.599 {
                    String::from("input/Bolmont_Mathis_2016/L06Z02r.dat")
                } else if mass <= 0.701 && mass >= 0.699 {
                    String::from("input/Bolmont_Mathis_2016/L07Z02r.dat")
                } else if mass <= 0.801 && mass >= 0.799 {
                    String::from("input/Bolmont_Mathis_2016/L08Z02r.dat")
                } else if mass <= 0.901 && mass >= 0.899 {
                    String::from("input/Bolmont_Mathis_2016/L09Z02r.dat")
                } else if mass <= 1.001 && mass >= 0.999 {
                    String::from("input/Bolmont_Mathis_2016/L10Z02r.dat")
                } else if mass <= 1.101 && mass >= 1.099 {
                    String::from("input/Bolmont_Mathis_2016/L11Z02r.dat")
                } else if mass <= 1.201 && mass >= 1.199 {
                    String::from("input/Bolmont_Mathis_2016/L12Z02r.dat")
                } else if mass <= 1.201 && mass >= 1.499 {
                    String::from("input/Bolmont_Mathis_2016/L13Z02r.dat")
                } else if mass <= 1.301 && mass >= 1.299 {
                    String::from("input/Bolmont_Mathis_2016/L14Z02r.dat")
                } else if mass <= 1.401 && mass >= 1.399 {
                    String::from("input/Bolmont_Mathis_2016/L15Z02r.dat")
                } else {
                    panic!("The evolution type BolmontMathis2016 does not support a mass of {} Msun!", mass);
                }
            },
            EvolutionType::LeconteChabrier2013(_) => {
                String::from("input/Leconte_Chabrier_2013/Jupiter.dat")
            },
            EvolutionType::NonEvolving => {
                String::from("input/empty.dat")
            }
        };
        //println!("Filename {}", filename);
        
        let mut rdr = csv::Reader::from_file(&filename).unwrap().has_headers(false).delimiter(b' ').flexible(true);
        for (i, row) in rdr.records().map(|r| r.unwrap()).enumerate() {
            let raw_time = row[0].parse::<f64>().unwrap_or(-1.);
            if i == 0 && raw_time < 0. {
                // First line, probably column headers
                continue
            }
            let raw_radius = match evolution_type {
                EvolutionType::Baraffe2015(_) => {
                    row[2].parse::<f64>().unwrap()
                },
                EvolutionType::GalletBolmont2017(_) => {
                    row[3].parse::<f64>().unwrap()
                },
                _ => {
                    row[1].parse::<f64>().unwrap()
                }
            };
            // All types have time and radius
            let (current_time, current_radius) = match evolution_type {
                EvolutionType::Leconte2011(_) => {
                    (raw_time, raw_radius * R_SUN)
                },
                EvolutionType::Baraffe1998(mass) => {
                    if (mass - 0.10).abs() <= 1.0e-7 {
                        (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                    } else if (mass - 1.0).abs() <= 1.0e-7 {
                        (raw_time * 365.25 - initial_time, raw_radius * M2AU)
                    } else {
                        (0., 0.)
                    }
                },
                EvolutionType::Baraffe2015(_) => {
                    (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                },
                EvolutionType::BolmontMathis2016(_) => {
                    (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                },
                EvolutionType::GalletBolmont2017(_) => {
                    (10_f64.powf(raw_time) * 365.25 - initial_time, raw_radius * R_SUN)
                },
                EvolutionType::LeconteChabrier2013(_) => {
                    (raw_time * 365.25 - initial_time, raw_radius * M2AU)
                },
                _ => {
                    (0., 0.)
                }
            };
            // Fields that only some types have...
            let current_radius_of_gyration_2 = match evolution_type {
                EvolutionType::LeconteChabrier2013(_) => {
                    let raw_radius_of_gyration_2 = row[3].parse::<f64>().unwrap();
                    raw_radius_of_gyration_2
                },
                EvolutionType::Baraffe2015(_) => {
                    let raw_radius_of_gyration_2 = row[3].parse::<f64>().unwrap();
                    raw_radius_of_gyration_2
                },
                _ => {
                    0.
                }
            };
            let current_love_number = match evolution_type {
                EvolutionType::LeconteChabrier2013(_) => {
                    let raw_love_number = row[2].parse::<f64>().unwrap();
                    raw_love_number
                },
                _ => {
                    0.
                }
            };
            let current_inverse_tidal_q_factor = match evolution_type {
                EvolutionType::BolmontMathis2016(_) => {
                    let raw_inverse_tidal_q_factor = row[2].parse::<f64>().unwrap();
                    raw_inverse_tidal_q_factor
                },
                EvolutionType::GalletBolmont2017(_) => {
                    let raw_inverse_tidal_q_factor = row[10].parse::<f64>().unwrap();
                    1./10_f64.powf(raw_inverse_tidal_q_factor)
                },
                EvolutionType::LeconteChabrier2013(true) => {
                    let raw_inverse_tidal_q_factor = row[5].parse::<f64>().unwrap();
                    raw_inverse_tidal_q_factor
                },
                _ => {
                    0.
                }
            };

            // Only save values needed for the evolution model to save memory
            // and computation time when cloning
            match evolution_type {
                EvolutionType::GalletBolmont2017(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    inverse_tidal_q_factor.push(current_inverse_tidal_q_factor);
                },
                EvolutionType::Baraffe2015(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    radius_of_gyration_2.push(current_radius_of_gyration_2);
                },
                EvolutionType::Leconte2011(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    // Read radius of gyration from an auxiliary file later on
                    //radius_of_gyration_2.push(current_radius_of_gyration_2);
                },
                EvolutionType::Baraffe1998(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                },
                EvolutionType::BolmontMathis2016(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    inverse_tidal_q_factor.push(current_inverse_tidal_q_factor);
                },
                EvolutionType::LeconteChabrier2013(dissipation_of_dynamical_tides) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    radius_of_gyration_2.push(current_radius_of_gyration_2);
                    love_number.push(current_love_number);
                    if dissipation_of_dynamical_tides {
                        inverse_tidal_q_factor.push(current_inverse_tidal_q_factor);
                    }
                },
                EvolutionType::NonEvolving => {
                }
            };
        }
        
        // Leconte2011 have a separate file for radius of gyration with a different time sampling:
        if let EvolutionType::Leconte2011(mass) = evolution_type {
            // Read separate file
            let mut aux_time: Vec<f64> = Vec::new();
            let mut aux_radius_of_gyration_2: Vec<f64> = Vec::new();
            // The column to read from the auxiliary file depends on the mass
            let aux_column = {
                if mass <= 0.0101 && mass >= 0.0099 {
                    1
                } else if mass <= 0.0121 && mass >= 0.0119 {
                    2
                } else if mass <= 0.0151 && mass >= 0.0149 {
                    3
                } else if mass <= 0.0201 && mass >= 0.0199 {
                    4
                } else if mass <= 0.0301 && mass >= 0.0299 {
                    5
                } else if mass <= 0.0401 && mass >= 0.0399 {
                    6
                } else if mass <= 0.0501 && mass >= 0.0499 {
                    7
                } else if mass <= 0.0601 && mass >= 0.0599 {
                    8
                } else if mass <= 0.0701 && mass >= 0.0699 {
                    9
                } else if mass <= 0.0721 && mass >= 0.0719 {
                    10
                } else if mass <= 0.0751 && mass >= 0.0749 {
                    11
                } else if mass <= 0.0801 && mass >= 0.0799 {
                    12
                } else {
                    panic!("The evolution type Leconte2011 does not support a mass of {} Msun!", mass);
                }
            };
            let mut rdr = csv::Reader::from_file("input/Leconte_2011/rg2BD.dat").unwrap().has_headers(false).delimiter(b' ').flexible(true);
            for row in rdr.records().map(|r| r.unwrap()) {
                let raw_time = row[0].parse::<f64>().unwrap();
                let raw_radius_of_gyration_2 = row[aux_column].parse::<f64>().unwrap();
                aux_time.push(raw_time);
                aux_radius_of_gyration_2.push(raw_radius_of_gyration_2);
            }
            // The time range from the main file and the auxiliary do not match
            // the main one is bigger, so we restrain to the smallest one (the aux one)
            if let Some(lower_limit) = time.iter().position(|&x| x >= *aux_time.first().unwrap()) {
                time.drain(0..lower_limit);
                radius.drain(0..lower_limit);
            }
            if let Some(upper_limit) = time.iter().position(|&x| x > *aux_time.last().unwrap()) {
                time.drain(upper_limit..);
                radius.drain(upper_limit..);
            }
            
            // Interpolate the radius of gyration to match the same time sampling shared by other parameters (e.g., radius)
            let precision = 7;
            for current_time in time.iter() {
                let (current_radius_of_gyration_2, _) = interpolate_b_spline(&aux_time, &aux_radius_of_gyration_2, *current_time);
                radius_of_gyration_2.push(math::round::half_up(current_radius_of_gyration_2, precision));
                //println!("Rg2 {} {}", current_time, current_radius_of_gyration_2);
            }

            // Convert time to be homogenous with the rest of models
            time = time.iter().map(|&t| t * 365.25 - initial_time).collect()
        };

        if time.len() > 0 && time[0] > 0. {
            panic!("Your initial time ({} days) is smaller than the minimum allowed age of the star ({} days)", initial_time, time[0]+initial_time);
        }
        if time.len() > 0 && time[time.len()-1] < time_limit {
            panic!("Your time limit ({} days) is greater than the maximum allowed age of the star ({} days)", time_limit, time[time.len()-1]);
        }

        Evolver { evolution_type:evolution_type, 
                time:time,
                radius:radius,
                radius_of_gyration_2:radius_of_gyration_2,
                love_number:love_number,
                inverse_tidal_q_factor:inverse_tidal_q_factor,
                left_index:0,
        }
    }

    // OPTIMIZATION: Skip first N elements which belong to the past
    fn idx(&self) -> usize {
        if self.left_index > 0 {
            return self.left_index - 1;
        } else {
            return 0;
        }
    }

    pub fn radius(&mut self, current_time: f64, current_radius: f64) -> f64 {
        let (new_radius, left_index) = match self.evolution_type {
            EvolutionType::NonEvolving => { (current_radius, 0) },
            _ => { interpolate_b_spline(&self.time[self.idx()..], &self.radius[self.idx()..], current_time) }
        };
        self.left_index += left_index;
        return new_radius;
    }

    pub fn radius_of_gyration_2(&mut self, current_time: f64, current_radius_of_gyration_2: f64) -> f64 {
        let (new_radius_of_gyration_2, left_index) = match self.evolution_type {
            EvolutionType::Baraffe2015(_) => { interpolate_b_spline(&self.time[self.idx()..], &self.radius_of_gyration_2[self.idx()..], current_time) },
            EvolutionType::Leconte2011(_) => { interpolate_b_spline(&self.time[self.idx()..], &self.radius_of_gyration_2[self.idx()..], current_time) },
            EvolutionType::LeconteChabrier2013(_) => { interpolate_b_spline(&self.time[self.idx()..], &self.radius_of_gyration_2[self.idx()..], current_time) },
            _ => { (current_radius_of_gyration_2, 0) }
        };
        self.left_index += left_index;
        return new_radius_of_gyration_2;
    }

    pub fn love_number(&mut self, current_time: f64, current_love_number: f64) -> f64 {
        let (new_love_number, left_index) = match self.evolution_type {
            EvolutionType::LeconteChabrier2013(_) => { interpolate_b_spline(&self.time[self.idx()..], &self.love_number[self.idx()..], current_time) },
            _ => { (current_love_number, 0) }
        };
        self.left_index += left_index;
        return new_love_number;
    }

    pub fn inverse_tidal_q_factor(&mut self, current_time: f64, current_inverse_tidal_q_factor: f64) -> f64 {
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
        let (new_inverse_tidal_q_factor, left_index) = match self.evolution_type {
            EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
                interpolate_b_spline(&self.time[self.idx()..], &self.inverse_tidal_q_factor[self.idx()..], current_time)
            },
            _ => (current_inverse_tidal_q_factor, 0),
        };
        self.left_index += left_index;
        return new_inverse_tidal_q_factor;
    }
}


pub fn calculate_particles_evolving_quantities(current_time: f64, particles: &mut [Particle], particles_evolvers: &mut Vec<Evolver>) {
    for (particle, evolver) in particles.iter_mut().zip(particles_evolvers.iter_mut()) {
        ////////////////////////////////////////////////////////////////////
        // Radius and radius of gyration 2
        ////////////////////////////////////////////////////////////////////
        // If particle radius/radius_of_gyration_2 evolves
        // - Compute ratio between previous and new moment of inertia 
        // - The ratio will be used in the spin integration
        let new_radius = evolver.radius(current_time, particle.radius);
        let new_radius_of_gyration_2 = evolver.radius_of_gyration_2(current_time, particle.radius_of_gyration_2);
        if new_radius != particle.radius || new_radius_of_gyration_2 != particle.radius_of_gyration_2 {
            // Update moment of inertia ratio only if it is not during first initialization
            if current_time > 0. {
                particle.moment_of_inertia_ratio = (particle.radius_of_gyration_2 * particle.radius.powi(2)) / (new_radius_of_gyration_2 * new_radius.powi(2));
            }
            particle.radius = new_radius;
            particle.radius_of_gyration_2 = new_radius_of_gyration_2;
            particle.moment_of_inertia = particle.mass * particle.radius_of_gyration_2 * particle.radius.powi(2);
        } else {
            particle.moment_of_inertia_ratio = 1.;
        }

        ////////////////////////////////////////////////////////////////////
        // Love number
        ////////////////////////////////////////////////////////////////////
        particle.tides.parameters.input.love_number = evolver.love_number(current_time, particle.tides.parameters.input.love_number);

        ////////////////////////////////////////////////////////////////////
        // Lag angle
        ////////////////////////////////////////////////////////////////////
        particle.tides.parameters.internal.lag_angle = match evolver.evolution_type {
            EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
                    let inverse_tidal_q_factor = evolver.inverse_tidal_q_factor(current_time, 0.);
                    let epsilon_squared = particle.norm_spin_vector_2/SUN_DYN_FREQ;
                    // Normal formula = 3.d0*epsilon_squared*Q_str_inv/(4.d0*k2s)
                    // but as for sigma it is necessary to divide by k2s, we do not divide here
                    let lag_angle = 3.0*epsilon_squared*inverse_tidal_q_factor/4.0;
                    lag_angle
                },
            _ => 0.,
        };
        //println!("[{}] Evolve Radius {:e} Gyration {:e} Love {:e} Lag {:e}", current_time, particle.radius, particle.radius_of_gyration_2, particle.love_number, particle.lag_angle);
    }

}
