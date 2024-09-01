extern crate time;
use serde::{Serialize, Deserialize};
use super::super::constants::{K2};
use super::{Axes};
use super::super::{Tides, RotationalFlattening, GeneralRelativity, Disk, Wind, EvolutionType};
use super::super::{TidesEffect, TidalModel, RotationalFlatteningEffect, RotationalFlatteningModel, GeneralRelativityEffect, DiskEffect, WindEffect};
use time::{OffsetDateTime, format_description};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum Reference {
    MostMassiveParticle,
    Particle(usize), // Index
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Particle {
    pub id: usize, // Unique internal identifier
    pub mass: f64,
    pub mass_g: f64,
    pub radius: f64,
    // Inertial frame where the center of mass of the system is at rest with respect to the origin of the coordinate system
    // (i.e., barycentric frame)
    pub inertial_position: Axes,
    pub inertial_velocity: Axes,
    pub inertial_acceleration: Axes,
    pub inertial_additional_acceleration: Axes,
    // Positions/velocities in a heliocentric frame: the host is at rest with respect to the origin of the coordinate system
    // where the host is the most massive particle in the universe
    pub heliocentric_position: Axes,
    pub heliocentric_velocity: Axes,
    pub heliocentric_distance: f64, // Optimization
    pub heliocentric_radial_velocity: f64, // Optimization
    pub heliocentric_norm_velocity_vector: f64, // Optimization
    pub heliocentric_norm_velocity_vector_2: f64, // Optimization
    // Spin
    pub spin: Axes,
    pub norm_spin_vector_2: f64,
    pub angular_momentum: Axes,
    pub dangular_momentum_dt: Axes, // Force
    pub radius_of_gyration_2: f64,  // radius of gyration square can be computed in terms of the mass moment of inertia, which 
                                    // depends on the shape of the body and determines the torque needed for a desired angular acceleration
    pub moment_of_inertia: f64, // Spin related
    //
    pub reference: Reference, // Particle of reference for computing keplerian orbital parameters
    //
    pub tides: Tides,
    pub rotational_flattening: RotationalFlattening,
    pub general_relativity: GeneralRelativity,
    pub wind: Wind,
    pub disk: Disk,
    pub evolution: EvolutionType,
}

impl Particle {

    pub fn new(mass: f64, radius: f64, radius_of_gyration: f64, position: Axes, velocity: Axes, spin: Axes) -> Particle {
        // Default effects: None
        let k_factor = 0.;
        let rotation_saturation = 0.;
        let tides = Tides::new(TidesEffect::Disabled);
        let rotational_flattening = RotationalFlattening::new(RotationalFlatteningEffect::Disabled);
        let general_relativity = GeneralRelativity::new(GeneralRelativityEffect::Disabled);
        let wind = Wind::new(WindEffect::Disabled, k_factor, rotation_saturation);
        let disk = Disk::new(DiskEffect::Disabled);
        let evolution = EvolutionType::NonEvolving;
        let radius_of_gyration_2 = radius_of_gyration.powi(2);
        let moment_of_inertia = mass * radius_of_gyration_2 * radius.powi(2);
        Particle { 
            id: 0, // Unique internal identifier, to be set by the universe
            mass: mass,
            mass_g: mass*K2,
            radius: radius,
            inertial_position: Axes{x: 0., y: 0., z: 0.}, // To be re-computed by the universe
            inertial_velocity: Axes{x: 0., y: 0., z: 0.}, // To be re-computed by the universe
            inertial_acceleration: Axes{x: 0., y: 0., z: 0.},
            inertial_additional_acceleration: Axes{x: 0., y: 0., z: 0.},
            heliocentric_position: position,
            heliocentric_velocity: velocity,
            heliocentric_distance: 0.,
            heliocentric_radial_velocity: 0.,
            heliocentric_norm_velocity_vector: 0.,
            heliocentric_norm_velocity_vector_2: 0.,
            spin: spin,
            norm_spin_vector_2: (spin.x.powi(2)) + (spin.y.powi(2)) + (spin.z.powi(2)),
            angular_momentum: Axes{x: moment_of_inertia*spin.x, y: moment_of_inertia*spin.y, z: moment_of_inertia*spin.z},
            dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
            radius_of_gyration_2: radius_of_gyration_2,
            moment_of_inertia: moment_of_inertia,
            reference: Reference::MostMassiveParticle,
            tides: tides,
            rotational_flattening: rotational_flattening,
            general_relativity: general_relativity,
            wind: wind,
            disk: disk,
            evolution: evolution,
        }
    }

    pub fn new_dummy() -> Particle {
        let k_factor = 0.;
        let rotation_saturation = 0.;
        Particle { 
            id: 0,
            mass: 0.,
            mass_g: 0.,
            radius: 0.,
            inertial_position: Axes{x: 0., y: 0., z: 0.},
            inertial_velocity: Axes{x: 0., y: 0., z: 0.},
            inertial_acceleration: Axes{x: 0., y: 0., z: 0.},
            inertial_additional_acceleration: Axes{x: 0., y: 0., z: 0.},
            heliocentric_position: Axes{x: 0., y: 0., z: 0.},
            heliocentric_velocity: Axes{x: 0., y: 0., z: 0.},
            heliocentric_distance: 0.,
            heliocentric_radial_velocity: 0.,
            heliocentric_norm_velocity_vector: 0.,
            heliocentric_norm_velocity_vector_2: 0.,
            spin: Axes{x: 0., y: 0., z: 0.},
            norm_spin_vector_2: 0.,
            angular_momentum: Axes{x: 0., y: 0., z: 0.},
            dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
            radius_of_gyration_2: 0.,
            moment_of_inertia: 0.,
            reference: Reference::MostMassiveParticle,
            tides: Tides::new(TidesEffect::Disabled),
            rotational_flattening: RotationalFlattening::new(RotationalFlatteningEffect::Disabled),
            general_relativity: GeneralRelativity::new(GeneralRelativityEffect::Disabled),
            wind: Wind::new(WindEffect::Disabled, k_factor, rotation_saturation),
            disk: Disk::new(DiskEffect::Disabled),
            evolution: EvolutionType::NonEvolving,
        }
    }

    pub fn set_tides(&mut self, tides: Tides) {
        self.tides = tides;
        self.check_uniform_viscosity_coefficient();
    }

    pub fn set_rotational_flattening(&mut self, rotational_flattening: RotationalFlattening) {
        self.rotational_flattening = rotational_flattening;
        self.check_uniform_viscosity_coefficient();
    }

    pub fn check_uniform_viscosity_coefficient(&mut self) {
        // If creep coplanar tides and rotational flattening are set, both need to have the same
        // uniform viscosity coefficient parameter
        let disabled_tides = match self.tides.effect {
            TidesEffect::Disabled => true,
            _ => false
        };
        let disabled_rotational_flattening = match self.rotational_flattening.effect {
            RotationalFlatteningEffect::Disabled => true,
            _ => false
        };
        if !disabled_tides && !disabled_rotational_flattening {
            let (creep_coplanar_tides, particle_uniform_viscosity_coefficient_for_tides) = match &self.tides.effect {
                TidesEffect::CentralBody(tidal_model) | TidesEffect::OrbitingBody(tidal_model) => {
                    if let TidalModel::CreepCoplanar(params) = tidal_model {
                        (true, params.uniform_viscosity_coefficient)
                    } else {
                        (false, 0.)
                    }
                },
                _ => (false, 0.)
            };
            let (creep_coplanar_rotational_flattening, particle_uniform_viscosity_coefficient_for_rotational_flattenning) = match self.rotational_flattening.effect {
                RotationalFlatteningEffect::CentralBody(rotational_flattening_model) | RotationalFlatteningEffect::OrbitingBody(rotational_flattening_model) => {
                    if let RotationalFlatteningModel::CreepCoplanar(params) = rotational_flattening_model {
                        (true, params.uniform_viscosity_coefficient)
                    } else {
                        (false, 0.)
                    }
                },
                _ => (false, 0.),
            };
            if (creep_coplanar_tides && !creep_coplanar_rotational_flattening) || (!creep_coplanar_tides && creep_coplanar_rotational_flattening) {
                panic!("[ERROR {} UTC] When using Creep Coplanar Tidal or rotational flattening effects, both effects need to be Creep Coplanar and not just one of them (e.g., it cannot be mixed with ConstantTimeLag or OblateSpheroid).", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            } else if creep_coplanar_tides && creep_coplanar_rotational_flattening {
                let diff_uniform_viscosity_coefficient = (particle_uniform_viscosity_coefficient_for_tides - particle_uniform_viscosity_coefficient_for_rotational_flattenning).abs();
                if diff_uniform_viscosity_coefficient > 1.0e-16 {
                    panic!("[ERROR {} UTC] When using Creep Coplanar Tidal and rotational flattening effects, the uniform viscosity coefficient must be identical {:.16}.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap(), diff_uniform_viscosity_coefficient);
                }
            }
        }
    }

    pub fn set_general_relativity(&mut self, general_relativity: GeneralRelativity) {
        self.general_relativity = general_relativity;
    }

    pub fn set_wind(&mut self, wind: Wind) {
        self.wind = wind;
    }

    pub fn set_disk(&mut self, disk: Disk) {
        self.disk = disk;
    }

    pub fn set_evolution(&mut self, evolution: EvolutionType) {
        evolution_warnings(evolution);
        self.evolution = evolution;
    }
}

fn evolution_warnings(evolution: EvolutionType) {
    match evolution {
        EvolutionType::GalletBolmont2017(_) => {
            println!("[WARNING {} UTC] Bodies with GalletBolmont2017 evolution will ignore initial radius and dissipation factor.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            println!("[WARNING {} UTC] GalletBolmont2017 prescription theoretically only works for circular orbits and non inclined orbits, use carefully.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap())
        },
        EvolutionType::BolmontMathis2016(_) => {
            println!("[WARNING {} UTC] Bodies with Baraffe2015 evolution will ignore initial radius and radius of gyration.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            println!("[WARNING {} UTC] BolmontMathis2016 prescription theoretically only works for circular orbits and non inclined orbits, use carefully. ", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap())
        },
        EvolutionType::Baraffe2015(_) => println!("[WARNING {} UTC]  ", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap()),
        EvolutionType::Leconte2011(_) => println!("[WARNING {} UTC] Bodies with Leconte2011 evolution will ignore initial radius and radius of gyration.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap()),
        EvolutionType::Baraffe1998(_) => println!("[WARNING {} UTC] Bodies with Baraffe1998 evolution will ignore initial radius. ", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap()),
        EvolutionType::LeconteChabrier2013(false) => println!("[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap()),
        EvolutionType::LeconteChabrier2013(true) => {
            println!("[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration, love number and dissipation factor.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            println!("[WARNING {} UTC] LeconteChabrier2013(true) prescription theoretically only works for circular orbits and non inclined orbits, use carefully.", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
        },
        EvolutionType::NonEvolving => {},
    }
}

