use super::super::constants::{K2};
use super::{Axes};
use super::{Tides, RotationalFlattening, GeneralRelativity, Disk, Wind, EvolutionType};
use super::{TidesEffect, RotationalFlatteningEffect, GeneralRelativityEffect, DiskEffect, WindEffect};
use time;

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
    pub dangular_momentum_dt: Axes, // Force
    pub dangular_momentum_dt_per_moment_of_inertia: Axes,
    pub radius_of_gyration_2: f64,  // radius of gyration square can be computed in terms of the mass moment of inertia, which 
                                    // depends on the shape of the body and determines the torque needed for a desired angular acceleration
    pub moment_of_inertia_ratio: f64, // Spin related
    pub moment_of_inertia: f64, // Spin related
    //
    pub tides: Tides,
    pub rotational_flattening: RotationalFlattening,
    pub general_relativity: GeneralRelativity,
    pub wind: Wind,
    pub disk: Disk,
    pub evolution_type: EvolutionType,
}

impl Particle {
    pub fn new(mass: f64, radius: f64, radius_of_gyration_2: f64, heliocentric_position: Axes, heliocentric_velocity: Axes, spin: Axes, tides: Tides, rotational_flattening: RotationalFlattening, general_relativity: GeneralRelativity, wind: Wind, disk: Disk, evolution_type: EvolutionType) -> Particle {
        evolution_warnings(evolution_type);
        Particle { 
            id: 0, // Unique internal identifier, to be set by the universe
            mass: mass,
            mass_g: mass*K2,
            radius: radius,
            inertial_position: Axes{x: 0., y: 0., z: 0.}, // To be re-computed by the universe
            inertial_velocity: Axes{x: 0., y: 0., z: 0.}, // To be re-computed by the universe
            inertial_acceleration: Axes{x: 0., y: 0., z: 0.},
            heliocentric_position: heliocentric_position,
            heliocentric_velocity: heliocentric_velocity,
            heliocentric_distance: 0.,
            heliocentric_radial_velocity: 0.,
            heliocentric_norm_velocity_vector: 0.,
            heliocentric_norm_velocity_vector_2: 0.,
            spin: spin,
            norm_spin_vector_2: (spin.x.powi(2)) + (spin.y.powi(2)) + (spin.z.powi(2)),
            dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
            dangular_momentum_dt_per_moment_of_inertia: Axes{x: 0., y: 0., z: 0.},
            radius_of_gyration_2: radius_of_gyration_2,
            moment_of_inertia_ratio: 1.,
            moment_of_inertia: mass * radius_of_gyration_2 * radius.powi(2),
            tides: tides,
            rotational_flattening: rotational_flattening,
            general_relativity: general_relativity,
            wind: wind,
            disk: disk,
            evolution_type: evolution_type,
        }
    }

    pub fn new_dummy() -> Particle {
        let dissipation_factor = 0.;
        let dissipation_factor_scale = 0.;
        let love_number = 0.;
        let fluid_love_number = 0.;
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
            heliocentric_position: Axes{x: 0., y: 0., z: 0.},
            heliocentric_velocity: Axes{x: 0., y: 0., z: 0.},
            heliocentric_distance: 0.,
            heliocentric_radial_velocity: 0.,
            heliocentric_norm_velocity_vector: 0.,
            heliocentric_norm_velocity_vector_2: 0.,
            spin: Axes{x: 0., y: 0., z: 0.},
            norm_spin_vector_2: 0.,
            dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
            dangular_momentum_dt_per_moment_of_inertia: Axes{x: 0., y: 0., z: 0.},
            radius_of_gyration_2: 0.,
            moment_of_inertia_ratio: 1.,
            moment_of_inertia: 0.,
            tides: Tides::new(TidesEffect::None, dissipation_factor, dissipation_factor_scale, love_number),
            rotational_flattening: RotationalFlattening::new(RotationalFlatteningEffect::None, fluid_love_number),
            general_relativity: GeneralRelativity::new(GeneralRelativityEffect::None),
            wind: Wind::new(WindEffect::None, k_factor, rotation_saturation),
            disk: Disk::new(DiskEffect::None),
            evolution_type: EvolutionType::NonEvolving,
        }
    }
}

fn evolution_warnings(evolution_type: EvolutionType) {
    match evolution_type {
        EvolutionType::GalletBolmont2017(_) => {
            println!("[WARNING {} UTC] Bodies with GalletBolmont2017 evolution will ignore initial radius and dissipation factor.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap());
            println!("[WARNING {} UTC] GalletBolmont2017 prescription theoretically only works for circular orbits and non inclined orbits, use carefully.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap())
        },
        EvolutionType::BolmontMathis2016(_) => {
            println!("[WARNING {} UTC] Bodies with Baraffe2015 evolution will ignore initial radius and radius of gyration.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap());
            println!("[WARNING {} UTC] BolmontMathis2016 prescription theoretically only works for circular orbits and non inclined orbits, use carefully. ", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap())
        },
        EvolutionType::Baraffe2015(_) => println!("[WARNING {} UTC]  ", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap()),
        EvolutionType::Leconte2011(_) => println!("[WARNING {} UTC] Bodies with Leconte2011 evolution will ignore initial radius and radius of gyration.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap()),
        EvolutionType::Baraffe1998(_) => println!("[WARNING {} UTC] Bodies with Baraffe1998 evolution will ignore initial radius. ", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap()),
        EvolutionType::LeconteChabrier2013(false) => println!("[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap()),
        EvolutionType::LeconteChabrier2013(true) => {
            println!("[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration, love number and dissipation factor.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap());
            println!("[WARNING {} UTC] LeconteChabrier2013(true) prescription theoretically only works for circular orbits and non inclined orbits, use carefully.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap());
        },
        EvolutionType::NonEvolving => {},
    }
}

