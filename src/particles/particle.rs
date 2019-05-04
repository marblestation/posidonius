use super::super::constants::{K2};
use super::{EvolutionType};
use super::{Axes};
use time;

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Particle {
    pub id: usize, // Unique internal identifier
    pub mass: f64,
    pub mass_g: f64,
    pub radius: f64,
    pub scaled_dissipation_factor: f64, // sigma
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub radius_of_gyration_2: f64,  // radius of gyration square can be computed in terms of the mass moment of inertia, which 
                                    // depends on the shape of the body and determines the torque needed for a desired angular acceleration
    pub love_number: f64,   // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the 
                            // susceptibility of its shape to change in response to a tidal potential.
    pub fluid_love_number: f64,   // love number for a completely fluid planet (used for rotational flattening effects)
    pub type_two_migration_time: f64, // type two migration time
    pub type_two_migration_inner_disk_edge_distance: f64, // type two migration stops at inner disk edge
    //
    // In the heliocentric frame the star is at rest with respect to the origin of the coordinate system
    pub position: Axes,
    pub velocity: Axes,
    pub acceleration: Axes,
    // In the inertial frame the center of mass of a system is at rest with respect to the origin of the coordinate system
    // (i.e., barycentric frame)
    pub inertial_position: Axes,
    pub inertial_velocity: Axes,
    pub inertial_acceleration: Axes,
    //
    pub radial_velocity: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub norm_spin_vector_2: f64,
    pub distance: f64,

    // Tides
    pub scalar_product_of_vector_position_with_stellar_spin: f64,
    pub scalar_product_of_vector_position_with_planetary_spin: f64,
    pub radial_component_of_the_force_induced_by_rotation: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub factor_for_the_force_induced_by_star_rotation: f64,
    pub factor_for_the_force_induced_by_planet_rotation: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub orthogonal_component_of_the_force_induced_by_star_rotation: f64,
    pub orthogonal_component_of_the_force_induced_by_planet_rotation: f64,
    pub denergy_dt: f64,
    pub radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: f64, // needed to compute denergy_dt
    pub dangular_momentum_dt_due_to_tides: Axes, // Force
    pub dangular_momentum_dt_induced_by_rotational_flattening: Axes, // Force
    pub dangular_momentum_dt: Axes, // Force
    pub spin: Axes,
    pub dangular_momentum_dt_per_moment_of_inertia: Axes,
    pub moment_of_inertia_ratio: f64, // Spin related
    pub moment_of_inertia: f64, // Spin related
    //
    pub wind_k_factor: f64,
    pub wind_rotation_saturation: f64,
    pub wind_rotation_saturation_2: f64,
    pub wind_factor: f64, // Spin related
    //
    pub tidal_acceleration: Axes,
    // Type II migration
    pub type_two_migration_acceleration: Axes,
    // Rotational flattening
    pub acceleration_induced_by_rotational_flattering: Axes,
    // General Relativity
    pub general_relativity_factor: f64,
    pub general_relativity_acceleration: Axes,
    // Evolution
    pub evolution_type: EvolutionType,
    pub lag_angle: f64, // MathisSolarLike
}

impl Particle {
    pub fn new(mass: f64, radius: f64, dissipation_factor: f64, dissipation_factor_scale: f64, radius_of_gyration_2: f64, love_number: f64, fluid_love_number: f64, position: Axes, velocity: Axes, acceleration: Axes, spin: Axes, evolution_type: EvolutionType, wind_k_factor: f64, wind_rotation_saturation: f64) -> Particle {
        let inertial_position = Axes{x: 0., y: 0., z: 0.};
        let inertial_velocity = Axes{x: 0., y: 0., z: 0.};
        let inertial_acceleration = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt_due_to_tides = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt_induced_by_rotational_flattening = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt_per_moment_of_inertia = Axes{x: 0., y: 0., z: 0.};
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        let moment_of_inertia = mass * radius_of_gyration_2 * radius*radius;
        let acceleration_induced_by_rotational_flattering = Axes{x: 0., y: 0., z: 0.};
        let general_relativity_acceleration = Axes{x: 0., y: 0., z: 0.};
        let norm_spin_vector_2 = (spin.x.powi(2)) + (spin.y.powi(2)) + (spin.z.powi(2));
        let id = 0; // To be set by the universe
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
            EvolutionType::LeconteChabrier2013 => println!("[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number.", time::now_utc().strftime("%Y.%m.%d %H:%M%S").unwrap()),
            EvolutionType::NonEvolving => {},
        }
        Particle { id:id, mass:mass, mass_g: mass*K2, radius: radius, 
                    scaled_dissipation_factor:dissipation_factor_scale*dissipation_factor, 
                    dissipation_factor_scale:dissipation_factor_scale,
                    radius_of_gyration_2:radius_of_gyration_2, 
                    love_number:love_number, fluid_love_number:fluid_love_number,
                    position:position, velocity:velocity, acceleration:acceleration, 
                    inertial_position: inertial_position, inertial_velocity: inertial_velocity, inertial_acceleration: inertial_acceleration,
                    radial_velocity: 0., norm_velocity_vector:0., norm_velocity_vector_2:0., norm_spin_vector_2: norm_spin_vector_2, distance:0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide:0., orthogonal_component_of_the_tidal_force_due_to_planetary_tide:0., 
                    factor_for_the_force_induced_by_star_rotation: 0.,
                    factor_for_the_force_induced_by_planet_rotation: 0.,
                    radial_component_of_the_tidal_force:0., 
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
                    orthogonal_component_of_the_force_induced_by_star_rotation: 0.,
                    orthogonal_component_of_the_force_induced_by_planet_rotation: 0.,
                    denergy_dt:0., 
                    radial_component_of_the_force_induced_by_rotation: 0.,
                    dangular_momentum_dt_due_to_tides: dangular_momentum_dt_due_to_tides,
                    dangular_momentum_dt_induced_by_rotational_flattening: dangular_momentum_dt_induced_by_rotational_flattening,
                    dangular_momentum_dt: dangular_momentum_dt,
                    spin:spin,
                    dangular_momentum_dt_per_moment_of_inertia: dangular_momentum_dt_per_moment_of_inertia,
                    moment_of_inertia_ratio: 1., wind_factor: 0.,
                    moment_of_inertia: moment_of_inertia,
                    wind_k_factor: wind_k_factor,
                    wind_rotation_saturation: wind_rotation_saturation,
                    wind_rotation_saturation_2: wind_rotation_saturation*wind_rotation_saturation,
                    tidal_acceleration:tidal_acceleration, 
                    acceleration_induced_by_rotational_flattering:acceleration_induced_by_rotational_flattering,
                    general_relativity_factor: 0., // To be setup by universe
                    general_relativity_acceleration:general_relativity_acceleration,
                    evolution_type:evolution_type,
                    lag_angle:0., // It will be initialized the first time evolve is called
        }
    }

    pub fn new_dummy() -> Particle {
        let mass = 0.;
        let radius = 0.;
        let dissipation_factor = 0.;
        let dissipation_factor_scale = 1.;
        let radius_of_gyration_2 = 0.;
        let love_number = 0.;
        let fluid_love_number = 0.;
        let position = Axes{x: 0., y: 0., z: 0.};
        let velocity = Axes{x: 0., y: 0., z: 0.};
        let acceleration = Axes{x: 0., y: 0., z: 0.};
        let spin = Axes{x: 0., y: 0., z: 0.};
        let evolution_type = EvolutionType::NonEvolving;
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        let acceleration_induced_by_rotational_flattering = Axes{x: 0., y: 0., z: 0.};
        let general_relativity_acceleration = Axes{x: 0., y: 0., z: 0.};
        let inertial_position = Axes{x: 0., y: 0., z: 0.};
        let inertial_velocity = Axes{x: 0., y: 0., z: 0.};
        let inertial_acceleration = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt_due_to_tides = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt_induced_by_rotational_flattening = Axes{x: 0., y: 0., z: 0.};
        let dangular_momentum_dt_per_moment_of_inertia = Axes{x: 0., y: 0., z: 0.};
        let norm_spin_vector_2 = (spin.x.powi(2)) + (spin.y.powi(2)) + (spin.z.powi(2));
        let id = 0;
        Particle { id:id, mass:mass, mass_g: mass*K2, radius: radius, 
                    scaled_dissipation_factor:dissipation_factor_scale*dissipation_factor, 
                    dissipation_factor_scale:dissipation_factor_scale,
                    radius_of_gyration_2:radius_of_gyration_2, 
                    love_number:love_number, fluid_love_number:fluid_love_number,
                    position:position, velocity:velocity, acceleration:acceleration, 
                    inertial_position: inertial_position, inertial_velocity: inertial_velocity, inertial_acceleration: inertial_acceleration,
                    radial_velocity: 0., norm_velocity_vector:0., norm_velocity_vector_2:0., norm_spin_vector_2: norm_spin_vector_2, distance:0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide:0., orthogonal_component_of_the_tidal_force_due_to_planetary_tide:0., 
                    factor_for_the_force_induced_by_star_rotation: 0.,
                    factor_for_the_force_induced_by_planet_rotation: 0.,
                    radial_component_of_the_tidal_force:0., 
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
                    orthogonal_component_of_the_force_induced_by_star_rotation: 0.,
                    orthogonal_component_of_the_force_induced_by_planet_rotation: 0.,
                    denergy_dt:0., 
                    radial_component_of_the_force_induced_by_rotation: 0.,
                    dangular_momentum_dt_due_to_tides: dangular_momentum_dt_due_to_tides,
                    dangular_momentum_dt_induced_by_rotational_flattening: dangular_momentum_dt_induced_by_rotational_flattening,
                    dangular_momentum_dt: dangular_momentum_dt,
                    spin:spin,
                    dangular_momentum_dt_per_moment_of_inertia: dangular_momentum_dt_per_moment_of_inertia,
                    moment_of_inertia_ratio: 1., wind_factor: 0.,
                    moment_of_inertia: 0.,
                    wind_k_factor: 0.,
                    wind_rotation_saturation: 0.,
                    wind_rotation_saturation_2: 0.,
                    tidal_acceleration:tidal_acceleration, 
                    acceleration_induced_by_rotational_flattering:acceleration_induced_by_rotational_flattering,
                    general_relativity_factor: 0.,
                    general_relativity_acceleration:general_relativity_acceleration,
                    evolution_type:evolution_type,
                    lag_angle:0., // It will be initialized the first time evolve is called
        }
    }
}
