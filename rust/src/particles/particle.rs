use super::super::constants::{K2, TWO_PI, R_SUN, R_EARTH, M2EARTH};
use super::{EvolutionType, SolarEvolutionType};
use super::{Axes};
use super::super::tools::{calculate_spin};

#[derive(Debug, Copy, Clone, RustcEncodable, RustcDecodable, PartialEq)]
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
    //
    pub position: Axes,
    pub velocity: Axes,
    pub acceleration: Axes,
    //
    pub radial_velocity: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub distance: f64,
    // Tides
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub denergy_dt: f64,
    pub torque: Axes, // Force
    pub spin: Axes,
    pub dspin_dt: Axes,
    pub tidal_acceleration: Axes,
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
    pub fn new(mass: f64, radius: f64, dissipation_factor: f64, dissipation_factor_scale: f64, radius_of_gyration_2: f64, love_number: f64, fluid_love_number: f64, position: Axes, velocity: Axes, acceleration: Axes, spin: Axes, evolution_type: EvolutionType) -> Particle {
        let torque = Axes{x: 0., y: 0., z: 0.};
        let dspin_dt = Axes{x: 0., y: 0., z: 0.};
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        let acceleration_induced_by_rotational_flattering = Axes{x: 0., y: 0., z: 0.};
        let general_relativity_acceleration = Axes{x: 0., y: 0., z: 0.};
        let id = 0; // To be set by the universe
        match evolution_type {
            EvolutionType::BrownDwarf(_) => println!("WARNING: Bodies with BrownDwarf evolution will ignore initial radius and radius of gyration."),
            EvolutionType::MDwarf => println!("WARNING: Bodies with MDwarf evolution will ignore initial radius."),
            EvolutionType::SolarLike(model) => {
                match model {
                    SolarEvolutionType::ConstantDissipation => println!("WARNING: Bodies with SolarLike evolution will ignore initial radius."),
                    SolarEvolutionType::EvolvingDissipation(_) => println!("WARNING: Bodies with MathisSolarLike evolution will ignore initial radius and dissipation factor."),
                }
            },
            EvolutionType::Jupiter => println!("WARNING: Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number."),
            EvolutionType::NonEvolving => {},
        }
        Particle { id:id, mass:mass, mass_g: mass*K2, radius: radius, 
                    scaled_dissipation_factor:dissipation_factor_scale*dissipation_factor, 
                    dissipation_factor_scale:dissipation_factor_scale,
                    radius_of_gyration_2:radius_of_gyration_2, 
                    love_number:love_number, fluid_love_number:fluid_love_number,
                    position:position, velocity:velocity, acceleration:acceleration, spin:spin,
                    radial_velocity: 0., norm_velocity_vector:0., norm_velocity_vector_2:0., distance:0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide:0., orthogonal_component_of_the_tidal_force_due_to_planetary_tide:0., radial_component_of_the_tidal_force:0., 
                    denergy_dt:0., torque:torque, dspin_dt:dspin_dt, 
                    tidal_acceleration:tidal_acceleration, 
                    acceleration_induced_by_rotational_flattering:acceleration_induced_by_rotational_flattering,
                    general_relativity_acceleration:general_relativity_acceleration,
                    general_relativity_factor: 0.,
                    evolution_type:evolution_type,
                    lag_angle:0., // It will be initialized the first time evolve is called
        }
    }

    pub fn new_dummy() -> Particle {
        let mass = 0.;
        let radius = 0.;
        let dissipation_factor = 0.;
        let dissipation_factor_scale = 0.;
        let radius_of_gyration_2 = 0.;
        let love_number = 0.;
        let fluid_love_number = 0.;
        let position = Axes{x: 0., y: 0., z: 0.};
        let velocity = Axes{x: 0., y: 0., z: 0.};
        let acceleration = Axes{x: 0., y: 0., z: 0.};
        let spin = Axes{x: 0., y: 0., z: 0.};
        let evolution_type = EvolutionType::NonEvolving;
        let torque = Axes{x: 0., y: 0., z: 0.};
        let dspin_dt = Axes{x: 0., y: 0., z: 0.};
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        let acceleration_induced_by_rotational_flattering = Axes{x: 0., y: 0., z: 0.};
        let general_relativity_acceleration = Axes{x: 0., y: 0., z: 0.};
        let id = 0;
        Particle { id:id, mass:mass, mass_g: mass*K2, radius: radius, 
                    scaled_dissipation_factor:dissipation_factor_scale*dissipation_factor, 
                    dissipation_factor_scale:dissipation_factor_scale,
                    radius_of_gyration_2:radius_of_gyration_2, 
                    love_number:love_number, fluid_love_number:fluid_love_number,
                    position:position, velocity:velocity, acceleration:acceleration, spin:spin,
                    radial_velocity: 0., norm_velocity_vector:0., norm_velocity_vector_2:0., distance:0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide:0., orthogonal_component_of_the_tidal_force_due_to_planetary_tide:0., radial_component_of_the_tidal_force:0., 
                    denergy_dt:0., torque:torque, dspin_dt:dspin_dt, 
                    tidal_acceleration:tidal_acceleration, 
                    acceleration_induced_by_rotational_flattering:acceleration_induced_by_rotational_flattering,
                    general_relativity_acceleration:general_relativity_acceleration,
                    general_relativity_factor: 0.,
                    evolution_type:evolution_type,
                    lag_angle:0., // It will be initialized the first time evolve is called
        }
    }

    pub fn new_brown_dwarf(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        let (rotation_period, love_number) = match evolution_type {
            EvolutionType::NonEvolving => { 
                let rotation_period: f64 = 70.0; // hours
                let love_number: f64 = 0.307; // BrownDwarf
                (rotation_period, love_number)
            },
            EvolutionType::BrownDwarf(mass) => { 
                let rotation_period; // hours
                let love_number;
                if mass <= 0.0101 && mass >= 0.0099 {
                    rotation_period = 8.0;
                    love_number = 0.3790;
                } else if mass <= 0.0121 && mass >= 0.0119 {
                    rotation_period = 13.0;
                    love_number = 0.3780;
                } else if mass <= 0.0151 && mass >= 0.0149 {
                    rotation_period = 19.0;
                    love_number = 0.3760;
                } else if mass <= 0.0201 && mass >= 0.0199 {
                    rotation_period = 24.0;
                    love_number = 0.3690;
                } else if mass <= 0.0301 && mass >= 0.0299 {
                    rotation_period = 30.0;
                    love_number = 0.3550;
                } else if mass <= 0.0401 && mass >= 0.0399 {
                    rotation_period = 36.0;
                    love_number = 0.3420;
                } else if mass <= 0.0501 && mass >= 0.0499 {
                    rotation_period = 41.0;
                    love_number = 0.3330;
                } else if mass <= 0.0601 && mass >= 0.0599 {
                    rotation_period = 47.0;
                    love_number = 0.3250;
                } else if mass <= 0.0701 && mass >= 0.0699 {
                    rotation_period = 53.0;
                    love_number = 0.3110;
                } else if mass <= 0.0721 && mass >= 0.0719 {
                    rotation_period = 58.0;
                    love_number = 0.3080;
                } else if mass <= 0.0751 && mass >= 0.0749 {
                    rotation_period = 64.0;
                    love_number = 0.3070;
                } else if mass <= 0.0801 && mass >= 0.0799 {
                    rotation_period = 70.0;
                    love_number = 0.3070;
                } else {
                    panic!("The evolution type BrownDwarf does not support a mass of {} Msun!", mass);
                }
                (rotation_period, love_number)
            },
            _ => { panic!("Evolution type should be solar like or non evolving to create a solar like body!"); }
        };

        let angular_frequency = TWO_PI/(rotation_period/24.); // days^-1
        let inclination = 0.;
        let obliquity = 0.;
        let spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity);

        let fluid_love_number = love_number;
        // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        let dissipation_factor: f64 = 2.006*3.845764e4; // -60+64

        let radius_factor: f64 = 0.845649342247916;
        let radius: f64 = radius_factor * R_SUN;
        let radius_of_gyration_2: f64 = 1.94e-1; // Brown dwarf
        let brown_dwarf = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, 
                                                evolution_type);
        brown_dwarf
    }

    pub fn new_solar_like(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        match evolution_type {
            EvolutionType::SolarLike(_) => { },
            EvolutionType::NonEvolving => { },
            _ => { panic!("Evolution type should be solar like or non evolving to create a solar like body!"); }
        }

        let rotation_period = 8.; // hours
        let love_number = 0.03; // SolarLike
        let angular_frequency = TWO_PI/(rotation_period/24.); // days^-1
        let inclination = 0.;
        let obliquity = 0.;
        let spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity);

        let fluid_love_number = love_number;
        // Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        let dissipation_factor: f64 = 4.992*3.845764e-2; // -66+64

        let radius_factor: f64 = 1.;
        let radius: f64 = radius_factor * R_SUN;
        let radius_of_gyration_2: f64 = 5.9e-2; // Sun
        let solarlike_star = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, 
                                                evolution_type);
        solarlike_star
    }

    pub fn new_m_dwarf(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        match evolution_type {
            EvolutionType::MDwarf => { },
            EvolutionType::NonEvolving => { },
            _ => { panic!("Evolution type should be M dwarf or non evolving to create a M dwarf body!"); }
        }

        let rotation_period = 70.; // hours
        let love_number: f64 = 0.307; // M Dwarf
        let fluid_love_number = love_number;

        let angular_frequency = TWO_PI/(rotation_period/24.); // days^-1
        let inclination = 0.;
        let obliquity = 0.;
        let spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity);

        // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        let dissipation_factor: f64 = 2.006*3.845764e4; // -60+64

        let radius_factor: f64 = 0.845649342247916;
        let radius: f64 = radius_factor * R_SUN;
        let radius_of_gyration_2: f64 = 2.0e-1; // M-dwarf
        let m_dwarf = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, 
                                                evolution_type);
        m_dwarf
    }

    pub fn new_jupiter_like(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        match evolution_type {
            EvolutionType::Jupiter => { },
            EvolutionType::NonEvolving => { },
            _ => { panic!("Evolution type should be jupyter or non evolving to create a jupiter body!"); }
        }

        let rotation_period = 9.8; // hours
        let love_number: f64 = 0.380; // Gas giant
        let fluid_love_number = love_number;

        let angular_frequency = TWO_PI/(rotation_period/24.); // days^-1
        let inclination = 0.;
        let obliquity = 0.;
        let spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity);

        let radius_factor: f64 = 10.9; // Jupiter in R_EARTH
        let radius: f64 = radius_factor * R_EARTH;

        // TODO: What k2pdelta/dissipation_factor is the recommended?
        //let k2pdelta: f64 = 8.101852e-9; // Gas giant
        let k2pdelta: f64 = 2.893519e-7; // Gas giant for Jupiter: 2-3d-2 s, here in day (Leconte)
        let dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * radius.powi(5));
        //let dissipation_factor: f64 = 2.006*3.845764e4; // Gas giant

        let radius_of_gyration_2: f64 = 2.54e-1; // Gas giant

        let jupiter = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, 
                                                evolution_type);
        jupiter
    }

    pub fn new_earth_like(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes) -> Particle {
        let evolution_type = EvolutionType::NonEvolving;

        let rotation_period = 24.; // hours
        let love_number: f64 = 0.299; // Earth
        let fluid_love_number = 0.9532; // Earth

        let angular_frequency = TWO_PI/(rotation_period/24.); // days^-1
        let inclination = 0.;
        let obliquity = 0.;
        let spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity);

        // Earth-like => mass-radius relationship from Fortney 2007
        let radius_factor : f64 = (0.0592*0.7+0.0975) * (mass.log10() + M2EARTH.log10() - K2.log10()).powi(2) 
                                 + (0.2337*0.7+0.4938) * (mass.log10() + M2EARTH.log10() - K2.log10()) 
                                 + 0.3102*0.7+0.7932;
        let radius: f64 = radius_factor * R_EARTH;
        let radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet
        let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets
        let dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * radius.powi(5));

        let earth_like = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, 
                                                evolution_type);
        earth_like
    }

    pub fn new_terrestrial(mass: f64, radius_factor: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes) -> Particle {
        let mut terrestrial = Particle::new_earth_like(mass, dissipation_factor_scale, position, velocity, acceleration);
        // Change radius
        terrestrial.radius = radius_factor * R_EARTH;
        // The dissipation factor depends on the radius and it has changed
        let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets
        let dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * terrestrial.radius.powi(5));
        terrestrial.scaled_dissipation_factor = terrestrial.dissipation_factor_scale*dissipation_factor;
        terrestrial
    }
    
    pub fn new_gas_giant(mass: f64, radius_factor: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes) -> Particle {
        let evolution_type = EvolutionType::NonEvolving;
        let mut gas_giant = Particle::new_jupiter_like(mass, dissipation_factor_scale, position, velocity, acceleration, evolution_type);
        // Change radius
        gas_giant.radius = radius_factor * R_EARTH;
        gas_giant
    }




}


