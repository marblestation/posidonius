use std::iter;
use super::super::constants::{G, SPEED_OF_LIGHT_2, MAX_PARTICLES, DBL_EPSILON_2};
use super::super::{Particle};
use super::super::{Axes};
use super::super::particles::universe::{IgnoreGravityTerms};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleInternalParameters {
    pub distance: f64,
    pub radial_velocity: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub factor: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleParameters {
    pub internal: GeneralRelativityParticleInternalParameters,
    pub output: GeneralRelativityParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleCoordinates {
    // Positions/velocities in a heliocentric frame 
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum GeneralRelativityImplementation {
    Kidder1995, // MercuryT
    Anderson1975, // REBOUNDx gr
    Newhall1983, // REBOUNDx gr full
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum GeneralRelativityEffect {
    CentralBody(GeneralRelativityImplementation),
    OrbitingBody,
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativity {
    pub effect: GeneralRelativityEffect,
    pub parameters: GeneralRelativityParticleParameters,
    pub coordinates: GeneralRelativityParticleCoordinates,
}

impl GeneralRelativity {
    pub fn new(effect: GeneralRelativityEffect) -> GeneralRelativity {
        GeneralRelativity {
            effect: effect,
            parameters: GeneralRelativityParticleParameters {
                internal: GeneralRelativityParticleInternalParameters {
                    distance: 0.,
                    radial_velocity: 0.,
                    norm_velocity_vector: 0.,
                    norm_velocity_vector_2: 0.,
                    factor: 0.,
                },
                output: GeneralRelativityParticleOutputParameters {
                    acceleration: Axes{x: 0., y: 0., z: 0.},
                    dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
                },
            },
            coordinates: GeneralRelativityParticleCoordinates {
                position: Axes{x: 0., y: 0., z: 0.},
                velocity: Axes{x: 0., y: 0., z: 0.},
            },
        }
    }
}

pub fn initialize(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) = host_particle.general_relativity.effect {
        if general_relativity_implementation != GeneralRelativityImplementation::Newhall1983 && general_relativity_implementation != GeneralRelativityImplementation::Disabled {
            host_particle.general_relativity.parameters.internal.factor = 0.;
            host_particle.general_relativity.parameters.output.acceleration.x = 0.;
            host_particle.general_relativity.parameters.output.acceleration.y = 0.;
            host_particle.general_relativity.parameters.output.acceleration.z = 0.;
            host_particle.general_relativity.parameters.output.dangular_momentum_dt.x = 0.;
            host_particle.general_relativity.parameters.output.dangular_momentum_dt.y = 0.;
            host_particle.general_relativity.parameters.output.dangular_momentum_dt.z = 0.;
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
                    particle.general_relativity.parameters.internal.factor = host_particle.mass_g*particle.mass_g / (host_particle.mass_g+particle.mass_g).powi(2);
                    particle.general_relativity.parameters.output.acceleration.x = 0.;
                    particle.general_relativity.parameters.output.acceleration.y = 0.;
                    particle.general_relativity.parameters.output.acceleration.z = 0.;
                    particle.general_relativity.parameters.output.dangular_momentum_dt.x = 0.;
                    particle.general_relativity.parameters.output.dangular_momentum_dt.y = 0.;
                    particle.general_relativity.parameters.output.dangular_momentum_dt.z = 0.;
                }
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) = host_particle.general_relativity.effect {
        if general_relativity_implementation != GeneralRelativityImplementation::Disabled {
            // Inertial to Heliocentric positions/velocities
            host_particle.general_relativity.coordinates.position.x = 0.;
            host_particle.general_relativity.coordinates.position.y = 0.;
            host_particle.general_relativity.coordinates.position.z = 0.;
            host_particle.general_relativity.coordinates.velocity.x = 0.;
            host_particle.general_relativity.coordinates.velocity.y = 0.;
            host_particle.general_relativity.coordinates.velocity.z = 0.;
            host_particle.general_relativity.parameters.internal.distance = 0.;
            host_particle.general_relativity.parameters.internal.radial_velocity = 0.;
            host_particle.general_relativity.parameters.internal.norm_velocity_vector = 0.;
            host_particle.general_relativity.parameters.internal.norm_velocity_vector_2 = 0.;
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
                    particle.general_relativity.coordinates.position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                    particle.general_relativity.coordinates.position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                    particle.general_relativity.coordinates.position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                    particle.general_relativity.coordinates.velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                    particle.general_relativity.coordinates.velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                    particle.general_relativity.coordinates.velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                    particle.general_relativity.parameters.internal.distance = (particle.general_relativity.coordinates.position.x.powi(2) 
                                                                        + particle.general_relativity.coordinates.position.y.powi(2)
                                                                        + particle.general_relativity.coordinates.position.z.powi(2)).sqrt();
                    particle.general_relativity.parameters.internal.radial_velocity = (particle.general_relativity.coordinates.position.x*particle.general_relativity.coordinates.velocity.x +
                                                                            particle.general_relativity.coordinates.position.y*particle.general_relativity.coordinates.velocity.y +
                                                                            particle.general_relativity.coordinates.position.z*particle.general_relativity.coordinates.velocity.z) 
                                                                            / particle.general_relativity.parameters.internal.distance;
                    particle.general_relativity.parameters.internal.norm_velocity_vector_2 = (particle.general_relativity.coordinates.velocity.x 
                                                                                       - host_particle.general_relativity.coordinates.velocity.x).powi(2) 
                                                                                    + (particle.general_relativity.coordinates.velocity.y 
                                                                                       - host_particle.general_relativity.coordinates.velocity.y).powi(2)
                                                                                    + (particle.general_relativity.coordinates.velocity.z 
                                                                                       - host_particle.general_relativity.coordinates.velocity.z).powi(2);
                    particle.general_relativity.parameters.internal.norm_velocity_vector = particle.general_relativity.parameters.internal.norm_velocity_vector_2.sqrt();
                }
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) = host_particle.general_relativity.effect {
        if general_relativity_implementation != GeneralRelativityImplementation::Disabled {
            host_particle.general_relativity.coordinates.position = host_particle.heliocentric_position;
            host_particle.general_relativity.coordinates.velocity = host_particle.heliocentric_velocity;
            host_particle.general_relativity.parameters.internal.distance = host_particle.heliocentric_distance;
            host_particle.general_relativity.parameters.internal.radial_velocity = host_particle.heliocentric_radial_velocity;
            host_particle.general_relativity.parameters.internal.norm_velocity_vector = host_particle.heliocentric_norm_velocity_vector;
            host_particle.general_relativity.parameters.internal.norm_velocity_vector_2 = host_particle.heliocentric_norm_velocity_vector_2;
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
                    particle.general_relativity.coordinates.position = particle.heliocentric_position;
                    particle.general_relativity.coordinates.velocity = particle.heliocentric_velocity;
                    particle.general_relativity.parameters.internal.distance = particle.heliocentric_distance;
                    particle.general_relativity.parameters.internal.radial_velocity = particle.heliocentric_radial_velocity;
                    particle.general_relativity.parameters.internal.norm_velocity_vector = particle.heliocentric_norm_velocity_vector;
                    particle.general_relativity.parameters.internal.norm_velocity_vector_2 = particle.heliocentric_norm_velocity_vector_2;
                }
            }
        }
    }
}


pub fn calculate_kidder1995_general_relativity_acceleration(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let mut host_particle = host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    calculate_kidder1995_first_order_general_relativity_acceleration(&mut host_particle, &mut particles, &mut more_particles);
    calculate_kidder1995_second_order_general_relativity_acceleration(&mut host_particle, &mut particles, &mut more_particles);
    calculate_kidder1995_spin_orbit_general_relativity_acceleration_and_dangular_momentum_dt(&mut host_particle, &mut particles, &mut more_particles);
}


fn calculate_kidder1995_first_order_general_relativity_acceleration(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let mut sum_total_general_relativity_acceleration = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            // Radial part of the GR force (Kidder 1995, Mardling & Lin 2002)
            // - Equation 11 from Bolmont et al. 2015
            let star_planet_mass_g = host_particle.mass_g + particle.mass_g;
            let distance_2 = particle.general_relativity.parameters.internal.distance.powi(2);
            let radial_velocity_2 = particle.general_relativity.parameters.internal.radial_velocity.powi(2);
            let radial_component_of_the_general_relativity_force = -star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                * ( (1.0 + 3.0 * particle.general_relativity.parameters.internal.factor) * particle.general_relativity.parameters.internal.norm_velocity_vector_2
                -2.0 * (2.0 + particle.general_relativity.parameters.internal.factor) * star_planet_mass_g/particle.general_relativity.parameters.internal.distance
                -1.5 * particle.general_relativity.parameters.internal.factor * radial_velocity_2);
            //println!("Radial component GR force {:e}", radial_component_of_the_general_relativity_force);
            // Orthoradial part of the GR force
            // - Equation 11 from Bolmont et al. 2015
            let orthogonal_component_of_the_general_relativity_force = star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                * 2.0 * (2.0 - particle.general_relativity.parameters.internal.factor) * particle.general_relativity.parameters.internal.radial_velocity * particle.general_relativity.parameters.internal.norm_velocity_vector;
            //println!("Ortho component GR force {:e}", orthogonal_component_of_the_general_relativity_force);
            // Total General Relativity force
            // - Equation 10 from Bolmont et al. 2015
            let total_general_relativity_acceleration_x = radial_component_of_the_general_relativity_force * particle.general_relativity.coordinates.position.x / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_general_relativity_force * particle.general_relativity.coordinates.velocity.x / particle.general_relativity.parameters.internal.norm_velocity_vector;
            let total_general_relativity_acceleration_y = radial_component_of_the_general_relativity_force * particle.general_relativity.coordinates.position.y / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_general_relativity_force * particle.general_relativity.coordinates.velocity.y / particle.general_relativity.parameters.internal.norm_velocity_vector;
            let total_general_relativity_acceleration_z = radial_component_of_the_general_relativity_force * particle.general_relativity.coordinates.position.z / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_general_relativity_force * particle.general_relativity.coordinates.velocity.z / particle.general_relativity.parameters.internal.norm_velocity_vector;
            
            sum_total_general_relativity_acceleration.x += particle.mass/host_particle.mass * total_general_relativity_acceleration_x;
            sum_total_general_relativity_acceleration.y += particle.mass/host_particle.mass * total_general_relativity_acceleration_y;
            sum_total_general_relativity_acceleration.z += particle.mass/host_particle.mass * total_general_relativity_acceleration_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.general_relativity.parameters.output.acceleration.x = total_general_relativity_acceleration_x;
            particle.general_relativity.parameters.output.acceleration.y = total_general_relativity_acceleration_y;
            particle.general_relativity.parameters.output.acceleration.z = total_general_relativity_acceleration_z;
            //println!("GR force {:e} {:e} {:e}", total_general_relativity_acceleration_x,
            //total_general_relativity_acceleration_y, total_general_relativity_acceleration_z);
        }
    }
        
    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.general_relativity.parameters.output.acceleration.x += sum_total_general_relativity_acceleration.x;
        //particle.general_relativity.parameters.output.acceleration.y += sum_total_general_relativity_acceleration.y;
        //particle.general_relativity.parameters.output.acceleration.z += sum_total_general_relativity_acceleration.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    host_particle.general_relativity.parameters.output.acceleration.x = -1.0 * sum_total_general_relativity_acceleration.x;
    host_particle.general_relativity.parameters.output.acceleration.y = -1.0 * sum_total_general_relativity_acceleration.y;
    host_particle.general_relativity.parameters.output.acceleration.z = -1.0 * sum_total_general_relativity_acceleration.z;
}

fn calculate_kidder1995_second_order_general_relativity_acceleration(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    // 2nd order Post-Newtonian
    let mut sum_total_second_order_general_relativity_acceleration = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let star_planet_mass_g = host_particle.mass_g + particle.mass_g;
            let distance_2 = particle.general_relativity.parameters.internal.distance.powi(2);
            let norm_velocity_vector_2 = particle.general_relativity.parameters.internal.norm_velocity_vector_2;
            let norm_velocity_vector_4 = norm_velocity_vector_2.powi(2);
            let radial_velocity_2 = particle.general_relativity.parameters.internal.radial_velocity.powi(2);
            let radial_velocity_4 = radial_velocity_2.powi(2);
            let general_relativity_factor_2 = particle.general_relativity.parameters.internal.factor.powi(2);

            // Radial part of the GR force (Kidder 1995, equation 2.2d)
            let radial_component_of_the_second_order_general_relativity_acceleration = -star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                    * (3.0/4.0*(12.0+29.0*particle.general_relativity.parameters.internal.factor)*(star_planet_mass_g.powi(2)/distance_2)
                    + particle.general_relativity.parameters.internal.factor*(3.0-4.0*particle.general_relativity.parameters.internal.factor)*norm_velocity_vector_4
                    + 15.0/8.0*particle.general_relativity.parameters.internal.factor*(1.0-3.0*particle.general_relativity.parameters.internal.factor)*radial_velocity_4
                    - 3.0/2.0*particle.general_relativity.parameters.internal.factor*(3.0-4.0*particle.general_relativity.parameters.internal.factor)*radial_velocity_2*norm_velocity_vector_2
                    - 0.5*particle.general_relativity.parameters.internal.factor*(13.0-4.0*particle.general_relativity.parameters.internal.factor)*(star_planet_mass_g/particle.general_relativity.parameters.internal.distance)*norm_velocity_vector_2
                    - (2.0 + 25.0*particle.general_relativity.parameters.internal.factor+2.0*general_relativity_factor_2)*(star_planet_mass_g/particle.general_relativity.parameters.internal.distance)*radial_velocity_2);

            let orthogonal_component_of_the_second_order_general_relativity_acceleration = -star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                    * (-0.5)*particle.general_relativity.parameters.internal.radial_velocity
                    * (particle.general_relativity.parameters.internal.factor*(15.0+4.0*particle.general_relativity.parameters.internal.factor)*norm_velocity_vector_2 
                    - (4.0+41.0*particle.general_relativity.parameters.internal.factor+8.0*general_relativity_factor_2)*(star_planet_mass_g/particle.general_relativity.parameters.internal.distance)
                    - 3.0*particle.general_relativity.parameters.internal.factor*(3.0+2.0*particle.general_relativity.parameters.internal.factor)*radial_velocity_2);

            let total_second_order_general_relativity_acceleration_x = radial_component_of_the_second_order_general_relativity_acceleration * particle.general_relativity.coordinates.position.x / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_second_order_general_relativity_acceleration * particle.general_relativity.coordinates.velocity.x;
            let total_second_order_general_relativity_acceleration_y = radial_component_of_the_second_order_general_relativity_acceleration * particle.general_relativity.coordinates.position.y / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_second_order_general_relativity_acceleration * particle.general_relativity.coordinates.velocity.y;
            let total_second_order_general_relativity_acceleration_z = radial_component_of_the_second_order_general_relativity_acceleration * particle.general_relativity.coordinates.position.z / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_second_order_general_relativity_acceleration * particle.general_relativity.coordinates.velocity.z;
            
            sum_total_second_order_general_relativity_acceleration.x += particle.mass/host_particle.mass * total_second_order_general_relativity_acceleration_x;
            sum_total_second_order_general_relativity_acceleration.y += particle.mass/host_particle.mass * total_second_order_general_relativity_acceleration_y;
            sum_total_second_order_general_relativity_acceleration.z += particle.mass/host_particle.mass * total_second_order_general_relativity_acceleration_z;

            particle.general_relativity.parameters.output.acceleration.x += total_second_order_general_relativity_acceleration_x;
            particle.general_relativity.parameters.output.acceleration.y += total_second_order_general_relativity_acceleration_y;
            particle.general_relativity.parameters.output.acceleration.z += total_second_order_general_relativity_acceleration_z;
            //println!("a {} {} {}", total_second_order_general_relativity_acceleration_x, total_second_order_general_relativity_acceleration_y, total_second_order_general_relativity_acceleration_z);
        }
    }
        
    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.general_relativity.parameters.output.acceleration.x += sum_total_second_order_general_relativity_acceleration.x;
        //particle.general_relativity.parameters.output.acceleration.y += sum_total_second_order_general_relativity_acceleration.y;
        //particle.general_relativity.parameters.output.acceleration.z += sum_total_second_order_general_relativity_acceleration.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    host_particle.general_relativity.parameters.output.acceleration.x += -1.0 * sum_total_second_order_general_relativity_acceleration.x;
    host_particle.general_relativity.parameters.output.acceleration.y += -1.0 * sum_total_second_order_general_relativity_acceleration.y;
    host_particle.general_relativity.parameters.output.acceleration.z += -1.0 * sum_total_second_order_general_relativity_acceleration.z;
}

fn calculate_kidder1995_spin_orbit_general_relativity_acceleration_and_dangular_momentum_dt(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    // - Equation 5 from https://arxiv.org/pdf/1102.5192.pdf
    // Spin effects are known for the dominant relativistic spin-orbit coupling term at 1.5PN
    // https://arxiv.org/pdf/gr-qc/0202016.pdf
    // Spin in Kidder is defined as angular_momentum
    let star_angular_momentum = Axes{
                                        x:host_particle.moment_of_inertia*host_particle.spin.x,
                                        y:host_particle.moment_of_inertia*host_particle.spin.y,
                                        z:host_particle.moment_of_inertia*host_particle.spin.z
    };
    let mut sum_total_general_relativity_spin_orbit_acceleration = Axes{x:0., y:0., z:0.};

    host_particle.general_relativity.parameters.output.dangular_momentum_dt.x = 0.;
    host_particle.general_relativity.parameters.output.dangular_momentum_dt.y = 0.;
    host_particle.general_relativity.parameters.output.dangular_momentum_dt.z = 0.;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            // - Equation 2.2c from Kidder 1995
            let star_planet_mass = host_particle.mass + particle.mass;
            let star_planet_diff_mass = host_particle.mass - particle.mass;
            let mass_factor = star_planet_diff_mass / star_planet_mass;

            // Spin in Kidder is defined as angular_momentum
            let particle_angular_momentum = Axes{
                                                x:particle.moment_of_inertia*particle.spin.x,
                                                y:particle.moment_of_inertia*particle.spin.y,
                                                z:particle.moment_of_inertia*particle.spin.z
            };

            let particle_normalized_position = Axes{
                                                x:particle.general_relativity.coordinates.position.x/particle.general_relativity.parameters.internal.distance,
                                                y:particle.general_relativity.coordinates.position.y/particle.general_relativity.parameters.internal.distance,
                                                z:particle.general_relativity.coordinates.position.z/particle.general_relativity.parameters.internal.distance
            };

            let mass_spin_factor_x = mass_factor*star_planet_mass*(particle_angular_momentum.x/particle.mass - star_angular_momentum.x/host_particle.mass);
            let mass_spin_factor_y = mass_factor*star_planet_mass*(particle_angular_momentum.y/particle.mass - star_angular_momentum.y/host_particle.mass);
            let mass_spin_factor_z = mass_factor*star_planet_mass*(particle_angular_momentum.z/particle.mass - star_angular_momentum.z/host_particle.mass);

            let element1_x: f64 = 6.*particle_normalized_position.x
                               * ((particle_normalized_position.y * particle.general_relativity.coordinates.velocity.z - particle_normalized_position.z * particle.general_relativity.coordinates.velocity.y)
                               * (2.*(star_angular_momentum.x+particle_angular_momentum.x) + mass_spin_factor_x));
            let element1_y :f64 = 6.*particle_normalized_position.y
                               * ((particle_normalized_position.z * particle.general_relativity.coordinates.velocity.x - particle_normalized_position.x * particle.general_relativity.coordinates.velocity.z)
                               * (2.*(star_angular_momentum.y+particle_angular_momentum.y) + mass_spin_factor_y));
            let element1_z: f64 = 6.*particle_normalized_position.z
                               * ((particle_normalized_position.x * particle.general_relativity.coordinates.velocity.y - particle_normalized_position.y * particle.general_relativity.coordinates.velocity.x)
                               * (2.*(star_angular_momentum.z+particle_angular_momentum.z) + mass_spin_factor_z));

            let element7s = Axes{
                                    x:7.*(star_angular_momentum.x+particle_angular_momentum.x) + 3.*mass_spin_factor_x,
                                    y:7.*(star_angular_momentum.y+particle_angular_momentum.y) + 3.*mass_spin_factor_y,
                                    z:7.*(star_angular_momentum.z+particle_angular_momentum.z) + 3.*mass_spin_factor_z
            };
            let element2_x: f64 = particle.general_relativity.coordinates.velocity.y * element7s.z - particle.general_relativity.coordinates.velocity.z * element7s.y;
            let element2_y: f64 = particle.general_relativity.coordinates.velocity.z * element7s.x - particle.general_relativity.coordinates.velocity.x * element7s.z;
            let element2_z: f64 = particle.general_relativity.coordinates.velocity.x * element7s.y - particle.general_relativity.coordinates.velocity.y * element7s.x;

            let element3s = Axes{
                                    x:3.*(star_angular_momentum.x+particle_angular_momentum.x) + mass_spin_factor_x,
                                    y:3.*(star_angular_momentum.y+particle_angular_momentum.y) + mass_spin_factor_y,
                                    z:3.*(star_angular_momentum.z+particle_angular_momentum.z) + mass_spin_factor_z
            };
            let element3_x: f64 = 3.*particle.general_relativity.parameters.internal.radial_velocity * (particle_normalized_position.y * element3s.z - particle_normalized_position.z * element3s.y);
            let element3_y: f64 = 3.*particle.general_relativity.parameters.internal.radial_velocity * (particle_normalized_position.z * element3s.x - particle_normalized_position.x * element3s.z);
            let element3_z: f64 = 3.*particle.general_relativity.parameters.internal.radial_velocity * (particle_normalized_position.x * element3s.y - particle_normalized_position.y * element3s.x);

            let factor_a = G / SPEED_OF_LIGHT_2;
            let total_general_relativity_spin_orbit_acceleration_x = factor_a * (element1_x - element2_x + element3_x);
            let total_general_relativity_spin_orbit_acceleration_y = factor_a * (element1_y - element2_y + element3_y);
            let total_general_relativity_spin_orbit_acceleration_z = factor_a * (element1_z - element2_z + element3_z);

            sum_total_general_relativity_spin_orbit_acceleration.x += particle.mass/host_particle.mass * total_general_relativity_spin_orbit_acceleration_x;
            sum_total_general_relativity_spin_orbit_acceleration.y += particle.mass/host_particle.mass * total_general_relativity_spin_orbit_acceleration_y;
            sum_total_general_relativity_spin_orbit_acceleration.z += particle.mass/host_particle.mass * total_general_relativity_spin_orbit_acceleration_z;

            particle.general_relativity.parameters.output.acceleration.x += total_general_relativity_spin_orbit_acceleration_x;
            particle.general_relativity.parameters.output.acceleration.y += total_general_relativity_spin_orbit_acceleration_y;
            particle.general_relativity.parameters.output.acceleration.z += total_general_relativity_spin_orbit_acceleration_z;
            //println!("{} {} {}", total_general_relativity_spin_orbit_acceleration_x, total_general_relativity_spin_orbit_acceleration_y, total_general_relativity_spin_orbit_acceleration_z);

            // Kidder 1995, equation 2.4a
            let mu = (host_particle.mass * particle.mass) / star_planet_mass;
            let newtonian_orbital_angular_momentum = Axes{
                                                x:mu * (particle.general_relativity.coordinates.position.y * particle.general_relativity.coordinates.velocity.z - particle.general_relativity.coordinates.position.z * particle.general_relativity.coordinates.velocity.y),
                                                y:mu * (particle.general_relativity.coordinates.position.z * particle.general_relativity.coordinates.velocity.x - particle.general_relativity.coordinates.position.x * particle.general_relativity.coordinates.velocity.z),
                                                z:mu * (particle.general_relativity.coordinates.position.x * particle.general_relativity.coordinates.velocity.y - particle.general_relativity.coordinates.position.y * particle.general_relativity.coordinates.velocity.x)
            };

            let factor_mass = 2. + 3./2. * particle.mass/host_particle.mass;
            let element1_x: f64 = factor_mass 
                * (newtonian_orbital_angular_momentum.y * star_angular_momentum.z - newtonian_orbital_angular_momentum.z * star_angular_momentum.y);
            let element1_y: f64 = factor_mass
                * (newtonian_orbital_angular_momentum.z * star_angular_momentum.x - newtonian_orbital_angular_momentum.x * star_angular_momentum.z);
            let element1_z: f64 = factor_mass
                * (newtonian_orbital_angular_momentum.x * star_angular_momentum.y - newtonian_orbital_angular_momentum.y * star_angular_momentum.x);

            let element2_x: f64 = particle_angular_momentum.y * star_angular_momentum.z - particle_angular_momentum.z * star_angular_momentum.y;
            let element2_y :f64 = particle_angular_momentum.z * star_angular_momentum.x - particle_angular_momentum.x * star_angular_momentum.z;
            let element2_z: f64 = particle_angular_momentum.x * star_angular_momentum.y - particle_angular_momentum.y * star_angular_momentum.x;

            let scalar_product_particle_normalized_position_with_particle_angular_momentum = 
                particle_normalized_position.x * particle_angular_momentum.x 
                + particle_normalized_position.y * particle_angular_momentum.y 
                + particle_normalized_position.z * particle_angular_momentum.z;
            let element3_x: f64 = 3. * scalar_product_particle_normalized_position_with_particle_angular_momentum 
                * (particle_normalized_position.y * star_angular_momentum.z - particle_normalized_position.z * star_angular_momentum.y);
            let element3_y :f64 = 3. * scalar_product_particle_normalized_position_with_particle_angular_momentum 
                * (particle_normalized_position.z * star_angular_momentum.x - particle_normalized_position.x * star_angular_momentum.z);
            let element3_z: f64 = 3. * scalar_product_particle_normalized_position_with_particle_angular_momentum 
                * (particle_normalized_position.x * star_angular_momentum.y - particle_normalized_position.y * star_angular_momentum.x);
            
            host_particle.general_relativity.parameters.output.dangular_momentum_dt.x += factor_a * (element1_x - element2_x + element3_x);
            host_particle.general_relativity.parameters.output.dangular_momentum_dt.y += factor_a * (element1_y - element2_y + element3_y);
            host_particle.general_relativity.parameters.output.dangular_momentum_dt.z += factor_a * (element1_z - element2_z + element3_z);
            //println!("{} {} {}", factor_a * (element1_x - element2_x + element3_x), factor_a * (element1_y - element2_y + element3_y), factor_a * (element1_z - element2_z + element3_z));
            
            // Kidder 1995, equation 2.4b
            let factor_mass = 2. + 3./2. * host_particle.mass/particle.mass;
            let element1_x: f64 = factor_mass 
                * (newtonian_orbital_angular_momentum.y * particle_angular_momentum.z - newtonian_orbital_angular_momentum.z * particle_angular_momentum.y);
            let element1_y: f64 = factor_mass
                * (newtonian_orbital_angular_momentum.z * particle_angular_momentum.x - newtonian_orbital_angular_momentum.x * particle_angular_momentum.z);
            let element1_z: f64 = factor_mass
                * (newtonian_orbital_angular_momentum.x * particle_angular_momentum.y - newtonian_orbital_angular_momentum.y * particle_angular_momentum.x);

            let element2_x: f64 = star_angular_momentum.y * particle_angular_momentum.z - star_angular_momentum.z * particle_angular_momentum.y;
            let element2_y :f64 = star_angular_momentum.z * particle_angular_momentum.x - star_angular_momentum.x * particle_angular_momentum.z;
            let element2_z: f64 = star_angular_momentum.x * particle_angular_momentum.y - star_angular_momentum.y * particle_angular_momentum.x;

            let scalar_product_particle_normalized_position_with_star_angular_momentum = 
                particle_normalized_position.x * star_angular_momentum.x 
                + particle_normalized_position.y * star_angular_momentum.y 
                + particle_normalized_position.z * star_angular_momentum.z;
            let element3_x: f64 = 3. * scalar_product_particle_normalized_position_with_star_angular_momentum 
                * (particle_normalized_position.y*particle_angular_momentum.z - particle_normalized_position.z*particle_angular_momentum.y);
            let element3_y :f64 = 3. * scalar_product_particle_normalized_position_with_star_angular_momentum 
                * (particle_normalized_position.z*particle_angular_momentum.x - particle_normalized_position.x*particle_angular_momentum.z);
            let element3_z: f64 = 3. * scalar_product_particle_normalized_position_with_star_angular_momentum
                * (particle_normalized_position.x*particle_angular_momentum.y - particle_normalized_position.y*particle_angular_momentum.x);
            
            particle.general_relativity.parameters.output.dangular_momentum_dt.x = factor_a * (element1_x - element2_x + element3_x);
            particle.general_relativity.parameters.output.dangular_momentum_dt.y = factor_a * (element1_y - element2_y + element3_y);
            particle.general_relativity.parameters.output.dangular_momentum_dt.z = factor_a * (element1_z - element2_z + element3_z);
        }
    }
    //for particle in particles.iter_mut() {
        //particle.general_relativity.parameters.output.acceleration.x += sum_total_general_relativity_spin_orbit_acceleration.x;
        //particle.general_relativity.parameters.output.acceleration.y += sum_total_general_relativity_spin_orbit_acceleration.y;
        //particle.general_relativity.parameters.output.acceleration.z += sum_total_general_relativity_spin_orbit_acceleration.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    host_particle.general_relativity.parameters.output.acceleration.x += -1.0 * sum_total_general_relativity_spin_orbit_acceleration.x;
    host_particle.general_relativity.parameters.output.acceleration.y += -1.0 * sum_total_general_relativity_spin_orbit_acceleration.y;
    host_particle.general_relativity.parameters.output.acceleration.z += -1.0 * sum_total_general_relativity_spin_orbit_acceleration.z;
}

////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------
// [start] General Relativity based on REBOUNDx gr.c
pub fn calculate_anderson1975_general_relativity_acceleration(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], ignored_gravity_terms: IgnoreGravityTerms) {
    // Calculate Newtonian accelerations in the current setup and considering all particles
    let mut host_particle = host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let (host_newtonian_inertial_accelerations, newtonian_inertial_accelerations) = get_anderson1975_newhall1983_newtonian_inertial_accelerations(&mut host_particle, &mut particles, &mut more_particles, ignored_gravity_terms);

    // Transform to Jacobi coordinates
    let (jacobi_star_mass, _jacobi_star_position, _jacobi_star_velocity, _jacobi_star_acceleration, jacobi_particles_positions, jacobi_particles_velocities, mut jacobi_particles_accelerations) = anderson1975_general_relativity_inertial_to_jacobi_posvelacc(&host_particle, &particles, &more_particles, host_newtonian_inertial_accelerations, newtonian_inertial_accelerations);

    let n_particles = particles.len() + more_particles.len();
    let mu = host_particle.mass_g;
    for (((jacobi_particle_acceleration, jacobi_particle_velocity), jacobi_particle_position), particle) in jacobi_particles_accelerations[..n_particles].iter_mut()
                                                                                                    .zip(jacobi_particles_velocities[..n_particles].iter())
                                                                                                    .zip(jacobi_particles_positions[..n_particles].iter())
                                                                                                    .zip(particles.iter().chain(more_particles.iter())){
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let mut vi = Axes{x: jacobi_particle_velocity.x, y: jacobi_particle_velocity.y, z: jacobi_particle_velocity.z};
            let mut vi2 = jacobi_particle_velocity.x.powi(2) + jacobi_particle_velocity.y.powi(2) + jacobi_particle_velocity.z.powi(2);
            let ri = (jacobi_particle_position.x.powi(2) + jacobi_particle_position.y.powi(2) + jacobi_particle_position.z.powi(2)).sqrt();
            let mut factor_a = (0.5*vi2 + 3.*mu/ri)/SPEED_OF_LIGHT_2;
            let mut old_v = Axes{x:0., y:0., z:0.};

            let max_iterations = 10;
            for q in 0..max_iterations {
                old_v.x = vi.x;
                old_v.y = vi.y;
                old_v.z = vi.z;
                vi.x = jacobi_particle_velocity.x/(1.-factor_a);
                vi.y = jacobi_particle_velocity.y/(1.-factor_a);
                vi.z = jacobi_particle_velocity.z/(1.-factor_a);
                vi2 =vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
                factor_a = (0.5*vi2 + 3.*mu/ri)/SPEED_OF_LIGHT_2;
                let dvx = vi.x - old_v.x;
                let dvy = vi.y - old_v.y;
                let dvz = vi.z - old_v.z;
                if (dvx*dvx + dvy*dvy + dvz*dvz)/vi2 < DBL_EPSILON_2 {
                    break;
                } else if q == max_iterations {
                    println!("[WARNING {} UTC] {} iterations in general relativity failed to converge. This is typically because the perturbation is too strong for the current implementation.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), max_iterations);
                }
            }

            let factor_b = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/SPEED_OF_LIGHT_2;
            let rdotrdot = jacobi_particle_position.x*jacobi_particle_velocity.x
                            + jacobi_particle_position.y*jacobi_particle_velocity.y
                            + jacobi_particle_position.z*jacobi_particle_velocity.z;
            let vidot = Axes{x: jacobi_particle_acceleration.x + factor_b*jacobi_particle_position.x,
                                y: jacobi_particle_acceleration.y + factor_b*jacobi_particle_position.y,
                                z: jacobi_particle_acceleration.z + factor_b*jacobi_particle_position.z};
            let vdotvdot = vi.x*vidot.x + vi.y*vidot.y + vi.z*vidot.z;
            let factor_d = (vdotvdot - 3.*mu/(ri*ri*ri)*rdotrdot)/SPEED_OF_LIGHT_2;
            jacobi_particle_acceleration.x = factor_b*(1.-factor_a)*jacobi_particle_position.x - factor_a*jacobi_particle_acceleration.x - factor_d*vi.x;
            jacobi_particle_acceleration.y = factor_b*(1.-factor_a)*jacobi_particle_position.y - factor_a*jacobi_particle_acceleration.y - factor_d*vi.y;
            jacobi_particle_acceleration.z = factor_b*(1.-factor_a)*jacobi_particle_position.z - factor_a*jacobi_particle_acceleration.z - factor_d*vi.z;
        }
    }


    let jacobi_star_acceleration = Axes{x:0., y:0., z:0.};
    let (star_acceleration, particles_accelerations) = anderson1975_general_relativity_jacobi_to_inertial_acc(&mut particles, &mut more_particles, jacobi_star_mass, jacobi_star_acceleration, jacobi_particles_accelerations);


    // This algorithm computes general_relativity.parameters.output.acceleration in the inertial frame,
    // which is the same coordinate system that is expressed all the rest of additional
    // effects
    for (particle, particle_acceleration) in particles.iter_mut().chain(more_particles.iter_mut()).zip(particles_accelerations.iter()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            particle.general_relativity.parameters.output.acceleration.x = particle_acceleration.x;
            particle.general_relativity.parameters.output.acceleration.y = particle_acceleration.y;
            particle.general_relativity.parameters.output.acceleration.z = particle_acceleration.z;
        }
    }
    host_particle.general_relativity.parameters.output.acceleration.x = star_acceleration.x;
    host_particle.general_relativity.parameters.output.acceleration.y = star_acceleration.y;
    host_particle.general_relativity.parameters.output.acceleration.z = star_acceleration.z;
}

fn anderson1975_general_relativity_inertial_to_jacobi_posvelacc(host_particle: &Particle, particles: &[Particle], more_particles: &[Particle], host_newtonian_inertial_accelerations: Axes, newtonian_inertial_accelerations: [Axes; MAX_PARTICLES-1]) -> (f64, Axes, Axes, Axes, [Axes; MAX_PARTICLES-1], [Axes; MAX_PARTICLES-1], [Axes; MAX_PARTICLES-1]) {
    let jacobi_star_mass;
    let mut jacobi_star_position = Axes{x:0., y:0., z:0. };
    let mut jacobi_star_velocity = Axes{x:0., y:0., z:0. };
    let mut jacobi_star_acceleration = Axes{x:0., y:0., z:0. };
    let mut jacobi_particles_positions = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];
    let mut jacobi_particles_velocities = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];
    let mut jacobi_particles_accelerations = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];

    let m0 = host_particle.mass;
    let mut eta = m0;
    let mut s_x = eta * host_particle.inertial_position.x;
    let mut s_y = eta * host_particle.inertial_position.y;
    let mut s_z = eta * host_particle.inertial_position.z;
    let mut s_vx = eta * host_particle.inertial_velocity.x;
    let mut s_vy = eta * host_particle.inertial_velocity.y;
    let mut s_vz = eta * host_particle.inertial_velocity.z;
    let mut s_ax = eta * host_newtonian_inertial_accelerations.x;
    let mut s_ay = eta * host_newtonian_inertial_accelerations.y;
    let mut s_az = eta * host_newtonian_inertial_accelerations.z;
    for ((((particle, particle_newtonian_inertial_accelerations), jacobi_particle_position), jacobi_particle_velocity), jacobi_particle_acceleration) in 
            particles.iter().chain(more_particles.iter()) // zip will pick the lowest common number of elements
                .zip(newtonian_inertial_accelerations.iter())
                .zip(jacobi_particles_positions.iter_mut())
                .zip(jacobi_particles_velocities.iter_mut())
                .zip(jacobi_particles_accelerations.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let ei = 1./eta;
            eta += particle.mass;
            let pme = eta*ei;
            jacobi_particle_position.x = particle.inertial_position.x - s_x*ei;
            jacobi_particle_position.y = particle.inertial_position.y - s_y*ei;
            jacobi_particle_position.z = particle.inertial_position.z - s_z*ei;
            jacobi_particle_velocity.x = particle.inertial_velocity.x - s_vx*ei;
            jacobi_particle_velocity.y = particle.inertial_velocity.y - s_vy*ei;
            jacobi_particle_velocity.z = particle.inertial_velocity.z - s_vz*ei;
            jacobi_particle_acceleration.x = particle_newtonian_inertial_accelerations.x - s_ax*ei;
            jacobi_particle_acceleration.y = particle_newtonian_inertial_accelerations.y - s_ay*ei;
            jacobi_particle_acceleration.z = particle_newtonian_inertial_accelerations.z - s_az*ei;
            s_x  = s_x  * pme + particle.mass*jacobi_particle_position.x ;
            s_y  = s_y  * pme + particle.mass*jacobi_particle_position.y ;
            s_z  = s_z  * pme + particle.mass*jacobi_particle_position.z ;
            s_vx = s_vx * pme + particle.mass*jacobi_particle_velocity.x;
            s_vy = s_vy * pme + particle.mass*jacobi_particle_velocity.y;
            s_vz = s_vz * pme + particle.mass*jacobi_particle_velocity.z;
            s_ax = s_ax * pme + particle.mass*jacobi_particle_acceleration.x;
            s_ay = s_ay * pme + particle.mass*jacobi_particle_acceleration.y;
            s_az = s_az * pme + particle.mass*jacobi_particle_acceleration.z;
        }
    }
    let mtot = eta;
    let mtot_i = 1./mtot;
    jacobi_star_mass = mtot;
    jacobi_star_position.x = s_x * mtot_i;
    jacobi_star_position.y = s_y * mtot_i;
    jacobi_star_position.z = s_z * mtot_i;
    jacobi_star_velocity.x = s_vx * mtot_i;
    jacobi_star_velocity.y = s_vy * mtot_i;
    jacobi_star_velocity.z = s_vz * mtot_i;
    jacobi_star_acceleration.x = s_ax * mtot_i;
    jacobi_star_acceleration.y = s_ay * mtot_i;
    jacobi_star_acceleration.z = s_az * mtot_i;
    return(jacobi_star_mass, jacobi_star_position, jacobi_star_velocity, jacobi_star_acceleration, jacobi_particles_positions, jacobi_particles_velocities, jacobi_particles_accelerations)
}

fn anderson1975_general_relativity_jacobi_to_inertial_acc(particles: &mut [Particle], more_particles: &mut [Particle], jacobi_star_mass: f64, jacobi_star_acceleration: Axes, jacobi_particles_accelerations: [Axes; MAX_PARTICLES-1]) -> (Axes, [Axes; MAX_PARTICLES-1]) {
    let mut star_acceleration = Axes{x:0., y:0., z:0. };
    let mut particles_accelerations = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1];
    let n_particles = particles.len() + more_particles.len();

    let mut eta = jacobi_star_mass;
    let mut s_ax = eta * jacobi_star_acceleration.x;
    let mut s_ay = eta * jacobi_star_acceleration.y;
    let mut s_az = eta * jacobi_star_acceleration.z;
    for ((particle, particle_acceleration), jacobi_particle_acceleration) in particles.iter().chain(more_particles.iter()).rev()
                                                                                .zip(particles_accelerations[..n_particles].iter_mut().rev())
                                                                                .zip(jacobi_particles_accelerations[..n_particles].iter().rev()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let ei = 1./eta;
            s_ax = (s_ax - particle.mass*jacobi_particle_acceleration.x) * ei;
            s_ay = (s_ay - particle.mass*jacobi_particle_acceleration.y) * ei;
            s_az = (s_az - particle.mass*jacobi_particle_acceleration.z) * ei;
            particle_acceleration.x = jacobi_particle_acceleration.x + s_ax;
            particle_acceleration.y = jacobi_particle_acceleration.y + s_ay;
            particle_acceleration.z = jacobi_particle_acceleration.z + s_az;
            eta -= particle.mass;
            s_ax *= eta;
            s_ay *= eta;
            s_az *= eta;
        }
    }
    let mtot = eta;
    let mtot_i = 1./mtot;
    star_acceleration.x = s_ax * mtot_i;
    star_acceleration.y = s_ay * mtot_i;
    star_acceleration.z = s_az * mtot_i;
    return(star_acceleration, particles_accelerations)
}
// [end] General Relativity based on REBOUNDx gr.c
//--------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

fn get_anderson1975_newhall1983_newtonian_inertial_accelerations(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], ignored_gravity_terms: IgnoreGravityTerms) -> (Axes, [Axes; MAX_PARTICLES-1]) {
    let mut host_newtonian_inertial_accelerations = Axes{x:0., y:0., z:0. };
    let mut newtonian_inertial_accelerations = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES-1]; 
    host_newtonian_inertial_accelerations.x = host_particle.inertial_acceleration.x;
    host_newtonian_inertial_accelerations.y = host_particle.inertial_acceleration.y;
    host_newtonian_inertial_accelerations.z = host_particle.inertial_acceleration.z;
    for (newtonian_acceleration, particle) in newtonian_inertial_accelerations.iter_mut().zip(particles.iter_mut().chain(more_particles.iter_mut())) {
        newtonian_acceleration.x = particle.inertial_acceleration.x;
        newtonian_acceleration.y = particle.inertial_acceleration.y;
        newtonian_acceleration.z = particle.inertial_acceleration.z;
    }
    // If some terms where ignored by the integrator, they should be added
    if ignored_gravity_terms == IgnoreGravityTerms::WHFastOne || ignored_gravity_terms == IgnoreGravityTerms::WHFastTwo {
        let n_particles;
        if ignored_gravity_terms == IgnoreGravityTerms::WHFastOne {
            n_particles = 2 - 1; // Only host - next particle interaction, and the host is already out of the particles vector
        } else {
            n_particles = particles.len() + more_particles.len();
        }
        for (newtonian_acceleration, particle) in newtonian_inertial_accelerations[..n_particles].iter_mut().zip(particles.iter_mut().chain(more_particles.iter_mut())) {
            let dx = host_particle.inertial_position.x - particle.general_relativity.coordinates.position.x;
            let dy = host_particle.inertial_position.y - particle.general_relativity.coordinates.position.y;
            let dz = host_particle.inertial_position.z - particle.general_relativity.coordinates.position.z;
            let r2 = dx.powi(2) + dy.powi(2) + dz.powi(2);
            let r = r2.sqrt();
            let prefac = G/(r2*r);
            let prefac_mass_star = prefac*host_particle.mass;
            let prefac_mass_particle = prefac*particle.mass;
            host_newtonian_inertial_accelerations.x -= prefac_mass_particle*dx;
            host_newtonian_inertial_accelerations.y -= prefac_mass_particle*dy;
            host_newtonian_inertial_accelerations.z -= prefac_mass_particle*dz;
            newtonian_acceleration.x += prefac_mass_star*dx;
            newtonian_acceleration.y += prefac_mass_star*dy;
            newtonian_acceleration.z += prefac_mass_star*dz;
        }
    }
    return (host_newtonian_inertial_accelerations, newtonian_inertial_accelerations);
}

////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------
// [start] General Relativity FULL based on REBOUNDx gr.c
pub fn calculate_newhall1983_general_relativity_acceleration(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], ignored_gravity_terms: IgnoreGravityTerms) {
    // host_particle is separated from particles for homogeneity with the rest of methods,
    // but this implementation of General Relativity does not uses a host, all
    // interactions between all the bodies are computed
    let mut host_a_const = [Axes{x:0., y:0., z:0. }; 1]; // array that stores the value of the constant term
    let mut host_a_new = [Axes{x:0., y:0., z:0. }; 1]; // stores the newly calculated term
    let mut host_rs = [[0.; MAX_PARTICLES]; 1];
    let mut host_drs = [[Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; 1];
    let mut a_const = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; // array that stores the value of the constant term
    let mut a_new = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; // stores the newly calculated term
    let mut rs = [[0.; MAX_PARTICLES]; MAX_PARTICLES];
    let mut drs = [[Axes{x:0., y:0., z:0. }; MAX_PARTICLES]; MAX_PARTICLES];
    
    let mut host_particle = host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let (host_newtonian_inertial_accelerations, newtonian_inertial_accelerations) = get_anderson1975_newhall1983_newtonian_inertial_accelerations(&mut host_particle, &mut particles, &mut more_particles, ignored_gravity_terms);

    for (i, ((particle_i, drs_i), rs_i)) in 
                                    iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                    .zip(host_drs.iter_mut().chain(drs.iter_mut()))
                                    .zip(host_rs.iter_mut().chain(rs.iter_mut()))
                                    .enumerate() {
        // compute distances
        for (j, (particle_j, drs_i_j)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                        .zip(drs_i.iter_mut())
                                        .enumerate() {
            if j != i{
                if particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled || particle_j.general_relativity.effect != GeneralRelativityEffect::Disabled {
                    drs_i_j.x = particle_i.inertial_position.x - particle_j.inertial_position.x;
                    drs_i_j.y = particle_i.inertial_position.y - particle_j.inertial_position.y;
                    drs_i_j.z = particle_i.inertial_position.z - particle_j.inertial_position.z;
                    rs_i[j] = (drs_i_j.x.powi(2) + drs_i_j.y.powi(2) + drs_i_j.z.powi(2)).sqrt();
                    //println!("i j: {} {} {:e}", i, j, rs_i[j]);
                }
            }
        }
    }

    for (i, (((particle_i, a_const_i), drs_i), rs_i)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                            .zip(host_a_const.iter_mut().chain(a_const.iter_mut()))
                                            .zip(host_drs.iter_mut().chain(drs.iter_mut()))
                                            .zip(host_rs.iter().chain(rs.iter()))
                                            .enumerate() {

        // then compute the constant terms:
        let mut a_constx = 0.;
        let mut a_consty = 0.;
        let mut a_constz = 0.;
        // 1st constant part
        for (j, (particle_j, drs_i_j)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                        .zip(drs_i.iter())
                                        .enumerate() {
            if j != i {
                if particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled || particle_j.general_relativity.effect != GeneralRelativityEffect::Disabled {
                    let dxij = drs_i_j.x;
                    let dyij = drs_i_j.y;
                    let dzij = drs_i_j.z;
                    let rij2 = rs_i[j].powi(2);
                    let rij3 = rij2*rs_i[j];

                    let mut a1 = 0.;
                    for (k, (particle_k, rs_i_k)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                                    .zip(rs_i.iter())
                                                        .enumerate() {
                        if k != i {
                            a1 += (4./(SPEED_OF_LIGHT_2)) * G*particle_k.mass/rs_i_k;
                        }
                    }

                    let mut a2 = 0.;
                    for (l, (particle_l, rs_l)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                                    .zip(host_rs.iter().chain(rs.iter()))
                                                        .enumerate() {
                        if l != j {
                            a2 += (1./(SPEED_OF_LIGHT_2)) * G*particle_l.mass/rs_l[j];
                        }
                    }

                    let vi2= particle_i.inertial_velocity.x.powi(2) + particle_i.inertial_velocity.y.powi(2) + particle_i.inertial_velocity.z.powi(2);
                    let a3 = -vi2/SPEED_OF_LIGHT_2;
                    
                    let vj2 = particle_j.inertial_velocity.x.powi(2) + particle_j.inertial_velocity.y.powi(2) + particle_j.inertial_velocity.z.powi(2);
                    let a4 = -2.*vj2/SPEED_OF_LIGHT_2;

                    let a5 = (4./SPEED_OF_LIGHT_2) * (particle_i.inertial_velocity.x*particle_j.inertial_velocity.x + particle_i.inertial_velocity.y*particle_j.inertial_velocity.y + particle_i.inertial_velocity.z*particle_j.inertial_velocity.z);

                    let a6_0 = dxij*particle_j.inertial_velocity.x + dyij*particle_j.inertial_velocity.y + dzij*particle_j.inertial_velocity.z;
                    let a6 = (3./(2.*SPEED_OF_LIGHT_2)) * a6_0.powi(2)/rij2;

                    let factor1 = a1 + a2 + a3 + a4 + a5 + a6;
                    //println!("factors {:e} {:e} {:e} {:e} {:e} {:e}", a1, a2, a3, a4, a5, a6);
                    a_constx += G*particle_j.mass*dxij*factor1/rij3;
                    a_consty += G*particle_j.mass*dyij*factor1/rij3;
                    a_constz += G*particle_j.mass*dzij*factor1/rij3;

                    // 2nd constant part
                    let dvxij = particle_i.inertial_velocity.x - particle_j.inertial_velocity.x; 
                    let dvyij = particle_i.inertial_velocity.y - particle_j.inertial_velocity.y; 
                    let dvzij = particle_i.inertial_velocity.z - particle_j.inertial_velocity.z; 

                    let factor2 = dxij*(4.*particle_i.inertial_velocity.x - 3.*particle_j.inertial_velocity.x) + dyij*(4.*particle_i.inertial_velocity.y -3.*particle_j.inertial_velocity.y) + dzij*(4.*particle_i.inertial_velocity.z - 3.*particle_j.inertial_velocity.z);

                    a_constx += G*particle_j.mass*factor2*dvxij/rij3/SPEED_OF_LIGHT_2;
                    a_consty += G*particle_j.mass*factor2*dvyij/rij3/SPEED_OF_LIGHT_2;
                    a_constz += G*particle_j.mass*factor2*dvzij/rij3/SPEED_OF_LIGHT_2;
                }
            }
        }
        a_const_i.x = a_constx;
        a_const_i.y = a_consty;
        a_const_i.z = a_constz;
        //println!("a_const_i {:?}", a_const_i);
    }


    let n_particles = particles.len() + more_particles.len();
    let dev_limit = 1.0e-30;
    let max_iterations = 10;
    // Now running the substitution again and again through the loop below
    for k in 0..max_iterations {
        let host_a_old = host_a_new.clone();
        let a_old = a_new.clone();
        // now add on the non-constant term
        for (i, ((((particle_i, a_new_i), drs_i), rs_i), a_const_i)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                        .zip(host_a_new.iter_mut().chain(a_new.iter_mut()))
                                        .zip(host_drs.iter_mut().chain(drs.iter_mut()))
                                        .zip(host_rs.iter_mut().chain(rs.iter_mut()))
                                        .zip(host_a_const.iter().chain(a_const.iter()))
                                        .enumerate() {
            let mut non_constx = 0.;
            let mut non_consty = 0.;
            let mut non_constz = 0.;
            for (j, ((((particle_j, a_old_j), a_newton_j), drs_i_j), rs_i_j)) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                                                                        .zip(host_a_old.iter().chain(a_old.iter()))
                                                                                        .zip(iter::once(&host_newtonian_inertial_accelerations).chain(newtonian_inertial_accelerations.iter()))
                                                                                        .zip(drs_i.iter())
                                                                                        .zip(rs_i.iter())
                                                                                        .enumerate() {
                if j != i {
                    if particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled || particle_j.general_relativity.effect != GeneralRelativityEffect::Disabled {
                        let dxij = drs_i_j.x;
                        let dyij = drs_i_j.y;
                        let dzij = drs_i_j.z;
                        let rij = rs_i_j;
                        let rij2 = rij.powi(2);
                        let rij3 = rij2*rij;
                        non_constx += (G*particle_j.mass*dxij/rij3)*(dxij*(a_newton_j.x+a_old_j.x)+dyij*(a_newton_j.y+a_old_j.y)+
                                    dzij*(a_newton_j.z+a_old_j.z))/(2.*SPEED_OF_LIGHT_2) + (7./(2.*SPEED_OF_LIGHT_2))*G*particle_j.mass*(a_newton_j.x+a_old_j.x)/rij;
                        non_consty += (G*particle_j.mass*dyij/rij3)*(dxij*(a_newton_j.x+a_old_j.x)+dyij*(a_newton_j.y+a_old_j.y)+
                                    dzij*(a_newton_j.z+a_old_j.z))/(2.*SPEED_OF_LIGHT_2) + (7./(2.*SPEED_OF_LIGHT_2))*G*particle_j.mass*(a_newton_j.y+a_old_j.y)/rij;
                        non_constz += (G*particle_j.mass*dzij/rij3)*(dxij*(a_newton_j.x+a_old_j.x)+dyij*(a_newton_j.y+a_old_j.y)+
                                    dzij*(a_newton_j.z+a_old_j.z))/(2.*SPEED_OF_LIGHT_2) + (7./(2.*SPEED_OF_LIGHT_2))*G*particle_j.mass*(a_newton_j.z+a_old_j.z)/rij;
                    }
                }
            }
            a_new_i.x = a_const_i.x + non_constx;
            a_new_i.y = a_const_i.y + non_consty;
            a_new_i.z = a_const_i.z + non_constz;
            //println!("non_constx {:?}", non_constx);
            //println!("non_consty {:?}", non_consty);
            //println!("non_constz {:?}", non_constz);
        }
        
        // break out loop if a_new is converging
        let mut maxdev = 0.;
        let mut dx = 0.;
        let mut dy = 0.;
        let mut dz = 0.;
        for ((particle_i, a_new_i), a_old_i) in iter::once(&*host_particle).chain(particles.iter()).chain(more_particles.iter()) // zip will pick the lowest common number of elements
                                .zip(host_a_new.iter_mut().chain(a_new[..n_particles].iter_mut()))
                                .zip(host_a_old.iter().chain(a_old.iter())) {
            if particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled {
                if a_new_i.x.abs() < dev_limit {
                    dx = (a_new_i.x - a_old_i.x).abs() / a_new_i.x;
                }
                if a_new_i.y.abs() < dev_limit {
                    dy = (a_new_i.y - a_old_i.y).abs() / a_new_i.y;
                }
                if a_new_i.z.abs() < dev_limit {
                    dz = (a_new_i.z - a_old_i.z).abs() / a_new_i.z;
                }
                if dx > maxdev { maxdev = dx; }
                if dy > maxdev { maxdev = dy; }
                if dz > maxdev { maxdev = dz; }
            }
        }

        if maxdev < dev_limit {
            break;
        } else if k == max_iterations {
            println!("[WARNING {} UTC] {} iterations in general relativity failed to converge.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), max_iterations);
        }

    }
    
    //// update acceleration in particles
    // This algorithm computes general_relativity.parameters.output.acceleration in the inertial frame,
    // which is the same coordinate system that is expressed all the rest of additional
    // effects
    for (particle, a_new_particle) in particles.iter_mut().chain(more_particles.iter_mut())
                        .zip(a_new.iter()){
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            particle.general_relativity.parameters.output.acceleration.x = a_new_particle.x;
            particle.general_relativity.parameters.output.acceleration.y = a_new_particle.y;
            particle.general_relativity.parameters.output.acceleration.z = a_new_particle.z;
        }
    }
    host_particle.general_relativity.parameters.output.acceleration.x = host_a_new[0].x;
    host_particle.general_relativity.parameters.output.acceleration.y = host_a_new[0].y;
    host_particle.general_relativity.parameters.output.acceleration.z = host_a_new[0].z;

}
// [end] General Relativity FULL based on REBOUNDx gr.c
//--------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////
