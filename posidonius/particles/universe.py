import os
import datetime
from axes import Axes
from posidonius.integrator import WHFast, Ias15, LeapFrog
from posidonius.constants import *
from posidonius.effects.evolution import NonEvolving, Leconte2011, Baraffe2015, Baraffe1998, LeconteChabrier2013, BolmontMathis2016, GalletBolmont2017
from posidonius.tools import calculate_spin, mass_radius_relation, get_center_of_mass_of_pair
import posidonius.effects as effects
from posidonius.particles.particle import Particle, DummyParticle

class ConsiderEffects(object):
    def __init__(self, input_properties):
        self._data = {
            "tides": False,
            "rotational_flattening": False,
            "general_relativity": False,
            "disk": False,
            "wind": False,
            "evolution": False,
        }
        # Update default values, ignore non-recognised keys
        for key, value in input_properties.iteritems():
            if key in self._data:
                self._data[key] = value

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()


class Universe(object):
    def __init__(self, initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects):
        self._time_step = time_step
        self._recovery_snapshot_period = recovery_snapshot_period
        self._historic_snapshot_period = historic_snapshot_period
        self._data = {
            "consider_effects": consider_effects.get(),
            "general_relativity_implementation": "Disabled",
            "hosts": {
                "index": {
                    "most_massive": MAX_PARTICLES+1,
                    "tides": MAX_PARTICLES+1,
                    "rotational_flattening": MAX_PARTICLES+1,
                    "general_relativity": MAX_PARTICLES+1,
                    "disk": MAX_PARTICLES+1,
                },
                "most_massive": {
                    "all": False,
                    "tides": False,
                    "rotational_flattening": False,
                    "general_relativity": False,
                    "disk": False,
                }
            },
            "initial_time": float(initial_time),
            "n_particles": 0,
            "particles": [],
            "particles_evolvers": [],
            "star_planet_dependent_dissipation_factors": {},
            "temporary_copied_particle_positions": [],
            "temporary_copied_particles_masses": [],
            "temporary_copied_particles_radiuses": [],
            "time_limit": float(time_limit),
        }


    def add_dummy_particle(self):
        self.add_particle(DummyParticle())
        self._data['n_particles'] -= 1 # Compensate the addition from the previous add_particle call

    def add_particle(self, particle):
        if self._data['n_particles'] > MAX_PARTICLES:
            raise Exception("Maximum number of particles reached: {}".format(MAX_PARTICLES))

        if effects.tides.CentralBody in particle.effects():
            if self._data['hosts']['index']['tides'] == MAX_PARTICLES+1:
                self._data['hosts']['index']['tides'] = self._data['n_particles']
            else:
                raise Exception("There can only be one central body for tidal effects!")

        if effects.rotational_flattening.CentralBody in particle.effects():
            if self._data['hosts']['index']['rotational_flattening'] == MAX_PARTICLES+1:
                self._data['hosts']['index']['rotational_flattening'] = self._data['n_particles']
            else:
                raise Exception("There can only be one central body for rotational flattening effects!")

        if effects.general_relativity.CentralBody in particle.effects():
            if self._data['hosts']['index']['general_relativity'] == MAX_PARTICLES+1:
                self._data['hosts']['index']['general_relativity'] = self._data['n_particles']
            else:
                raise Exception("There can only be one central body for rotational flattening effects!")
            self._data['general_relativity_implementation'] = particle.general_relativity_implementation()

        if effects.disk.CentralBody in particle.effects():
            if self._data['hosts']['index']['disk'] == MAX_PARTICLES+1:
                self._data['hosts']['index']['disk'] = self._data['n_particles']
            else:
                raise Exception("Only one body with a disk is allowed!")

        evolver = particle.get_evolver(self._data['initial_time'])
        if len(evolver['time']) > 0 and evolver['time'][0] > 0.:
            raise Exception("Your initial time ({} days | {:.2e} years) is smaller than the minimum allowed age of the star ({} days | {:.2e} years)".format(self._data['initial_time'], self._data['initial_time']/365.25, evolver['time'][0]+self._data['initial_time'], (evolver['time'][0]+self._data['initial_time'])/365.25));
        if len(evolver['time']) > 0 and evolver['time'][-1] < self._data['initial_time']:
            raise Exception("Your time limit ({} days | {:.2e} years) is greater than the maximum allowed age of the star ({} days | {:.2e} years)", self._data['initial_time'], self._data['initial_time']/365.25, evolver['time'][0], evolver['time'][0]/365.25)

        self._data['particles'].append(particle.get())
        self._data['particles_evolvers'].append(evolver)
        self._data['temporary_copied_particles_radiuses'].append(0.0)
        self._data['temporary_copied_particles_masses'].append(0.0)
        self._data['temporary_copied_particle_positions'].append(Axes(0.0, 0.0, 0.0).get())
        self._data['n_particles'] += 1

    def add_brown_dwarf(self, mass, dissipation_factor_scale, position, velocity, general_relativity_implementation, evolution_type, wind_k_factor=0., wind_rotation_saturation=0.):
        rotation_period = None
        love_number = None
        if type(evolution_type) == NonEvolving:
            rotation_period = 70.0 # hours
            love_number = 0.307 # BrownDwarf
        elif type(evolution_type) in (Leconte2011, Baraffe2015):
            mass = evolution_type._data[evolution_type.__class__.__name__]

            if mass <= 0.0101 and mass >= 0.0099:
                rotation_period = 8.0
                love_number = 0.3790
            elif mass <= 0.0121 and mass >= 0.0119:
                rotation_period = 13.0
                love_number = 0.3780
            elif mass <= 0.0151 and mass >= 0.0149:
                rotation_period = 19.0
                love_number = 0.3760
            elif mass <= 0.0201 and mass >= 0.0199:
                rotation_period = 24.0
                love_number = 0.3690
            elif mass <= 0.0301 and mass >= 0.0299:
                rotation_period = 30.0
                love_number = 0.3550
            elif mass <= 0.0401 and mass >= 0.0399:
                rotation_period = 36.0
                love_number = 0.3420
            elif mass <= 0.0501 and mass >= 0.0499:
                rotation_period = 41.0
                love_number = 0.3330
            elif mass <= 0.0601 and mass >= 0.0599:
                rotation_period = 47.0
                love_number = 0.3250
            elif mass <= 0.0701 and mass >= 0.0699:
                rotation_period = 53.0
                love_number = 0.3110
            elif mass <= 0.0721 and mass >= 0.0719:
                rotation_period = 58.0
                love_number = 0.3080
            elif mass <= 0.0751 and mass >= 0.0749:
                rotation_period = 64.0
                love_number = 0.3070
            elif mass <= 0.0801 and mass >= 0.0799:
                rotation_period = 70.0
                love_number = 0.3070
            else:
                raise Exception("The evolution type Leconte2011 does not support a mass of {} Msun and Baraffe2015 with higher masses is possible but then it should not be considered as a Brown Dwarf!".format(mass))
        else:
            raise Exception("Evolution type should be Leconte2011, Baraffe2015 or NonEvolving!")

        angular_frequency = TWO_PI/(rotation_period/24.) # days^-1
        inclination = 0.
        obliquity = 0.
        spin = calculate_spin(angular_frequency, inclination, obliquity)

        fluid_love_number = love_number
        # BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        dissipation_factor = 2.006*3.845764e4 # -60+64

        radius_factor = 0.845649342247916
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 1.94e-1 # Brown dwarf

        tides = effects.tides.CentralBody({
            "dissipation_factor_scale": dissipation_factor_scale,
            "dissipation_factor": dissipation_factor,
            "love_number": love_number,
        })
        rotational_flattening = effects.rotational_flattening.CentralBody({"fluid_love_number": fluid_love_number})
        general_relativity = effects.general_relativity.CentralBody(general_relativity_implementation)
        if wind_k_factor == 0:
            wind = effects.wind.Disabled()
        else:
            wind = effects.wind.Interaction({
                "k_factor": wind_k_factor,
                "rotation_saturation": wind_rotation_saturation,
            })
        disk = effects.disk.Disabled()
        self.add_particle(Particle(mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type))

    def add_solar_like(self, mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution_type, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842):
        """
        Wind parametrisation (Bouvier 1997):

        wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
        """
        disk = effects.disk.Disabled()
        self._add_solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution_type, disk, wind_k_factor, wind_rotation_saturation)

    def add_solar_like_with_disk(self, mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution_type, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842):
        """
        Wind parametrisation (Bouvier 1997):

        wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
        """
        disk_surface_density_normalization_gcm = 1000. # g.cm^-2
        disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
        disk_properties = {
            'inner_edge_distance': 0.01,  # AU
            'outer_edge_distance': 100.0, # AU
            'lifetime': 1.0e5 * 365.25e0, # days
            'alpha': 1.0e-2,
            'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/M_SUN) * AU**2, # Msun.AU^-2
            'mean_molecular_weight': 2.4,
        }
        disk = effects.disk.CentralBody(disk_properties)
        self._add_solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution_type, disk, wind_k_factor, wind_rotation_saturation)

    def _add_solar_like_with_disk(self, mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution_type, disk, wind_k_factor, wind_rotation_saturation):
        if type(evolution_type) not in (BolmontMathis2016, Leconte2011, Baraffe2015, GalletBolmont2017, NonEvolving) and not (type(evolution_type) is Baraffe1998 and mass == 1.0):
            raise Exception("Evolution type should be BolmontMathis2016 Leconte2011 Baraffe1998 (mass = 0.10) Baraffe2015 GalletBolmont2017 or NonEvolving!")

        # Typical rotation period: 24 hours
        angular_frequency = TWO_PI/(rotation_period/24.) # days^-1
        inclination = 0.
        obliquity = 0.
        spin = calculate_spin(angular_frequency, inclination, obliquity)

        love_number = 0.03
        fluid_love_number = love_number
        # Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        dissipation_factor = 4.992*3.845764e-2 # -66+64

        radius_factor = 1.
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 5.9e-2 # Sun

        tides = effects.tides.CentralBody({
            "dissipation_factor_scale": dissipation_factor_scale,
            "dissipation_factor": dissipation_factor,
            "love_number": love_number,
        })
        rotational_flattening = effects.rotational_flattening.CentralBody({"fluid_love_number": fluid_love_number})
        general_relativity = effects.general_relativity.CentralBody(general_relativity_implementation)
        if wind_k_factor == 0:
            wind = effects.wind.Disabled()
        else:
            wind = effects.wind.Interaction({
                "k_factor": wind_k_factor,
                "rotation_saturation": wind_rotation_saturation,
            })
        self.add_particle(Particle(mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type))


    def add_m_dwarf(self, mass, dissipation_factor_scale, position, velocity, rotation_period, general_relativity_implementation, evolution_type, wind_k_factor=0., wind_rotation_saturation=0.):
        if type(evolution_type) not in (Baraffe2015, NonEvolving) and not (type(evolution_type) is Baraffe1998 and mass == 0.10):
            raise Exception("Evolution type should be Baraffe2015 Baraffe1998 (mass = 0.10) or NonEvolving!")

        # Typical rotation period: 70 hours
        angular_frequency = TWO_PI/(rotation_period/24.) # days^-1
        inclination = 0.
        obliquity = 0.
        spin = calculate_spin(angular_frequency, inclination, obliquity)

        love_number = 0.307 # M Dwarf
        fluid_love_number = love_number

        # BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        dissipation_factor = 2.006*3.845764e4 # -60+64

        radius_factor = 0.845649342247916
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 2.0e-1 # M-dwarf

        tides = effects.tides.CentralBody({
            "dissipation_factor_scale": dissipation_factor_scale,
            "dissipation_factor": dissipation_factor,
            "love_number": love_number,
        })
        rotational_flattening = effects.rotational_flattening.CentralBody({"fluid_love_number": fluid_love_number})
        general_relativity = effects.general_relativity.CentralBody(general_relativity_implementation)
        if wind_k_factor == 0:
            wind = effects.wind.Disabled()
        else:
            wind = effects.wind.Interaction({
                "k_factor": wind_k_factor,
                "rotation_saturation": wind_rotation_saturation,
            })
        disk = effects.disk.Disabled()
        self.add_particle(Particle(mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type))


    def add_jupiter_like(self, mass, dissipation_factor_scale, position, velocity, spin, evolution_type):
        if type(evolution_type) not in (LeconteChabrier2013, NonEvolving):
            raise Exception("Evolution type should be LeconteChabrier2013 or NonEvolving!")

        # Typical rotation period: 9.8 hours
        love_number = 0.380 # Gas giant
        fluid_love_number = love_number

        radius_factor = 10.9 # Jupiter in R_EARTH
        radius = radius_factor * R_EARTH

        # TODO: What k2pdelta/dissipation_factor is the recommended?
        #k2pdelta = 8.101852e-9 # Gas giant
        k2pdelta = 2.893519e-7 # Gas giant for Jupiter: 2-3d-2 s, here in day (Leconte)
        dissipation_factor = 2. * K2 * k2pdelta/(3. * np.power(radius, 5))
        #dissipation_factor = 2.006*3.845764e4 // Gas giant

        radius_of_gyration_2 = 2.54e-1 # Gas giant

        tides = effects.tides.OrbitingBody({
            "dissipation_factor_scale": dissipation_factor_scale,
            "dissipation_factor": dissipation_factor,
            "love_number": love_number,
        })
        rotational_flattening = effects.rotational_flattening.OrbitingBody({"fluid_love_number": fluid_love_number})
        general_relativity = effects.general_relativity.OrbitingBody()
        wind = effects.wind.Disabled()
        disk = effects.disk.OrbitingBody()
        self.add_particle(Particle(mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type))


    def add_earth_like(self, mass, dissipation_factor_scale, position, velocity, spin, evolution_type):
        if type(evolution_type) not in (NonEvolving,):
            raise Exception("Evolution type should be NonEvolving!")

        # Typical rotation period: 24 hours
        love_number = 0.299 # Earth
        fluid_love_number = 0.9532 # Earth

        # Earth-like => mass-radius relationship from Fortney 2007
        radius_factor = mass_radius_relation(mass, planet_mass_type='AU', planet_percent_rock=0.70)
        radius = radius_factor * R_EARTH
        radius_of_gyration_2 = 3.308e-1 # Earth type planet
        k2pdelta = 2.465278e-3 # Terrestrial planets
        dissipation_factor = 2. * K2 * k2pdelta/(3. * np.power(radius, 5))

        tides = effects.tides.OrbitingBody({
            "dissipation_factor_scale": dissipation_factor_scale,
            "dissipation_factor": dissipation_factor,
            "love_number": love_number,
        })
        rotational_flattening = effects.rotational_flattening.OrbitingBody({"fluid_love_number": fluid_love_number})
        general_relativity = effects.general_relativity.OrbitingBody()
        wind = effects.wind.Disabled()
        disk = effects.disk.OrbitingBody()
        self.add_particle(Particle(mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type))


    def populate_inertial_frame(self):
        # Compute center of mass
        center_of_mass_position = Axes(0., 0., 0.)
        center_of_mass_velocity = Axes(0., 0., 0.)
        center_of_mass_mass = 0.;

        for particle in self._data['particles']:
            center_of_mass_mass = get_center_of_mass_of_pair(center_of_mass_position,
                                                                    center_of_mass_velocity,
                                                                    center_of_mass_mass,
                                                                    particle)

        for particle in self._data['particles']:
            particle['inertial_position']['x'] = particle['heliocentric_position']['x'] - center_of_mass_position.x()
            particle['inertial_position']['y'] = particle['heliocentric_position']['y'] - center_of_mass_position.y()
            particle['inertial_position']['z'] = particle['heliocentric_position']['z'] - center_of_mass_position.z()
            particle['inertial_velocity']['x'] = particle['heliocentric_velocity']['x'] - center_of_mass_velocity.x()
            particle['inertial_velocity']['y'] = particle['heliocentric_velocity']['y'] - center_of_mass_velocity.y()
            particle['inertial_velocity']['z'] = particle['heliocentric_velocity']['z'] - center_of_mass_velocity.z()



    def assign_id(self):
        for i, particle in enumerate(self._data['particles']):
            particle['id'] = i

    def find_most_massive_particle(self):
        max_mass_found = 0.
        for i, particle in enumerate(self._data['particles']):
            if particle['mass'] > max_mass_found:
                max_mass_found = particle['mass']
                self._data['hosts']['index']['most_massive'] = i
        if self._data['hosts']['index']['general_relativity'] != MAX_PARTICLES+1 and self._data['hosts']['index']['general_relativity'] != self._data['hosts']['index']['most_massive']:
                raise Exception("The most massive body should be the central body for general relativity effects!")

    def disable_unnecessary_effects(self):
        if self._data["consider_effects"]["tides"] and self._data['hosts']['index']['tides'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled tides because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["tides"] = False
        if self._data["consider_effects"]["rotational_flattening"] and self._data['hosts']['index']['rotational_flattening'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled rotational flattening because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["rotational_flattening"] = False
        if self._data["consider_effects"]["general_relativity"] and self._data['hosts']['index']['general_relativity'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled general relativity because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["general_relativity"] = False
        if self._data["consider_effects"]["disk"] and self._data['hosts']['index']['disk'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled disk because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["disk"] = False
        found_evolving_body = False
        found_wind = False
        for particle in self._data['particles']:
            if particle["evolution_type"] != "NonEvolving":
                found_evolving_body = True
            if particle["wind"]["effect"] != "Disabled":
                found_wind = True
        if self._data["consider_effects"]["wind"] and not found_wind:
            print("[INFO {} UTC] Disabled wind because no wind was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["wind"] = False
        if self._data["consider_effects"]["evolution"] and not found_evolving_body:
            print("[INFO {} UTC] Disabled evolution because no evolving body was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["evolution"] = False


    def update_hosts(self):
        self._data["hosts"]["most_massive"] = {
            "all": False,
            "tides": self._data["consider_effects"]["tides"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['tides'],
            "rotational_flattening": self._data["consider_effects"]["rotational_flattening"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['rotational_flattening'],
            "general_relativity": self._data["consider_effects"]["general_relativity"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['general_relativity'],
            "disk": self._data["consider_effects"]["disk"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['disk'],
        }
        self._data["hosts"]["most_massive"]["all"] = ((self._data["consider_effects"]["tides"] and self._data["hosts"]["most_massive"]["tides"]) or not self._data["consider_effects"]["tides"]) \
            and ((self._data["consider_effects"]["rotational_flattening"] and self._data["hosts"]["most_massive"]["rotational_flattening"]) or not self._data["consider_effects"]["rotational_flattening"]) \
            and ((self._data["consider_effects"]["general_relativity"] and self._data["hosts"]["most_massive"]["general_relativity"]) or not self._data["consider_effects"]["general_relativity"]) \
            and ((self._data["consider_effects"]["disk"] and self._data["hosts"]["most_massive"]["disk"]) or not self._data["consider_effects"]["disk"])

    def compute_general_relativity_factor(self):
        most_massive_particle_index = self._data['hosts']['index']['most_massive']
        central_body_mass_g = self._data['particles'][most_massive_particle_index]['mass_g']
        for i, particle in enumerate(self._data['particles']):
            if not self._data["consider_effects"]["general_relativity"] \
                    or "CentralBody" in particle['general_relativity']['effect'] \
                    or particle['general_relativity']['effect'] == "Disabled" \
                    or particle['mass'] == 0:
                # Disabled GR or Enabled GR and central body
                particle['general_relativity']['parameters']['internal']['factor'] = 0.0
            else:
                particle['general_relativity']['parameters']['internal']['factor'] = central_body_mass_g*particle['mass_g'] / np.power(central_body_mass_g + particle['mass_g'], 2)

    def get(self):
        self.disable_unnecessary_effects()
        self.populate_inertial_frame()
        self.assign_id()
        self.find_most_massive_particle()
        self.update_hosts()
        self.compute_general_relativity_factor()

        # Add dummy particles to fill the vector
        n_dummy_particles = MAX_PARTICLES - self._data['n_particles']
        for i in range(n_dummy_particles):
            self.add_dummy_particle()

        data = self._data.copy()

        # Reset indices according to enabled effects
        if data["consider_effects"]["rotational_flattening"] and not data["consider_effects"]["tides"]:
            # In practice, rotational_flattening will use tidal_host_particle_index
            # make sure it is the correct index, even if tides are disabled
            data['hosts']['index']['tides'] = data['hosts']['index']['rotational_flattening']
        elif not data["consider_effects"]["tides"]:
            data['hosts']['index']['tides'] = MAX_PARTICLES+1
        if not data["consider_effects"]["rotational_flattening"]:
            data['hosts']['index']['rotational_flattening'] = MAX_PARTICLES+1
        if not data["consider_effects"]["general_relativity"]:
            data['hosts']['index']['general_relativity'] = MAX_PARTICLES+1
        if not data["consider_effects"]["disk"]:
            data['hosts']['index']['disk'] = MAX_PARTICLES+1

        # Forget the dummy particles
        if n_dummy_particles > 0:
            self._data['particles'] = self._data['particles'][:-n_dummy_particles]
            self._data['particles_evolvers'] = self._data['particles_evolvers'][:-n_dummy_particles]
            #self._data['n_particles'] -= n_dummy_particles
        return data

    def write(self, filename, integrator="WHFast", whfast_alternative_coordinates="DemocraticHeliocentric"):
        if integrator.lower() == "whfast":
            universe_integrator = WHFast(whfast_alternative_coordinates, self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        elif integrator.lower() == "ias15":
            universe_integrator = Ias15(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        elif integrator.lower() == "leapfrog":
            universe_integrator = LeapFrog(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        else:
            raise Exception("Unknown integtrator '{}'".format(integrator))
        base_filename = os.path.splitext(filename)[0]
        print("[INFO {} UTC] Start the simulation with:".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        print("posidonius start {} {} {}".format(filename, base_filename+".bin", base_filename+"_history.bin"))


