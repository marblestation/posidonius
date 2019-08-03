import os
import datetime
from axes import Axes
from disk_type import DiskHost, NoDisk, DiskInteraction
from integrator import WHFast, Ias15, LeapFrog
from constants import *
from evolution_type import NonEvolving, Leconte2011, Baraffe2015, Baraffe1998, LeconteChabrier2013, BolmontMathis2016, GalletBolmont2017
from tools import calculate_spin, mass_radius_relation, get_center_of_mass_of_pair

class ConsiderGeneralRelativity(object):
    def __init__(self, variant):
        self._data = {}
        if variant in ("Kidder1995", "Anderson1975", "Newhall1983", "None"):
            self._data = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        return self._data

class Universe(object):
    def __init__(self, initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity):
        self._time_step = time_step
        self._recovery_snapshot_period = recovery_snapshot_period
        self._historic_snapshot_period = historic_snapshot_period
        self._data = {}
        self._data['time_limit'] = float(time_limit)
        self._data['initial_time'] = float(initial_time)
        self._data['consider_tides'] = consider_tides
        self._data['consider_rotational_flattening'] = consider_rotational_flattening
        self._data['consider_disk_interaction'] = consider_disk_interaction
        if consider_general_relativity == True:
            consider_general_relativity = "Kidder1995" # MercuryT
        elif consider_general_relativity == False:
            consider_general_relativity = "None"
        self._data['consider_general_relativity'] = ConsiderGeneralRelativity(consider_general_relativity).get()
        self._data['particles'] = []
        self._data['particles_evolvers'] = []
        self._data['n_particles'] = 0
        self._data['evolving_particles_exist'] = False
        self._data['wind_effects_exist'] = False
        self._data['star_planet_dependent_dissipation_factors'] = {}
        self._data['host_particle_index'] = 0 # Most massive particle
        self._data['tidal_host_particle_index'] = 0 # Particle that is the main one for tidal effects
        self._data['disk_host_particle_index'] = MAX_PARTICLES+1
        self._data['temporary_copied_particles_radiuses'] = []
        self._data['temporary_copied_particles_masses'] = []
        self._data['temporary_copied_particle_velocities'] = []
        self._data['temporary_copied_particle_positions'] = []



    def add_dummy_particle(self):
        mass = 0.0
        radius = 0.0
        dissipation_factor = 0.0
        dissipation_factor_scale = 1.0
        radius_of_gyration_2 = 0.0
        love_number = 0.0
        fluid_love_number = 0.0
        disk = NoDisk()
        position = Axes(0., 0., 0.)
        velocity = Axes(0., 0., 0.)
        spin = Axes(0., 0., 0.)
        evolution_type = NonEvolving()
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk)
        self._data['n_particles'] -= 1 # Compensate the addition from the previous add_particle call

    def add_particle(self, mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk, wind_k_factor=0., wind_rotation_saturation=0.):
        if self._data['n_particles'] > MAX_PARTICLES:
            raise Exception("Maximum number of particles reached: {}".format(MAX_PARTICLES))
        particle = {}
        particle['mass'] = float(mass)
        particle['radius'] = float(radius)
        particle['scaled_dissipation_factor'] = float(dissipation_factor)*float(dissipation_factor_scale)
        particle['dissipation_factor_scale'] = float(dissipation_factor_scale)
        particle['radius_of_gyration_2'] = float(radius_of_gyration_2)
        particle['love_number'] = float(love_number)
        particle['fluid_love_number'] = float(fluid_love_number)
        particle['disk'] = disk.get()
        particle['migration_timescale'] = 0.0
        if type(disk) == DiskHost and (disk._data['Host']['inner_edge_distance'] != 0 or disk._data['Host']['outer_edge_distance'] != 0):
            if self._data['disk_host_particle_index'] == MAX_PARTICLES+1:
                self._data['disk_host_particle_index'] = self._data['n_particles']
            else:
                raise Exception("Only one body with a disk is allowed!")
        particle['position'] = position.get()
        particle['velocity'] = velocity.get()
        particle['spin'] = spin.get()
        particle['evolution_type'] = evolution_type.get()

        particle['inertial_position'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['inertial_velocity'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['inertial_acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}

        if type(evolution_type) != NonEvolving:
            self._data['evolving_particles_exist'] = True;

        particle['id'] = 0
        particle['acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['dangular_momentum_dt_due_to_tides'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['dangular_momentum_dt_induced_by_rotational_flattening'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['dangular_momentum_dt'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['lag_angle'] = 0.0
        particle['mass_g'] = particle['mass'] * K2
        particle['general_relativity_factor'] = 0.0
        particle['norm_velocity_vector'] = 0.0
        particle['tidal_acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['disk_interaction_acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['scalar_product_of_vector_position_with_stellar_spin'] = 0.0
        particle['scalar_product_of_vector_position_with_planetary_spin'] = 0.0
        particle['radial_component_of_the_force_induced_by_rotation'] = 0.0
        particle['orthogonal_component_of_the_force_induced_by_star_rotation'] = 0.0
        particle['orthogonal_component_of_the_force_induced_by_planet_rotation'] = 0.0
        particle['factor_for_the_force_induced_by_star_rotation'] = 0.0
        particle['factor_for_the_force_induced_by_planet_rotation'] = 0.0
        particle['radial_component_of_the_tidal_force'] = 0.0
        particle['orthogonal_component_of_the_tidal_force_due_to_stellar_tide'] = 0.0
        particle['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'] = 0.0
        particle['general_relativity_acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['radial_velocity'] = 0.0
        particle['wind_k_factor'] = wind_k_factor
        particle['wind_rotation_saturation'] = wind_rotation_saturation
        particle['wind_rotation_saturation_2'] = particle['wind_rotation_saturation']*particle['wind_rotation_saturation']
        if wind_k_factor != 0.:
            self._data['wind_effects_exist'] = True;

        particle['moment_of_inertia_ratio'] = 1.0
        particle['moment_of_inertia'] = particle['mass'] * particle['radius_of_gyration_2'] * particle['radius']*particle['radius']
        particle['dangular_momentum_dt_per_moment_of_inertia'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['wind_factor'] = 0.
        particle['distance'] = 0.0
        particle['acceleration_induced_by_rotational_flattering'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['norm_velocity_vector_2'] = 0.0
        particle['norm_spin_vector_2'] = particle['spin']['x']**2 + particle['spin']['y']**2 + particle['spin']['z']**2
        particle['denergy_dt'] = 0.
        particle['radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass'] = 0.

        evolver = evolution_type.get_evolver(self._data['initial_time'])
        if len(evolver['time']) > 0 and evolver['time'][0] > 0.:
            raise Exception("Your initial time ({} days | {:.2e} years) is smaller than the minimum allowed age of the star ({} days | {:.2e} years)".format(self._data['initial_time'], self._data['initial_time']/365.25, evolver['time'][0]+self._data['initial_time'], (evolver['time'][0]+self._data['initial_time'])/365.25));
        if len(evolver['time']) > 0 and evolver['time'][-1] < self._data['initial_time']:
            raise Exception("Your time limit ({} days | {:.2e} years) is greater than the maximum allowed age of the star ({} days | {:.2e} years)", self._data['initial_time'], self._data['initial_time']/365.25, evolver['time'][0], evolver['time'][0]/365.25)

        self._data['particles'].append(particle)
        self._data['particles_evolvers'].append(evolver)
        self._data['temporary_copied_particles_radiuses'].append(0.0)
        self._data['temporary_copied_particles_masses'].append(0.0)
        self._data['temporary_copied_particle_velocities'].append({u'x': 0.0, u'y': 0.0, u'z': 0.0})
        self._data['temporary_copied_particle_positions'].append({u'x': 0.0, u'y': 0.0, u'z': 0.0})
        self._data['n_particles'] += 1

    def add_brown_dwarf(self, mass, dissipation_factor_scale, position, velocity,  evolution_type, wind_k_factor=0., wind_rotation_saturation=0.):
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

        disk = NoDisk()

        radius_factor = 0.845649342247916
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 1.94e-1 # Brown dwarf

        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk, wind_k_factor=wind_k_factor, wind_rotation_saturation=wind_rotation_saturation)

    def add_solar_like(self, mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842):
        """
        Wind parametrisation (Bouvier 1997):

        wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
        """
        disk = NoDisk()
        self._add_solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type, disk, wind_k_factor, wind_rotation_saturation)

    def add_solar_like_with_disk(self, mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842):
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
        disk = DiskHost(disk_properties)
        self._add_solar_like_with_disk(mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type, disk, wind_k_factor, wind_rotation_saturation)

    def _add_solar_like_with_disk(self, mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type, disk, wind_k_factor, wind_rotation_saturation):
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
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk, wind_k_factor=wind_k_factor, wind_rotation_saturation=wind_rotation_saturation)


    def add_m_dwarf(self, mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type, wind_k_factor=0., wind_rotation_saturation=0.):
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

        disk = NoDisk()

        radius_factor = 0.845649342247916
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 2.0e-1 # M-dwarf
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk, wind_k_factor=wind_k_factor, wind_rotation_saturation=wind_rotation_saturation)


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

        disk = DiskInteraction(True)

        radius_of_gyration_2 = 2.54e-1 # Gas giant
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk)


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

        disk = DiskInteraction(True)

        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type, disk)


    def populate_inertial_frame(self):
        # Compute center of mass
        center_of_mass_position = Axes(0., 0., 0.)
        center_of_mass_velocity = Axes(0., 0., 0.)
        center_of_mass_acceleration = Axes(0., 0., 0.)
        center_of_mass_mass = 0.;

        for particle in self._data['particles']:
            center_of_mass_mass = get_center_of_mass_of_pair(center_of_mass_position,
                                                                    center_of_mass_velocity,
                                                                    center_of_mass_acceleration,
                                                                    center_of_mass_mass,
                                                                    particle)

        for particle in self._data['particles']:
            particle['inertial_position']['x'] = particle['position']['x'] - center_of_mass_position.x()
            particle['inertial_position']['y'] = particle['position']['y'] - center_of_mass_position.y()
            particle['inertial_position']['z'] = particle['position']['z'] - center_of_mass_position.z()
            particle['inertial_velocity']['x'] = particle['velocity']['x'] - center_of_mass_velocity.x()
            particle['inertial_velocity']['y'] = particle['velocity']['y'] - center_of_mass_velocity.y()
            particle['inertial_velocity']['z'] = particle['velocity']['z'] - center_of_mass_velocity.z()



    def assign_id(self):
        for i, particle in enumerate(self._data['particles']):
            particle['id'] = i

    def compute_general_relativity_factor(self):
        central_body_mass_g = self._data['particles'][0]['mass_g']
        for i, particle in enumerate(self._data['particles']):
            if self._data['consider_general_relativity'] == "None" or i == 0:
                # Disabled GR or Enabled GR and central body
                particle['general_relativity_factor'] = 0.0
            else:
                particle['general_relativity_factor'] =  central_body_mass_g*particle['mass_g'] / np.power(central_body_mass_g + particle['mass_g'], 2)

    def get(self):
        self.populate_inertial_frame()
        self.assign_id()
        self.compute_general_relativity_factor()

        # Add dummy particles to fill the vector
        n_dummy_particles = MAX_PARTICLES - self._data['n_particles']
        for i in xrange(n_dummy_particles):
            self.add_dummy_particle()

        data = self._data.copy()

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


