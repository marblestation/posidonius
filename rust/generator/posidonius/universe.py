from axes import Axes
from integrator import WHFast
from constants import *
from evolution_type import NonEvolving, BrownDwarf, MDwarf, Jupiter, SolarLikeEvolvingDissipation, SolarLikeConstantDissipation
from tools import calculate_spin

class Universe(object):
    def __init__(self, initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_general_relativy):
        self._time_step = time_step
        self._recovery_snapshot_period = recovery_snapshot_period
        self._historic_snapshot_period = historic_snapshot_period
        self._data = {}
        self._data['time_limit'] = float(time_limit)
        self._data['initial_time'] = float(initial_time)
        self._data['consider_tides'] = consider_tides
        self._data['consider_rotational_flattening'] = consider_rotational_flattening
        self._data['consider_general_relativy'] = consider_general_relativy
        self._data['particles'] = []
        self._data['particles_evolvers'] = []
        self._data['n_particles'] = 0
        self._data['evolving_particles_exist'] = False
        self._data['star_planet_dependent_dissipation_factors'] = {}
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
        position = Axes(0., 0., 0.)
        velocity = Axes(0., 0., 0.)
        spin = Axes(0., 0., 0.)
        evolution_type = NonEvolving()
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)
        self._data['n_particles'] -= 1

    def add_particle(self, mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type):
        if self._data['n_particles'] == MAX_PARTICLES:
            raise Exception("Maximum number of particles reached: {}".format(MAX_PARTICLES))
        particle = {}
        particle['mass'] = float(mass)
        particle['radius'] = float(radius)
        particle['scaled_dissipation_factor'] = float(dissipation_factor)*float(dissipation_factor_scale)
        particle['dissipation_factor_scale'] = float(dissipation_factor_scale)
        particle['radius_of_gyration_2'] = float(radius_of_gyration_2)
        particle['love_number'] = float(love_number)
        particle['fluid_love_number'] = float(fluid_love_number)
        particle['position'] = position.get()
        particle['velocity'] = velocity.get()
        particle['spin'] = spin.get()
        particle['evolution_type'] = evolution_type.get()

        if type(evolution_type) != NonEvolving:
            self._data['evolving_particles_exist'] = True;

        particle['id'] = self._data['n_particles']
        particle['acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['torque'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['lag_angle'] = 0.0
        particle['general_relativity_factor'] = 0.0
        particle['norm_velocity_vector'] = 0.0
        particle['tidal_acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
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
        particle['moment_of_inertia_ratio'] = 1.0
        particle['mass_g'] = particle['mass'] * K2
        particle['dspin_dt'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['wind_factor'] = 0.
        particle['distance'] = 0.0
        particle['acceleration_induced_by_rotational_flattering'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle['norm_velocity_vector_2'] = 0.0
        particle['norm_spin_vector_2'] = 0.0
        particle['denergy_dt'] = 0.

        evolver = evolution_type.get_evolver(self._data['initial_time'])
        if len(evolver['time']) > 0 and evolver['time'][0] > 0.:
            raise Exception("Your initial time ({} days) is smaller than the minimum allowed age of the star ({} days)".format(self._data['initial_time'], evolver['time'][0]+self._data['initial_time']));
        if len(evolver['time']) > 0 and evolver['time'][-1] < self._data['initial_time']:
            raise Exception("Your time limit ({} days) is greater than the maximum allowed age of the star ({} days)", self._data['initial_time'], evolver['time'][0])

        self._data['particles'].append(particle)
        self._data['particles_evolvers'].append(evolver)
        self._data['temporary_copied_particles_radiuses'].append(0.0)
        self._data['temporary_copied_particles_masses'].append(0.0)
        self._data['temporary_copied_particle_velocities'].append({u'x': 0.0, u'y': 0.0, u'z': 0.0})
        self._data['temporary_copied_particle_positions'].append({u'x': 0.0, u'y': 0.0, u'z': 0.0})
        self._data['n_particles'] += 1

    def add_brown_dwarf(self, mass, dissipation_factor_scale, position, velocity,  evolution_type):
        rotation_period = None
        love_number = None
        if type(evolution_type) == NonEvolving:
            rotation_period = 70.0 # hours
            love_number = 0.307 # BrownDwarf
        elif type(evolution_type) == BrownDwarf:
            mass = evolution_type._data['fields'][0]

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
                raise Exception("The evolution type BrownDwarf does not support a mass of {} Msun!".format(mass))
        else:
            raise Exception("Evolution type should be BrownDwarf or NonEvolving!")

        angular_frequency = TWO_PI/(rotation_period/24.) # days^-1
        inclination = 0.
        obliquity = 0.
        spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity)

        fluid_love_number = love_number
        # BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        dissipation_factor = 2.006*3.845764e4 # -60+64

        radius_factor = 0.845649342247916
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 1.94e-1 # Brown dwarf

        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)


    def add_solar_like(self, mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type):
        if type(evolution_type) not in (SolarLikeEvolvingDissipation, SolarLikeConstantDissipation):
            raise Exception("Evolution type should be SolarLikeEvolvingDissipation or SolarLikeConstantDissipation!")

        # Typical rotation period: 24 hours
        angular_frequency = TWO_PI/(rotation_period/24.) # days^-1
        inclination = 0.
        obliquity = 0.
        spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity)

        fluid_love_number = love_number
        # Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        dissipation_factor = 4.992*3.845764e-2 # -66+64

        radius_factor = 1.
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 5.9e-2 # Sun
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)


    def add_m_dwarf(self, mass, dissipation_factor_scale, position, velocity, rotation_period, evolution_type):
        if type(evolution_type) not in (MDwarf, NonEvolving):
            raise Exception("Evolution type should be MDwarf or NonEvolving!")

        # Typical rotation period: 70 hours
        angular_frequency = TWO_PI/(rotation_period/24.) # days^-1
        inclination = 0.
        obliquity = 0.
        spin = calculate_spin(angular_frequency, inclination, obliquity, position, velocity)

        love_number = 0.307 # M Dwarf
        fluid_love_number = love_number

        # BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        dissipation_factor = 2.006*3.845764e4 # -60+64

        radius_factor = 0.845649342247916
        radius = radius_factor * R_SUN
        radius_of_gyration_2 = 2.0e-1 # M-dwarf
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)


    def add_jupiter_like(self, mass, dissipation_factor_scale, position, velocity, spin, evolution_type):
        if type(evolution_type) not in (Jupiter, NonEvolving):
            raise Exception("Evolution type should be Jupiter or NonEvolving!")

        # Typical rotation period: 9.8 hours
        love_number = 0.380 # Gas giant
        fluid_love_number = love_number

        radius_factor = 10.9 # Jupiter in R_EARTH
        radius = radius_factor * R_EARTH

        # TODO: What k2pdelta/dissipation_factor is the recommended?
        #k2pdelta = 8.101852e-9 # Gas giant
        k2pdelta = 2.893519e-7 # Gas giant for Jupiter: 2-3d-2 s, here in day (Leconte)
        dissipation_factor = 2. * K2 * k2pdelta/(3. * np.power(radius))
        #dissipation_factor = 2.006*3.845764e4 // Gas giant

        radius_of_gyration_2 = 2.54e-1 # Gas giant
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)


    def add_earth_like(self, mass, dissipation_factor_scale, position, velocity, spin, evolution_type):
        if type(evolution_type) not in (NonEvolving,):
            raise Exception("Evolution type should be NonEvolving!")

        # Typical rotation period: 24 hours
        love_number = 0.299 # Earth
        fluid_love_number = 0.9532 # Earth

        # Earth-like => mass-radius relationship from Fortney 2007
        radius_factor = (0.0592*0.7+0.0975) * np.power(np.log10(mass) + np.log10(M2EARTH), 2) \
                                 + (0.2337*0.7+0.4938) * (np.log10(mass) + np.log10(M2EARTH)) \
                                 + 0.3102*0.7+0.7932
        radius = radius_factor * R_EARTH
        radius_of_gyration_2 = 3.308e-1 # Earth type planet
        k2pdelta = 2.465278e-3 # Terrestrial planets
        dissipation_factor = 2. * K2 * k2pdelta/(3. * np.power(radius, 5))

        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)

    def add_terrestrial(self, mass, radius_factor, dissipation_factor_scale, position, velocity, spin, evolution_type):
        if type(evolution_type) not in (NonEvolving,):
            raise Exception("Evolution type should be NonEvolving!")

        # Typical rotation period: 24 hours
        love_number = 0.299 # Earth
        fluid_love_number = 0.9532 # Earth

        radius = radius_factor * R_EARTH
        radius_of_gyration_2 = 3.308e-1 # Earth type planet
        k2pdelta = 2.465278e-3 # Terrestrial planets
        dissipation_factor = 2. * K2 * k2pdelta/(3. * np.power(radius, 5))

        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)


    def add_gas_giant(self, mass, radius_factor, dissipation_factor_scale, position, velocity, spin, evolution_type):
        if type(evolution_type) not in (NonEvolving):
            raise Exception("Evolution type should be NonEvolving!")

        # Typical rotation period: 9.8 hours
        love_number = 0.380 # Gas giant
        fluid_love_number = love_number

        radius = radius_factor * R_EARTH

        # TODO: What k2pdelta/dissipation_factor is the recommended?
        #k2pdelta = 8.101852e-9 # Gas giant
        k2pdelta = 2.893519e-7 # Gas giant for Jupiter: 2-3d-2 s, here in day (Leconte)
        dissipation_factor = 2. * K2 * k2pdelta/(3. * np.power(radius))
        #dissipation_factor = 2.006*3.845764e4 // Gas giant

        radius_of_gyration_2 = 2.54e-1 # Gas giant
        self.add_particle(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number, position, velocity, spin, evolution_type)






    def get(self):
        # Add dummy particles to fill the vector
        n_dummy_particles = MAX_PARTICLES - self._data['n_particles']
        for i in xrange(n_dummy_particles):
            self.add_dummy_particle()

        data = self._data.copy()

        # Forget the dummy particles
        if n_dummy_particles > 0:
            self._data['particles'] = self._data['particles'][:-n_dummy_particles]
            self._data['particles_evolvers'] = self._data['particles_evolvers'][:-n_dummy_particles]
            self._data['n_particles'] -= n_dummy_particles
        return data

    def write(self, filename, integrator="WHFast"):
        if integrator == "WHFast":
            universe_integrator = WHFast(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        elif integrator == "Ias15":
            universe_integrator = Ias15(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        elif integrator == "LeapFrog":
            universe_integrator = LeapFrog(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        else:
            raise Exception("Unknown integtrator '{}'".format(integrator))


