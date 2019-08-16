from axes import Axes
import posidonius.effects as effects
from posidonius.constants import K2

class Particle(object):
    def __init__(self, mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type):
        self._data = {
            "dangular_momentum_dt": Axes(0.0, 0.0, 0.0).get(),
            "dangular_momentum_dt_per_moment_of_inertia": Axes(0.0, 0.0, 0.0).get(),
            "disk": disk.get(),
            "evolution_type": evolution_type.get(),
            "general_relativity": general_relativity.get(),
            "heliocentric_distance": 0.0,
            "heliocentric_norm_velocity_vector": 0.0,
            "heliocentric_norm_velocity_vector_2": 0.0,
            "heliocentric_position": position.get(),
            "heliocentric_radial_velocity": 0.0,
            "heliocentric_velocity": velocity.get(),
            "id": 0,
            "inertial_acceleration": Axes(0.0, 0.0, 0.0).get(),
            "inertial_position": Axes(0.0, 0.0, 0.0).get(),
            "inertial_velocity": Axes(0.0, 0.0, 0.0).get(),
            "mass": float(mass),
            "mass_g": float(mass)*K2,
            "moment_of_inertia": float(mass)*float(radius_of_gyration_2)*float(radius)*float(radius),
            "moment_of_inertia_ratio": 1.0,
            "norm_spin_vector_2": spin.x()**2 + spin.y()**2 + spin.z()**2,
            "radius": float(radius),
            "radius_of_gyration_2": float(radius_of_gyration_2),
            "rotational_flattening": rotational_flattening.get(),
            "spin": spin.get(),
            "tides": tides.get(),
            "wind": wind.get(),
        }
        self._effects = {
            'tides': tides,
            'rotational_flattening': rotational_flattening,
            'general_relativity': general_relativity,
            'wind': wind,
            'disk': disk,
        }
        self._evolution_type = evolution_type

    def effects(self):
        return [type(value) for key, value in self._effects.iteritems()]

    def general_relativity_implementation(self):
        return self._effects['general_relativity'].implementation()

    def evolution_type(self):
        return type(self._evolution_type)

    def get_evolver(self, initial_time):
        return self._evolution_type.get_evolver(initial_time)

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class DummyParticle(Particle):
    def __init__(self):
        mass = 0.0
        radius = 0.0
        radius_of_gyration_2 = 0.0
        position = Axes(0.0, 0.0, 0.0)
        velocity = Axes(0.0, 0.0, 0.0)
        spin = Axes(0.0, 0.0, 0.0)
        tides = effects.tides.Disabled()
        rotational_flattening = effects.rotational_flattening.Disabled()
        general_relativity = effects.general_relativity.Disabled()
        wind = effects.wind.Disabled()
        disk = effects.disk.Disabled()
        evolution_type = effects.evolution.NonEvolving()
        super(DummyParticle, self).__init__(mass, radius, radius_of_gyration_2, position, velocity, spin, tides, rotational_flattening, general_relativity, wind, disk, evolution_type)
