import six
import numpy as np
from posidonius.particles.axes import Axes
import posidonius.effects as effects
from posidonius.constants import K2
from posidonius.effects.evolution import NonEvolving
from posidonius import effects

class Reference(object):
    def __init__(self, variant, index=None):
        self._data = {
            "effect": "Disabled",
            "parameters": {
                "input": {
                    "k_factor": 0.0,
                    "rotation_saturation": 0.0,
                },
                "internal": {
                    "rotation_saturation_2": 0.0,
                },
                "output": {
                    "factor": 0.0,
                },
            },
        }
        if variant in ("Particle", ):
            self._data[variant] = int(index)
        if variant in ("MostMassiveParticle", ):
            self._data = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class MostMassiveParticle(Reference):
    def __init__(self):
        super(MostMassiveParticle, self).__init__("MostMassiveParticle")

class ReferenceParticle(Reference):
    def __init__(self, index):
        super(ReferenceParticle, self).__init__("ReferenceParticle", index=index)


class Particle(object):
    def __init__(self, mass, radius, radius_of_gyration, position, velocity, spin):
        radius_of_gyration_2 = float(radius_of_gyration)**2
        reference = MostMassiveParticle()
        tides = effects.tides.Disabled()
        rotational_flattening = effects.rotational_flattening.Disabled()
        general_relativity = effects.general_relativity.Disabled()
        wind = effects.wind.Disabled()
        disk = effects.disk.Disabled()
        evolution = NonEvolving()
        moment_of_inertia = float(mass)*float(radius_of_gyration_2)*float(radius)*float(radius)
        self._data = {
            "angular_momentum": Axes(moment_of_inertia*spin.x(), moment_of_inertia*spin.y(), moment_of_inertia*spin.z()).get(),
            "dangular_momentum_dt": Axes(0.0, 0.0, 0.0).get(),
            "disk": disk.get(),
            "evolution": evolution.get(),
            "general_relativity": general_relativity.get(),
            "heliocentric_distance": 0.0,
            "heliocentric_norm_velocity_vector": 0.0,
            "heliocentric_norm_velocity_vector_2": 0.0,
            "heliocentric_position": position.get(),
            "heliocentric_radial_velocity": 0.0,
            "heliocentric_velocity": velocity.get(),
            "id": 0,
            "inertial_acceleration": Axes(0.0, 0.0, 0.0).get(),
            "inertial_additional_acceleration": Axes(0.0, 0.0, 0.0).get(),
            "inertial_position": Axes(0.0, 0.0, 0.0).get(),
            "inertial_velocity": Axes(0.0, 0.0, 0.0).get(),
            "mass": float(mass),
            "mass_g": float(mass)*K2,
            "moment_of_inertia": moment_of_inertia,
            "norm_spin_vector_2": spin.x()**2 + spin.y()**2 + spin.z()**2,
            "radius": float(radius),
            "radius_of_gyration_2": float(radius_of_gyration_2),
            "reference": reference.get(),
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
        self._evolution = evolution

    def effects(self):
        return [type(value) for key, value in six.iteritems(self._effects)]

    def general_relativity_implementation(self):
        return self._effects['general_relativity'].implementation()

    def evolution(self):
        return type(self._evolution)

    def get_evolver(self, initial_time):
        return self._evolution.get_evolver(initial_time)

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

    def set_tides(self, tides):
        self._data["tides"] = tides.get()
        self._effects["tides"] = tides
        self.check_uniform_viscosity_coefficient()

    def set_rotational_flattening(self, rotational_flattening):
        self._data["rotational_flattening"] = rotational_flattening.get()
        self._effects["rotational_flattening"] = rotational_flattening
        self.check_uniform_viscosity_coefficient()

    def check_uniform_viscosity_coefficient(self):
        # If creep coplanar tides and rotational flattening are set, both need to have the same
        # uniform viscosity coefficient parameter
        if self._effects["tides"] and self._effects["rotational_flattening"]:
            disabled_tides = isinstance(self._effects["tides"], effects.tides.Disabled)
            disabled_rotational_flattening = isinstance(self._effects["rotational_flattening"], effects.rotational_flattening.Disabled)
            if not disabled_tides and not disabled_rotational_flattening:
                creep_coplanar_tides = isinstance(self._effects["tides"]._model, effects.tides.CreepCoplanar)
                creep_coplanar_rotational_flattening = isinstance(self._effects["rotational_flattening"]._model, effects.rotational_flattening.CreepCoplanar)
                if (creep_coplanar_tides and not creep_coplanar_rotational_flattening) or (not creep_coplanar_tides and creep_coplanar_rotational_flattening):
                    raise Exception("When using Creep Coplanar Tidal or rotational flattening effects, both effects need to be Creep Coplanar and not just one of them (e.g., it cannot be mixed with ConstantTimeLag or OblateSpheroid).")
                elif creep_coplanar_tides and creep_coplanar_rotational_flattening:
                    uniform_viscosity_coefficient_diff = self._effects["rotational_flattening"]._model._data['CreepCoplanar']['uniform_viscosity_coefficient'] - self._effects["tides"]._model._data['CreepCoplanar']['uniform_viscosity_coefficient']
                    if abs(uniform_viscosity_coefficient_diff) > 1.e-16:
                        raise Exception("When using Creep Coplanar Tidal and rotational flattening effects, the uniform viscosity coefficient must be identical.")

    def set_general_relativity(self, general_relativity):
        self._data["general_relativity"] = general_relativity.get()
        self._effects["general_relativity"] = general_relativity

    def set_wind(self, wind):
        self._data["wind"] = wind.get()
        self._effects["wind"] = wind

    def set_disk(self, disk):
        self._data["disk"] = disk.get()
        self._effects["disk"] = disk

    def set_evolution(self, evolution):
        self._data["evolution"] = evolution.get()
        self._evolution = evolution

    def set_reference(self, reference):
        self._data["reference"] = reference.get()


class DummyParticle(Particle):
    def __init__(self):
        mass = 0.0
        radius = 0.0
        radius_of_gyration_2 = 0.0
        position = Axes(0.0, 0.0, 0.0)
        velocity = Axes(0.0, 0.0, 0.0)
        spin = Axes(0.0, 0.0, 0.0)
        reference = MostMassiveParticle()
        tides = effects.tides.Disabled()
        rotational_flattening = effects.rotational_flattening.Disabled()
        general_relativity = effects.general_relativity.Disabled()
        wind = effects.wind.Disabled()
        disk = effects.disk.Disabled()
        evolution = effects.evolution.NonEvolving()
        super(DummyParticle, self).__init__(mass, radius, radius_of_gyration_2, position, velocity, spin)
        super(DummyParticle, self).set_tides(tides)
        super(DummyParticle, self).set_rotational_flattening(rotational_flattening)
        super(DummyParticle, self).set_general_relativity(general_relativity)
        super(DummyParticle, self).set_wind(wind)
        super(DummyParticle, self).set_disk(disk)
        super(DummyParticle, self).set_evolution(evolution)
        super(DummyParticle, self).set_reference(reference)
