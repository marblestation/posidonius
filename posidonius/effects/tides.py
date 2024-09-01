import six
from posidonius.particles.axes import Axes

class Tides(object):
    def __init__(self, variant, tidal_model=None):
        self._model = tidal_model
        self._data = {
            "effect": "Disabled",
            "parameters": {
                "internal": {
                    "denergy_dt": 0.0,
                    "distance": 0.0,
                    "lag_angle": 0.0,
                    "orthogonal_component_of_the_tidal_force_due_to_planetary_tide": 0.0,
                    "orthogonal_component_of_the_tidal_force_due_to_stellar_tide": 0.0,
                    "radial_component_of_the_tidal_force": 0.0,
                    "radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass": 0.0,
                    "radial_velocity": 0.0,
                    "scalar_product_of_vector_position_with_planetary_spin": 0.0,
                    "scalar_product_of_vector_position_with_stellar_spin": 0.0,
                    "scaled_dissipation_factor": 0.0,
                    "shape": Axes(0.0, 0.0, 0.0).get(),
                },
                "output": {
                    "acceleration": Axes(0.0, 0.0, 0.0).get(),
                    "dangular_momentum_dt": Axes(0.0, 0.0, 0.0).get(),
                },
            },
            "coordinates": {
                "position": Axes(0.0, 0.0, 0.0).get(),
                "velocity": Axes(0.0, 0.0, 0.0).get(),
            },
        }
        if variant in ("CentralBody", "OrbitingBody", ):
            self._data["effect"] = {variant: tidal_model.get()}
            if isinstance(tidal_model, ConstantTimeLag):
                self._data["parameters"]["internal"]["scaled_dissipation_factor"] = self._data["effect"][variant]["ConstantTimeLag"]["dissipation_factor"] * self._data["effect"][variant]["ConstantTimeLag"]["dissipation_factor_scale"]
        elif variant in ("Disabled", ):
            self._data["effect"] = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Disabled(Tides):
    def __init__(self):
        super(Disabled, self).__init__("Disabled")

class OrbitingBody(Tides):
    def __init__(self, tidal_model):
        super(OrbitingBody, self).__init__("OrbitingBody", tidal_model=tidal_model)

class CentralBody(Tides):
    def __init__(self, tidal_model):
        super(CentralBody, self).__init__("CentralBody", tidal_model=tidal_model)

class ConstantTimeLag(object):
    def __init__(self, input_parameters):
        self._data = {
            "ConstantTimeLag": {
                "dissipation_factor_scale": 0.0,
                "dissipation_factor": 0.0,
                "love_number": 0.0,
            },
        }
        # Update default values, ignore non-recognised keys
        for key, value in six.iteritems(input_parameters):
            if key in self._data["ConstantTimeLag"]:
                self._data["ConstantTimeLag"][key] = float(value)
            else:
                print("Ignored parameter: {}".format(key))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class CreepCoplanar(object):
    def __init__(self, input_parameters):
        self._data = {
            "CreepCoplanar": {
                "uniform_viscosity_coefficient": 0.0,
            },
        }
        # Update default values, ignore non-recognised keys
        for key, value in six.iteritems(input_parameters):
            if key in self._data["CreepCoplanar"]:
                self._data["CreepCoplanar"][key] = float(value)
            else:
                print("Ignored parameter: {}".format(key))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Kaula(object):
    def __init__(self, input_parameters):
        self._data = {
            "Kaula": {
                "love_number_excitation_frequency": [0.] * 1024,
                "imaginary_part_love_number": [0.] * 1024,
                "real_part_love_number": [0.] * 1024,
                "num_datapoints": 0.0,
                "kaula_tidal_force": Axes( 0.0, 0.0, 0.0).get(),
            },
        }
        # Update default values, ignore non-recognised keys
        for key, value in six.iteritems(input_parameters):
            if key in self._data["Kaula"]:
                if isinstance(value, (tuple, list)):
                    self._data["Kaula"][key] = [float(v) for v in value]
                else:
                    self._data["Kaula"][key] = float(value)
            else:
                print("Ignored parameter: {}".format(key))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()
