import six
from posidonius.particles.axes import Axes

class RotationalFlattening(object):
    def __init__(self, variant, rotational_flattening_model=None):
        self._model = rotational_flattening_model
        self._data = {
            "effect": "Disabled",
            "parameters": {
                "internal": {
                    "distance": 0.0,
                    "factor_for_the_force_induced_by_planet_rotation": 0.0,
                    "factor_for_the_force_induced_by_star_rotation": 0.0,
                    "orthogonal_component_of_the_force_induced_by_planet_rotation": 0.0,
                    "orthogonal_component_of_the_force_induced_by_star_rotation": 0.0,
                    "radial_component_of_the_force_induced_by_rotation": 0.0,
                    "scalar_product_of_vector_position_with_planetary_spin": 0.0,
                    "scalar_product_of_vector_position_with_stellar_spin": 0.0,
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
            self._data["effect"] = {variant: rotational_flattening_model.get()}
        elif variant in ("Disabled", ):
            self._data["effect"] = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Disabled(RotationalFlattening):
    def __init__(self):
        super(Disabled, self).__init__("Disabled")

class OrbitingBody(RotationalFlattening):
    def __init__(self, rotational_flattening_model):
        super(OrbitingBody, self).__init__("OrbitingBody", rotational_flattening_model=rotational_flattening_model)

class CentralBody(RotationalFlattening):
    def __init__(self, rotational_flattening_model):
        super(CentralBody, self).__init__("CentralBody", rotational_flattening_model=rotational_flattening_model)

class OblateSpheroid(object):
    def __init__(self, input_parameters):
        self._data = {
            "OblateSpheroid": {
                "love_number": 0.0,
            },
        }
        # Update default values, ignore non-recognised keys
        for key, value in six.iteritems(input_parameters):
            if key in self._data["OblateSpheroid"]:
                self._data["OblateSpheroid"][key] = float(value)
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

