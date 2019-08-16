from ..axes import Axes

class Tides(object):
    def __init__(self, variant, input_parameters=None):
        self._data = {
            "effect": "None",
            "parameters": {
                "input": {
                    "dissipation_factor_scale": 0.0,
                    "dissipation_factor": 0.0,
                    "love_number": 0.0,
                },
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
            self._data["effect"] = variant
            # Update default values, ignore non-recognised keys
            for key, value in input_parameters.iteritems():
                if key in self._data["parameters"]["input"]:
                    self._data["parameters"]["input"][key] = float(value)
                else:
                    print("Ignored parameter: {}".format(key))
            self._data["parameters"]["internal"]["scaled_dissipation_factor"] = self._data["parameters"]["input"]["dissipation_factor"] * self._data["parameters"]["input"]["dissipation_factor_scale"]
        elif variant in ("None", ):
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
        super(Disabled, self).__init__("None")

class OrbitingBody(Tides):
    def __init__(self, input_parameters):
        super(OrbitingBody, self).__init__("OrbitingBody", input_parameters=input_parameters)

class CentralBody(Tides):
    def __init__(self, input_parameters):
        super(CentralBody, self).__init__("CentralBody", input_parameters=input_parameters)


