import six
from posidonius.particles.axes import Axes

class Wind(object):
    def __init__(self, variant, input_parameters=None):
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
                    "dangular_momentum_dt": Axes(0.0, 0.0, 0.0).get(),
                },
            },
        }
        if variant in ("Interaction", ):
            self._data["effect"] = variant
            # Update default values, ignore non-recognised keys
            for key, value in six.iteritems(input_parameters):
                if key in self._data["parameters"]["input"]:
                    self._data["parameters"]["input"][key] = float(value)
                else:
                    print("Ignored parameter: {}".format(key))
            self._data["parameters"]["internal"]["rotation_saturation_2"] = self._data["parameters"]["input"]["rotation_saturation"]*self._data["parameters"]["input"]["rotation_saturation"]
        elif variant in ("Disabled", ):
            self._data["effect"] = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Disabled(Wind):
    def __init__(self):
        super(Disabled, self).__init__("Disabled")

class Interaction(Wind):
    def __init__(self, input_parameters):
        super(Interaction, self).__init__("Interaction", input_parameters=input_parameters)

