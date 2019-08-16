from posidonius.particles.axes import Axes

class GeneralRelativity(object):
    def __init__(self, variant, implementation=None):
        self._data = {
            "effect": "Disabled",
            "parameters": {
                "internal": {
                    "distance": 0.0,
                    "factor": 0.0,
                    "norm_velocity_vector": 0.0,
                    "norm_velocity_vector_2": 0.0,
                    "radial_velocity": 0.0,
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
        if variant in ("CentralBody",):
            if implementation is None:
                raise Exception("General relativity central body requires an implementation")
            self._data["effect"] = {
                variant: implementation,
            }
        elif variant in ("OrbitingBody", "Disabled", ):
            self._data["effect"] = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def implementation(self):
        return self._data['effect'].get('CentralBody')

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Disabled(GeneralRelativity):
    def __init__(self):
        super(Disabled, self).__init__("Disabled")

class OrbitingBody(GeneralRelativity):
    def __init__(self):
        super(OrbitingBody, self).__init__("OrbitingBody")

class CentralBody(GeneralRelativity):
    def __init__(self, implementation):
        if implementation not in ("Kidder1995", "Anderson1975", "Newhall1983"):
            raise Exception("Unknown general relativity implementation: {}".format(implementation))
        super(CentralBody, self).__init__("CentralBody", implementation=implementation)

