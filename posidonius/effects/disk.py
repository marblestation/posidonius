from posidonius.particles.axes import Axes

class Disk(object):
    def __init__(self, variant, properties=None):
        self._data = {
            "effect": "Disabled",
            "parameters": {
                "internal": {
                    "distance": 0.0,
                    "migration_timescale": 0.0,
                    "norm_velocity_vector": 0.0,
                    "norm_velocity_vector_2": 0.0,
                },
                "output": {
                    "acceleration": Axes(0.0, 0.0, 0.0).get(),
                },
            },
            "coordinates": {
                "position": Axes(0.0, 0.0, 0.0).get(),
                "velocity": Axes(0.0, 0.0, 0.0).get(),
            },
        }
        if variant in ("CentralBody",):
            if properties is None:
                raise Exception("Disk central body requires properties")
            self._data["effect"] = {
                variant: properties,
            }
        elif variant in ("OrbitingBody", "Disabled", ):
            self._data["effect"] = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Disabled(Disk):
    def __init__(self):
        super(Disabled, self).__init__("Disabled", properties=None)

class OrbitingBody(Disk):
    def __init__(self):
        super(OrbitingBody, self).__init__("OrbitingBody")

class CentralBody(Disk):
    def __init__(self, input_properties):
        properties = {
            "alpha": 0.0,
            "inner_edge_distance": 0.0,
            "lifetime": 0.0,
            "mean_molecular_weight": 0.0,
            "outer_edge_distance": 0.0,
            "surface_density_normalization": 0.0,
        }
        # Update default values, ignore non-recognised keys
        for key, value in input_properties.iteritems():
            if key in properties:
                properties[key] = value
        super(CentralBody, self).__init__("CentralBody", properties)

