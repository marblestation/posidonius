
class DiskType(object):
    def __init__(self, variant, properties=None, enabled=False):
        self._data = {}
        if variant in ("Host"):
            self._data[variant] = properties
        elif variant in ("Interaction", ):
            self._data[variant] = enabled
        elif variant in ("None",):
            self._data = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class NoDisk(DiskType):
    def __init__(self):
        super(NoDisk, self).__init__("None", properties=None)

class DiskInteraction(DiskType):
    def __init__(self, enabled):
        super(DiskInteraction, self).__init__("Interaction", enabled=enabled)

class DiskHost(DiskType):
    def __init__(self, input_properties):
        properties = {
            'inner_edge_distance': 0.0,
            'outer_edge_distance': 0.0,
            'lifetime': 0.0,
            'alpha': 0.0,
            'surface_density_normalization': 0.0,
            'mean_molecular_weight':  0.0,
        }
        # Update default values, ignore non-recognised keys
        for key, value in input_properties.iteritems():
            if key in properties:
                properties[key] = value
        super(DiskHost, self).__init__("Host", properties)

