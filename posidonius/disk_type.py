
class DiskType(object):
    def __init__(self, variant, properties=None):
        self._data = {}
        if variant in ("Properties"):
            self._data[variant] = properties
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

class Disk(DiskType):
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
        super(Disk, self).__init__("Properties", properties)

