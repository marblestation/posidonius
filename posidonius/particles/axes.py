import numpy as np

class Axes(object):
    def __init__(self, x, y, z):
        if type(x) is np.ndarray:
            if len(x) != len(y) or len(x) != len(z):
                raise Exception("Arrays length do not match!")
            self._data = {u'x': x, u'y': y, u'z': z}
        else:
            self._data = {u'x': float(x), u'y': float(y), u'z': float(z)}

    def get(self):
        return self._data.copy()

    def x(self):
        return self._data['x']

    def y(self):
        return self._data['y']

    def z(self):
        return self._data['z']

    def set_x(self, x):
        if type(x) is np.ndarray:
            if len(x) != len(self._data['y']) or len(x) != len(self._data['z']):
                raise Exception("Arrays length do not match!")
            self._data['x'] = x
        else:
            self._data['x'] = float(x)

    def set_y(self, y):
        if type(y) is np.ndarray:
            if len(y) != len(self._data['x']) or len(y) != len(self._data['z']):
                raise Exception("Arrays length do not match!")
            self._data['y'] = y
        else:
            self._data['y'] = float(y)

    def set_z(self, z):
        if type(z) is np.ndarray:
            if len(z) != len(self._data['x']) or len(z) != len(self._data['y']):
                raise Exception("Arrays length do not match!")
            self._data['z'] = z
        else:
            self._data['z'] = float(z)

