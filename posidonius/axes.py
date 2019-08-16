
class Axes(object):
    def __init__(self, x, y, z):
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
        self._data['x'] = float(x)

    def set_y(self, y):
        self._data['y'] = float(y)

    def set_z(self, z):
        self._data['z'] = float(z)

