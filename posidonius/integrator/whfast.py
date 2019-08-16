from posidonius.constants import *
from posidonius.integrator.common import Integrator

class CoordinatesType(object):
    def __init__(self, variant):
        self._data = {}
        if variant in ("Jacobi", "DemocraticHeliocentric", "WHDS"):
            self._data = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

    def get(self):
        return self._data

class WHFast(Integrator):

    def __init__(self, alternative_coordinates, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(WHFast, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)
        self._data['half_time_step'] = self._data['time_step']*0.5
        self._data['last_spin'] = [{u'x': 0.0, u'y': 0.0, u'z': 0.0}]*MAX_PARTICLES # For spin integration with the midpoint method
        self._data['particle_spin_errors'] = [{u'x': 0.0, u'y': 0.0, u'z': 0.0}]*MAX_PARTICLES # For spin integration with the midpoint method
        self._data['timestep_warning'] = 0
        self._data['alternative_coordinates_type'] = CoordinatesType(alternative_coordinates).get()
        self._data['particles_alternative_coordinates'] = []
        particle_alternative_coordinates = {}
        particle_alternative_coordinates['mass'] = 0.
        particle_alternative_coordinates['mass_g'] = 0.
        particle_alternative_coordinates['position'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle_alternative_coordinates['velocity'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        particle_alternative_coordinates['acceleration'] = {u'x': 0.0, u'y': 0.0, u'z': 0.0}
        self._data['particles_alternative_coordinates'] = [particle_alternative_coordinates] * MAX_PARTICLES

