import json
from constants import *

class Integrator(object):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        self._data = {}
        self._data['time_step'] = float(time_step)
        self._data['half_time_step'] = self._data['time_step']*0.5
        self._data['universe'] = universe.get()
        self._data['last_spin'] = [{u'x': 0.0, u'y': 0.0, u'z': 0.0}]*MAX_PARTICLES # For spin integration with the midpoint method
        self._data['current_time'] = 0.0
        self._data['current_iteration'] = 0
        self._data['recovery_snapshot_period'] = float(recovery_snapshot_period)
        self._data['historic_snapshot_period'] = float(historic_snapshot_period)
        self._data['last_recovery_snapshot_time'] = -1.0
        self._data['last_historic_snapshot_time'] = -1.0
        self._data['n_historic_snapshots'] = 0
        self._data['hash'] = 0

    def write(self, filename):
        json.dump(self._data, open(filename, "w"))

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

class LeapFrog(Integrator):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(LeapFrog, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)

class Ias15(Integrator):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(Ias15, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)
        self._data['n_particles'] = universe._data['n_particles']
        self._data['integrator_iterations_max_exceeded'] = 0
        self._data['time_step_last_success'] = 0.
        self._data['b'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['br'] = [[0.,] * 3*MAX_PARTICLES,] * 7
        self._data['g'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['e'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['er'] = [[0.,] * 3*MAX_PARTICLES,] * 7
        self._data['at'] = [0.,] *  3*MAX_PARTICLES
        self._data['x0'] = [0.,] *  3*MAX_PARTICLES
        self._data['v0'] = [0.,] *  3*MAX_PARTICLES
        self._data['a0'] = [0.,] *  3*MAX_PARTICLES
        self._data['radius0'] = [0.,] *  MAX_PARTICLES
        self._data['radius_of_gyration_2_0'] = [0.,] *  MAX_PARTICLES
        self._data['sb'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['sbr'] = [[0.,] * 3*MAX_PARTICLES,] * 7
        self._data['sg'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['se'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['ser'] = [[0.,] * 3*MAX_PARTICLES,] * 7
        self._data['dangular_momentum_dtt'] = [0.,] *  3*MAX_PARTICLES
        self._data['spin0'] = [0.,] *  3*MAX_PARTICLES
        self._data['dangular_momentum_dt0'] = [0.,] *  3*MAX_PARTICLES
        self._data['angular_momentum0'] = [0.,] *  3*MAX_PARTICLES
        self._data['moment_of_inertia0'] = [0.,] * MAX_PARTICLES
        self._data['csx'] = [0.,] * 3*MAX_PARTICLES
        self._data['csv'] = [0.,] * 3*MAX_PARTICLES
        self._data['css'] = [0.,] * 3*MAX_PARTICLES
        self._data['s'] = [0.,] * 9

