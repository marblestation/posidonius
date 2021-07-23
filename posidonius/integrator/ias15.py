from posidonius.constants import *
from posidonius.integrator.common import Integrator

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
        self._data['sb'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['sbr'] = [[0.,] * 3*MAX_PARTICLES,] * 7
        self._data['sg'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['se'] = [[0.,] *  3*MAX_PARTICLES,] * 7
        self._data['ser'] = [[0.,] * 3*MAX_PARTICLES,] * 7
        self._data['dangular_momentum_dtt'] = [0.,] *  3*MAX_PARTICLES
        self._data['dangular_momentum_dt0'] = [0.,] *  3*MAX_PARTICLES
        self._data['angular_momentum0'] = [0.,] *  3*MAX_PARTICLES
        self._data['csx'] = [0.,] * 3*MAX_PARTICLES
        self._data['csv'] = [0.,] * 3*MAX_PARTICLES
        self._data['css'] = [0.,] * 3*MAX_PARTICLES
        self._data['s'] = [0.,] * 9

