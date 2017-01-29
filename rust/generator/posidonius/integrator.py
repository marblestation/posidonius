import json

class Integrator(object):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        self._data = {}
        self._data['time_step'] = float(time_step)
        self._data['half_time_step'] = self._data['time_step']/2.
        self._data['universe'] = universe.get()
        self._data['current_time'] = 0.0
        self._data['current_iteration'] = 0
        self._data['recovery_snapshot_period'] = 3652500.0
        self._data['historic_snapshot_period'] = 36525.0
        self._data['last_recovery_snapshot_time'] = -1.0
        self._data['last_historic_snapshot_time'] = -1.0
        self._data['n_historic_snapshots'] = 0
        self._data['hash'] = 0

    def write(self, filename):
        json.dump(self._data, open(filename, "w"))


class WHFast(Integrator):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(WHFast, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)
        self._data['timestep_warning'] = 0
        self._data['safe_mode'] = True
        self._data['universe_heliocentric'] = self._data['universe'].copy()
        self._data['recalculate_heliocentric_this_timestep'] = True
        self._data['recalculate_heliocentric_but_not_synchronized_warning'] = 0
        self._data['set_to_center_of_mass'] = False
        self._data['corrector'] = 0
        self._data['is_synchronized'] = True

class LeapFrog(Integrator):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(LeapFrog, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)

class Ias15(Integrator):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(Ias15, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)
        self._data['n_particles'] = universe._data['n_particles']
        self._data['integrator_iterations_max_exceeded'] = 0
        self._data['time_step_last_success'] = 0.
        self._data['b'] = [[0.,] * self._data['n_particles'],] * 7
        self._data['br'] = [[0.,] * self._data['n_particles'],] * 7
        self._data['g'] = [[0.,] * self._data['n_particles'],] * 7
        self._data['e'] = [[0.,] * self._data['n_particles'],] * 7
        self._data['er'] = [[0.,] * self._data['n_particles'],] * 7
        self._data['at'] = [0.,] * self._data['n_particles']
        self._data['x0'] = [0.,] * self._data['n_particles']
        self._data['v0'] = [0.,] * self._data['n_particles']
        self._data['a0'] = [0.,] * self._data['n_particles']
        self._data['csx'] = [0.,] * self._data['n_particles']
        self._data['csv'] = [0.,] * self._data['n_particles']
        self._data['s'] = [0.,] * 9

