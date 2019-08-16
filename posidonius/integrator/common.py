import json

class Integrator(object):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        self._data = {}
        self._data['time_step'] = float(time_step)
        self._data['universe'] = universe.get()
        self._data['current_time'] = 0.0
        self._data['current_iteration'] = 0
        self._data['recovery_snapshot_period'] = float(recovery_snapshot_period)
        self._data['historic_snapshot_period'] = float(historic_snapshot_period)
        self._data['last_recovery_snapshot_time'] = -1.0
        self._data['last_historic_snapshot_time'] = -1.0
        self._data['n_historic_snapshots'] = 0
        self._data['hash'] = 0

    def write(self, filename):
        json.dump(self._data, open(filename, "w"), indent=2, sort_keys=True)

