from posidonius.constants import *
from posidonius.integrator.common import Integrator

class LeapFrog(Integrator):
    def __init__(self, time_step, recovery_snapshot_period, historic_snapshot_period, universe):
        super(LeapFrog, self).__init__(time_step, recovery_snapshot_period, historic_snapshot_period, universe)
        self._data['half_time_step'] = self._data['time_step']*0.5

