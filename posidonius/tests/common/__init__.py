import os
import errno
from posidonius.tests.common.stars import *
from posidonius.tests.common.planets import *


def _mkdir_p(path):
    """
    Creates a directory. Same behaviour as 'mkdir -p'.
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

def setup(current_dirname, current_filename, current_function_name):
    test_name = "{}-{}".format(current_filename, current_function_name)
    tmp_data_dirname = os.path.join(current_dirname, "data/tmp/{0}/".format(test_name))
    expected_data_dirname = os.path.join(current_dirname, "data/{0}/".format(test_name))
    _mkdir_p(expected_data_dirname)
    _mkdir_p(tmp_data_dirname)
    json_filename = os.path.join(tmp_data_dirname, "case.json")
    expected_json_filename = os.path.join(expected_data_dirname, "case.json")
    if os.path.exists(json_filename):
        os.remove(json_filename)
    return expected_json_filename, json_filename

def simulation_properties():
    initial_time = 1.2e8*365.25 # time [days] where simulation starts
    time_step = 0.08 # days
    time_limit   = time_step*200. # days
    historic_snapshot_period = 100.*365.25 # days
    recovery_snapshot_period = 10.*historic_snapshot_period # days
    return initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period
