from tools import interpolate_b_spline
import numpy as np
from constants import *

class EvolutionType(object):
    def __init__(self, variant, subvariant=None, mass=None):
        self._data = {}
        if variant == "BrownDwarf":
            self._data['fields'] = [float(mass)]
            self._data['variant'] = variant
        elif variant == "SolarLike":
            self._data['variant'] = variant
            if subvariant == "EvolvingDissipation":
                self._data['fields'] = [{ "variant": subvariant, "fields": [float(mass)] }]
            else:
                # ConstantDissipation
                self._data['fields'] = [subvariant]
            pass
        elif variant in ("MDwarf", "Jupiter", "NonEvolving"):
            self._data = variant

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class BrownDwarf(EvolutionType):
    def __init__(self, mass):
        super(BrownDwarf, self).__init__("BrownDwarf", mass=mass)

    def get_evolver(self, initial_time):
        mass = self._data['fields'][0]
        if mass <= 0.0101 and mass >= 0.0099:
            filename = "input/data_host_body/mass_10.0000.dat"
        elif mass <= 0.0121 and mass >= 0.0119:
            filename = "input/data_host_body/mass_12.0000.dat"
        elif mass <= 0.0151 and mass >= 0.0149:
            filename = "input/data_host_body/mass_15.0000.dat"
        elif mass <= 0.0201 and mass >= 0.0199:
            filename = "input/data_host_body/mass_20.0000.dat"
        elif mass <= 0.0301 and mass >= 0.0299:
            filename = "input/data_host_body/mass_30.0000.dat"
        elif mass <= 0.0401 and mass >= 0.0399:
            filename = "input/data_host_body/mass_40.0000.dat"
        elif mass <= 0.0501 and mass >= 0.0499:
            filename = "input/data_host_body/mass_50.0000.dat"
        elif mass <= 0.0601 and mass >= 0.0599:
            filename = "input/data_host_body/mass_60.0000.dat"
        elif mass <= 0.0701 and mass >= 0.0699:
            filename = "input/data_host_body/mass_70.0000.dat"
        elif mass <= 0.0721 and mass >= 0.0719:
            filename = "input/data_host_body/mass_72.0000.dat"
        elif mass <= 0.0751 and mass >= 0.0749:
            filename = "input/data_host_body/mass_75.0000.dat"
        elif mass <= 0.0801 and mass >= 0.0799:
            filename = "input/data_host_body/mass_80.0000.dat"
        else:
            raise Exception("The evolution type BrownDwarf does not support a mass of {} Msun!".format(mass))
        aux_filename = "input/data_host_body/rg2BD.dat"
        data = np.loadtxt(BASE_DIR+filename)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,1] * R_SUN
        aux_data = np.loadtxt(BASE_DIR+aux_filename)
        # BrownDwarf have a separate file for radius of gyration with a different time sampling
        # that should be homogenized:
        aux_time = aux_data[:, 0] * 365.25 - initial_time
        aux_radius_of_gyration_2 = aux_data[:, 2]
        resampled_radius_of_gyration_2 = []
        for current_time in time:
            current_radius_of_gyration_2, ignore = interpolate_b_spline(aux_time, aux_radius_of_gyration_2, current_time)
            resampled_radius_of_gyration_2.append(current_radius_of_gyration_2)

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = resampled_radius_of_gyration_2
        return evolver

class MDwarf(EvolutionType):
    def __init__(self):
        super(MDwarf, self).__init__("MDwarf")

    def get_evolver(self, initial_time):
        filename = "input/data_host_body/01Msun.dat"
        data = np.loadtxt(BASE_DIR+filename)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,1] * R_SUN

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = []
        return evolver

class Jupiter(EvolutionType):
    def __init__(self):
        super(Jupiter, self).__init__("Jupiter")

    def get_evolver(self, initial_time):
        filename = "input/data_host_body/Jupiter.dat"
        data = np.loadtxt(BASE_DIR+filename)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,1] * M2AU
        love_number = data[:,2]
        radius_of_gyration_2 = data[:,3]

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = love_number.tolist()
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = radius_of_gyration_2.tolist()
        return evolver

class SolarLikeEvolvingDissipation(EvolutionType):
    def __init__(self, mass):
        super(SolarLikeEvolvingDissipation, self).__init__("SolarLike", subvariant="EvolvingDissipation", mass=mass)

    def get_evolver(self, initial_time):
        if mass <= 0.401 and mass >= 0.399:
            filename = "input/data_host_body/L04Z02r.dat"
        elif mass <= 0.501 and mass >= 0.499:
            filename = "input/data_host_body/L05Z02r.dat"
        elif mass <= 0.601 and mass >= 0.599:
            filename = "input/data_host_body/L06Z02r.dat"
        elif mass <= 0.701 and mass >= 0.699:
            filename = "input/data_host_body/L07Z02r.dat"
        elif mass <= 0.801 and mass >= 0.799:
            filename = "input/data_host_body/L08Z02r.dat"
        elif mass <= 0.901 and mass >= 0.899:
            filename = "input/data_host_body/L09Z02r.dat"
        elif mass <= 1.001 and mass >= 0.999:
            filename = "input/data_host_body/L10Z02r.dat"
        elif mass <= 1.101 and mass >= 1.099:
            filename = "input/data_host_body/L11Z02r.dat"
        elif mass <= 1.201 and mass >= 1.199:
            filename = "input/data_host_body/L12Z02r.dat"
        elif mass <= 1.201 and mass >= 1.499:
            filename = "input/data_host_body/L13Z02r.dat"
        elif mass <= 1.301 and mass >= 1.299:
            filename = "input/data_host_body/L14Z02r.dat"
        elif mass <= 1.401 and mass >= 1.399:
            filename = "input/data_host_body/L15Z02r.dat"
        else:
            raise Exception("The evolution type MathisSolarLike does not support a mass of {} Msun!".format(mass))

        data = np.loadtxt(BASE_DIR+filename)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,1] * R_SUN
        inverse_tidal_q_factor = data[:,2]

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = inverse_tidal_q_factor.tolist()
        evolver['radius_of_gyration_2'] = []
        return evolver

class SolarLikeConstantDissipation(EvolutionType):
    def __init__(self, mass):
        super(SolarLikeEvolvingDissipation, self).__init__("SolarLike", subvariant="ConstantDissipation")

    def get_evolver(self, initial_time):
        filename = "input/data_host_body/SRad_Spli_M-1_0000.dat"
        data = np.loadtxt(BASE_DIR+filename)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,1] * M2AU

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = []
        return evolver

class NonEvolving(EvolutionType):
    def __init__(self):
        super(NonEvolving, self).__init__("NonEvolving")

    def get_evolver(self, initial_time):
        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = []
        evolver['time'] = []
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = []
        return evolver
