import datetime
from tools import interpolate_b_spline
import numpy as np
from constants import *

class EvolutionType(object):
    def __init__(self, variant, mass=None):
        self._data = {}
        if variant in ("BolmontMathis2016", "Baraffe2015", "Leconte2011", "Baraffe1998", "GalletBolmont2017"):
            self._data[variant] = float(mass)
        elif variant in ("LeconteChabrier2013", "NonEvolving"):
            self._data = variant
        else:
            raise Exception("Unknown variant '{}'".format(variant))

        if variant == "GalletBolmont2017":
            print("[WARNING {} UTC] Bodies with GalletBolmont2017 evolution will ignore initial radius and dissipation factor.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            print("[WARNING {} UTC] GalletBolmont2017 prescription theoretically only works for circular orbits and non inclined orbits, use carefully.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if variant == "BolmontMathis2016":
            print("[WARNING {} UTC] Bodies with BolmontMathis2016 evolution will ignore initial radius and dissipation factor.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            print("[WARNING {} UTC] BolmontMathis2016 prescription theoretically only works for circular orbits and non inclined orbits, use carefully.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        elif variant == "Baraffe2015":
            print("[WARNING {} UTC] Bodies with Baraffe2015 evolution will ignore initial radius and radius of gyration.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        elif variant == "Leconte2011":
            print("[WARNING {} UTC] Bodies with Leconte2011 evolution will ignore initial radius and radius of gyration.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        elif variant == "Baraffe1998":
            print("[WARNING {} UTC] Bodies with Baraffe1998 evolution will ignore initial radius.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        elif variant == "LeconteChabrier2013":
            print("[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number.".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()

class Leconte2011(EvolutionType):
    def __init__(self, mass):
        super(Leconte2011, self).__init__("Leconte2011", mass=mass)

    def get_evolver(self, initial_time):
        mass = self._data[self.__class__.__name__]
        if mass <= 0.0101 and mass >= 0.0099:
            filename = "input/Leconte_2011/mass_10.0000.dat"
        elif mass <= 0.0121 and mass >= 0.0119:
            filename = "input/Leconte_2011/mass_12.0000.dat"
        elif mass <= 0.0151 and mass >= 0.0149:
            filename = "input/Leconte_2011/mass_15.0000.dat"
        elif mass <= 0.0201 and mass >= 0.0199:
            filename = "input/Leconte_2011/mass_20.0000.dat"
        elif mass <= 0.0301 and mass >= 0.0299:
            filename = "input/Leconte_2011/mass_30.0000.dat"
        elif mass <= 0.0401 and mass >= 0.0399:
            filename = "input/Leconte_2011/mass_40.0000.dat"
        elif mass <= 0.0501 and mass >= 0.0499:
            filename = "input/Leconte_2011/mass_50.0000.dat"
        elif mass <= 0.0601 and mass >= 0.0599:
            filename = "input/Leconte_2011/mass_60.0000.dat"
        elif mass <= 0.0701 and mass >= 0.0699:
            filename = "input/Leconte_2011/mass_70.0000.dat"
        elif mass <= 0.0721 and mass >= 0.0719:
            filename = "input/Leconte_2011/mass_72.0000.dat"
        elif mass <= 0.0751 and mass >= 0.0749:
            filename = "input/Leconte_2011/mass_75.0000.dat"
        elif mass <= 0.0801 and mass >= 0.0799:
            filename = "input/Leconte_2011/mass_80.0000.dat"
        else:
            raise Exception("The evolution type Leconte2011 does not support a mass of {} Msun!".format(mass))
        aux_filename = "input/Leconte_2011/rg2BD.dat"
        data = np.loadtxt(BASE_DIR+filename)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,1] * R_SUN
        aux_data = np.loadtxt(BASE_DIR+aux_filename)
        # Leconte2011 has a separate file for radius of gyration with a different time sampling
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


class Baraffe2015(EvolutionType):
    """
    Constant dissipation
    """
    def __init__(self, mass):
        super(Baraffe2015, self).__init__("Baraffe2015", mass=mass)

    def get_evolver(self, initial_time):
        mass = self._data[self.__class__.__name__]
        if mass <= 0.0101 and mass >= 0.0099:
            filename = "input/Baraffe_2015/0010_Msun.dat"
        elif mass <= 0.0151 and mass >= 0.0149:
            filename = "input/Baraffe_2015/0015_Msun.dat"
        elif mass <= 0.0201 and mass >= 0.0199:
            filename = "input/Baraffe_2015/0020_Msun.dat"
        elif mass <= 0.0301 and mass >= 0.0299:
            filename = "input/Baraffe_2015/0030_Msun.dat"
        elif mass <= 0.0401 and mass >= 0.0399:
            filename = "input/Baraffe_2015/0040_Msun.dat"
        elif mass <= 0.0501 and mass >= 0.0499:
            filename = "input/Baraffe_2015/0050_Msun.dat"
        elif mass <= 0.0601 and mass >= 0.0599:
            filename = "input/Baraffe_2015/0060_Msun.dat"
        elif mass <= 0.0701 and mass >= 0.0699:
            filename = "input/Baraffe_2015/0070_Msun.dat"
        elif mass <= 0.0721 and mass >= 0.0719:
            filename = "input/Baraffe_2015/0072_Msun.dat"
        elif mass <= 0.0751 and mass >= 0.0749:
            filename = "input/Baraffe_2015/0075_Msun.dat"
        elif mass <= 0.0801 and mass >= 0.0799:
            filename = "input/Baraffe_2015/0080_Msun.dat"
        elif np.abs(mass - 0.09) < 1e-7:
            filename = "input/Baraffe_2015/0090_Msun.dat"
        elif np.abs(mass - 0.11) < 1e-7:
            filename = "input/Baraffe_2015/0110_Msun.dat"
        elif np.abs(mass - 0.13) < 1e-7:
            filename = "input/Baraffe_2015/0130_Msun.dat"
        elif np.abs(mass - 0.15) < 1e-7:
            filename = "input/Baraffe_2015/0150_Msun.dat"
        elif np.abs(mass - 0.17) < 1e-7:
            filename = "input/Baraffe_2015/0170_Msun.dat"
        elif np.abs(mass - 0.20) < 1e-7:
            filename = "input/Baraffe_2015/0200_Msun.dat"
        elif np.abs(mass - 0.30) < 1e-7:
            filename = "input/Baraffe_2015/0300_Msun.dat"
        elif np.abs(mass - 0.40) < 1e-7:
            filename = "input/Baraffe_2015/0400_Msun.dat"
        elif np.abs(mass - 0.50) < 1e-7:
            filename = "input/Baraffe_2015/0500_Msun.dat"
        elif np.abs(mass - 0.60) < 1e-7:
            filename = "input/Baraffe_2015/0600_Msun.dat"
        elif np.abs(mass - 0.70) < 1e-7:
            filename = "input/Baraffe_2015/0700_Msun.dat"
        elif np.abs(mass - 0.80) < 1e-7:
            filename = "input/Baraffe_2015/0800_Msun.dat"
        elif np.abs(mass - 0.90) < 1e-7:
            filename = "input/Baraffe_2015/0900_Msun.dat"
        elif np.abs(mass - 1.00) < 1e-7:
            filename = "input/Baraffe_2015/1000_Msun.dat"
        elif np.abs(mass - 1.10) < 1e-7:
            filename = "input/Baraffe_2015/1100_Msun.dat"
        elif np.abs(mass - 1.20) < 1e-7:
            filename = "input/Baraffe_2015/1200_Msun.dat"
        elif np.abs(mass - 1.30) < 1e-7:
            filename = "input/Baraffe_2015/1300_Msun.dat"
        elif np.abs(mass - 1.40) < 1e-7:
            filename = "input/Baraffe_2015/1400_Msun.dat"
        else:
            raise Exception("The evolution type Baraffe2015 does not support a mass of {} Msun!".format(mass))
        data = np.loadtxt(BASE_DIR+filename, skiprows=1)
        time = data[:,0] * 365.25 - initial_time
        radius = data[:,2] * R_SUN
        radius_of_gyration_2 = data[:,3]

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = radius_of_gyration_2.tolist()
        return evolver


class Baraffe1998(EvolutionType):
    def __init__(self, mass):
        super(Baraffe1998, self).__init__("Baraffe1998", mass=mass)

    def get_evolver(self, initial_time):
        mass = self._data[self.__class__.__name__]
        if np.abs(mass - 0.10) < 1e-7:
            # M Dwarf
            filename = "input/Baraffe_1998/01Msun.dat"
            data = np.loadtxt(BASE_DIR+filename)
            time = data[:,0] * 365.25 - initial_time
            radius = data[:,1] * R_SUN
        elif np.abs(mass - 1.00) < 1e-7:
            # Sun
            filename = "input/Baraffe_1998/SRad_Spli_M-1_0000.dat"
            data = np.loadtxt(BASE_DIR+filename)
            time = data[:,0] * 365.25 - initial_time
            radius = data[:,1] * M2AU
        else:
            raise Exception("The evolution type Baraffe1998 does not support a mass of {} Msun!".format(mass))

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = []
        evolver['radius_of_gyration_2'] = []
        return evolver

class LeconteChabrier2013(EvolutionType):
    """
    Jupiter type
    """
    def __init__(self):
        super(LeconteChabrier2013, self).__init__("LeconteChabrier2013")

    def get_evolver(self, initial_time):
        filename = "input/Leconte_Chabrier_2013/Jupiter.dat"
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

class BolmontMathis2016(EvolutionType):
    """
    Evolving dissipation
    """
    def __init__(self, mass):
        super(BolmontMathis2016, self).__init__("BolmontMathis2016", mass=mass)

    def get_evolver(self, initial_time):
        mass = self._data[self.__class__.__name__]
        if mass <= 0.401 and mass >= 0.399:
            filename = "input/Bolmont_Mathis_2016/L04Z02r.dat"
        elif mass <= 0.501 and mass >= 0.499:
            filename = "input/Bolmont_Mathis_2016/L05Z02r.dat"
        elif mass <= 0.601 and mass >= 0.599:
            filename = "input/Bolmont_Mathis_2016/L06Z02r.dat"
        elif mass <= 0.701 and mass >= 0.699:
            filename = "input/Bolmont_Mathis_2016/L07Z02r.dat"
        elif mass <= 0.801 and mass >= 0.799:
            filename = "input/Bolmont_Mathis_2016/L08Z02r.dat"
        elif mass <= 0.901 and mass >= 0.899:
            filename = "input/Bolmont_Mathis_2016/L09Z02r.dat"
        elif mass <= 1.001 and mass >= 0.999:
            filename = "input/Bolmont_Mathis_2016/L10Z02r.dat"
        elif mass <= 1.101 and mass >= 1.099:
            filename = "input/Bolmont_Mathis_2016/L11Z02r.dat"
        elif mass <= 1.201 and mass >= 1.199:
            filename = "input/Bolmont_Mathis_2016/L12Z02r.dat"
        elif mass <= 1.301 and mass >= 1.299:
            filename = "input/Bolmont_Mathis_2016/L13Z02r.dat"
        elif mass <= 1.401 and mass >= 1.399:
            filename = "input/Bolmont_Mathis_2016/L14Z02r.dat"
        else:
            raise Exception("The evolution type Bolmont_Mathis_2016 does not support a mass of {} Msun!".format(mass))

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


class GalletBolmont2017(EvolutionType):
    """
    Evolving dissipation
    """
    def __init__(self, mass):
        super(GalletBolmont2017, self).__init__("GalletBolmont2017", mass=mass)

    def get_evolver(self, initial_time):
        mass = self._data[self.__class__.__name__]
        if mass <= 0.301 and mass >= 0.299:
            filename = "input/Gallet_Bolmont_2017/M_03_Z_0134.dat"
        elif mass <= 0.401 and mass >= 0.399:
            filename = "input/Gallet_Bolmont_2017/M_04_Z_0134.dat"
        elif mass <= 0.601 and mass >= 0.599:
            filename = "input/Gallet_Bolmont_2017/M_06_Z_0134.dat"
        elif mass <= 0.701 and mass >= 0.699:
            filename = "input/Gallet_Bolmont_2017/M_07_Z_0134.dat"
        elif mass <= 0.801 and mass >= 0.799:
            filename = "input/Gallet_Bolmont_2017/M_08_Z_0134.dat"
        elif mass <= 0.901 and mass >= 0.899:
            filename = "input/Gallet_Bolmont_2017/M_09_Z_0134.dat"
        elif mass <= 1.001 and mass >= 0.999:
            filename = "input/Gallet_Bolmont_2017/M_10_Z_0134.dat"
        elif mass <= 1.101 and mass >= 1.099:
            filename = "input/Gallet_Bolmont_2017/M_11_Z_0134.dat"
        elif mass <= 1.201 and mass >= 1.199:
            filename = "input/Gallet_Bolmont_2017/M_12_Z_0134.dat"
        elif mass <= 1.401 and mass >= 1.399:
            filename = "input/Gallet_Bolmont_2017/M_14_Z_0134.dat"
        else:
            raise Exception("The evolution type Gallet_Bolmont_2017 does not support a mass of {} Msun!".format(mass))

        data = np.loadtxt(BASE_DIR+filename)
        time = np.power(10., data[:,0]) * 365.25 - initial_time
        radius = data[:,3] * R_SUN
        inverse_tidal_q_factor = 1./np.power(10., data[:,10])

        evolver = {}
        evolver['left_index'] = 0
        evolver['evolution_type'] = self.get()
        evolver['love_number'] = []
        evolver['radius'] = radius.tolist()
        evolver['time'] = time.tolist()
        evolver['inverse_tidal_q_factor'] = inverse_tidal_q_factor.tolist()
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
