import numpy as np
import os

MAX_PARTICLES = 10
BASE_DIR = os.path.dirname(os.path.realpath(__file__)) + "/../"
if not os.path.exists(BASE_DIR+"/input/"):
    BASE_DIR = os.path.dirname(os.path.realpath(__file__)) + "/../../"
    if not os.path.exists(BASE_DIR+"/input/"):
        raise Exception("Input directory with stellar models not found!")

K = 0.01720209895 # Gaussian constant
K2 = K*K
G  = K2  # Gravitational constant with Mercury units
G_SI = 6.6742367e-11  # m^3.kg^-1.s^-2 (S.I. units)
PI = np.pi
TWO_PI = PI*2.
DEG2RAD = PI / 180. # conversion factor from degrees to radians

R_SUN = 4.67920694e-3 # AU
R_EARTH = 4.25874677e-5 # AU
AU = 1.49598e11 # m
HR = 3600. # s
M_EARTH = 3.0e-6 # M_SUN
M_SUN = 1.98892e30 # Kg

M2EARTH = 1.9891e6/5.9794 # Factor for mass-radius relation (valid only for earth type planets)
M2AU = 6.684587153547e-12 # 1 meter in AU
