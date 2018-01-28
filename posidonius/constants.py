import numpy as np
import os

MAX_PARTICLES = 10
BASE_DIR = os.path.dirname(os.path.realpath(__file__)) + "/../"
if not os.path.exists(BASE_DIR+"/input/"):
    BASE_DIR = os.path.dirname(os.path.realpath(__file__)) + "/../../"
    if not os.path.exists(BASE_DIR+"/input/"):
        raise Exception("Input directory with stellar models not found!")

PI = np.pi
TWO_PI = PI*2.
DEG2RAD = PI / 180. # conversion factor from degrees to radians
RAD2DEG = 180./PI

################################################################################
# The IAU 2009 system of astronomical constants: the report of the IAU working group on numerical standards for Fundamental Astronomy
# https://en.wikipedia.org/wiki/Astronomical_constant
M_SUN = 1.9818e30 # kg
M2EARTH = 332946.050895 # ratio sun/earth mass
#M_EARTH = M_SUN/M2EARTH # kg
M_EARTH = 1./M2EARTH # 3.0034895963231186e-06 M_SUN
G_SI = 6.67428e-11 # m^3.kg^-1.s^-2 (S.I. units)
AU = 1.49597870700e11 # m
M2AU = 1./AU # 6.684587122268445e-12 AU (1 meter in AU)

## Resolution B3 on recommended nominal conversion constants for selected solar and planetary properties
# http://adsabs.harvard.edu/abs/2015arXiv151007674M
R_SUN = 6.957e8 # m
R_SUN = R_SUN/AU # 4.650467260962157e-3 AU
R_EARTH = 6.3781e6 # m
R_EARTH = R_EARTH/AU # 4.263496512454037e-05 AU
################################################################################

#-------------------------------------------------------------------------------
# Constants in S.I. units
#-------------------------------------------------------------------------------
HOUR      =  3600.                    # s
DAY       =  24.*HOUR                 # s
YEAR      =  365.25*DAY               # s

#-------------------------------------------------------------------------------
# Constants in the units used in the Mercury code
#-------------------------------------------------------------------------------
G_MERCURY = G_SI/(np.power(AU,3))*M_SUN*(DAY*DAY) # Gravitational constant with Mercury units
G = G_MERCURY
K2 = G_MERCURY
K = np.sqrt(G_MERCURY) # Gaussian constant (in Mercury is 0.01720209895 but it does not follow from the recommended IAU constants)

