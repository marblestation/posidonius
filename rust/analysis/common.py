import os
import pandas as pd
import numpy as np
from fractions import Fraction
import struct


#-------------------------------------------------------------------------------
# Constants in S.I. units
#-------------------------------------------------------------------------------
Msun      =  1.98892e30               # kg
Rsun      =  6.96e8                   # m
AU        =  1.49598e11               # m
day       =  24*3600.                 # s
Rearth    =  6371.0e3                 # m

G         =  6.6742367e-11            # m^3.kg^-1.s^-2

yr        =  365.25*24*3600           # s
hr        =  3600.                    # s

Mjup      =  9.5511e-4 * Msun         # kg
Mearth    =  3.e-6 * Msun             # kg
Rjup      =  69173.e3                 # m

RAD2DEG = 180./np.pi



#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def read_history(filename):
    f = open(filename, "rb")
    # (np.floor(np.log10(np.max((100., 10.)))) - 2.)*10.

    if not os.path.exists(filename):
        raise Exception("File does not exists!")

    fields = ('current_time', 'time_step', 'particle', 'position_x', 'position_y', 'position_z', 'spin_x', 'spin_y', 'spin_z', 'velocity_x', 'velocity_y', 'velocity_z', 'semi-major_axis', 'perihelion_distance', 'eccentricity', 'inclination', 'longitude_of_perihelion', 'longitude_of_ascending_node', 'mean_anomaly', 'orbital_angular_momentum_x', 'orbital_angular_momentum_y', 'orbital_angular_momentum_z', 'orbital_angular_momentum', 'denergy_dt', 'total_energy', 'total_angular_momentum', 'mass', 'radius', 'radius_of_gyration_2', 'scaled_dissipation_factor', 'love_number', 'lag_angle')

    data = []
    while True:
        try:
            row = f.read(8+8+4+8*(len(fields)-3))
            vrow = struct.unpack('> d d i' + ' d'*(len(fields)-3), row)
        except:
            break
        else:
            data.append(vrow)

    data = pd.DataFrame(data, columns=fields, index=np.arange(len(data)))
    if len(data) == 0:
        raise Exception("Empty file!")

    # Force to always have N lines per snapshot corresponding to N particles
    n_particles = int(data['particle'].max())+1
    outer_particles = n_particles-1
    last_particle = int(data.iloc[-1]['particle'])
    excess = (n_particles - (outer_particles - last_particle)) % n_particles
    if excess > 0:
        data = data[:-1*excess]
    data = data.to_records()
    return n_particles, data

def filter_history(n_particles, data, discard_first_hundred_years=False):
    # Ignore first 100 years
    data['current_time'] /= 362.25 # From days to years
    if discard_first_hundred_years:
        data = data[data['current_time'] >= 100.]

    star = data['particle'] == 0
    star_data = data[star]

    planets_data = {}
    planets_keys = [] # To ensure the order
    for i in xrange(n_particles-1):
        planets_data["{}".format(i+1)] = data[data['particle'] == i+1]
        planets_keys.append("{}".format(i+1))

    return star_data, planets_data, planets_keys

#-------------------------------------------------------------------------------
# Main process
#-------------------------------------------------------------------------------
# The 2 eccentricity dependant factors in the equation in a
def Na1(e):
    return (1.0+31.0/2.0*(e*e)+255.0/8.0*np.power(e, 4) \
            + 185.0/16.0* np.power(e, 6) + 25.0/64.0* np.power(e,8)) / np.power((1.0- (e*e)), (15.0/2.0))

def Na2(e):
    return (1.+15./2.0*(e*e) + 45./8.0 * np.power(e, 4) + 5./16.0* np.power(e, 6)) / np.power((1.0- (e*e)), 6)

def No2(e):
    return (1.0 + 3.0 * (e*e) + 3.0/8.0 * np.power(e, 4)) / np.power((1.0- (e*e)), 5)

def energydot(a, e, rotp, oblp, G, Mp, Ms, Rp, k2deltat_plan):
    return 2. * Kplan(k2deltat_plan, G, Mp, Ms, Rp, a) \
            * (Na1(e) - 2.0*Na2(e) * np.cos(oblp) * (rotp/(norb(G, Mp, Ms)*np.power(a, -1.5))) \
            + (1.0 + np.power(np.cos(oblp),2))/2.0 * No2(e) * np.sqrt(1.0-(e*e)) * np.power(rotp/(norb(G, Mp, Ms)* np.power(a, -1.5)), 2))

# Jeremys Ki factor :
def Kplan(k2deltat_plan, G, Mp, Ms, Rp, a):
    return 3./2. * k2deltat_plan * (G*(Mp*Mp)/Rp) * np.power(Ms/Mp,2) \
            * np.power(Rp/a, 6) * np.power(norb(G,Mp,Ms) * np.power(a, -1.5), 2)

# Mean orbital angular velocity without the a dependance         (m^3/2.s-1)
def norb(G, Mp, Ms):
    return np.sqrt(G) * np.sqrt(Mp+Ms)

#-------------------------------------------------------------------------------
# Resonances
#-------------------------------------------------------------------------------


def get_possible_resonances(periodRatio, uncertainty=0.05, denominator_limit=12, numerator_limit=20, sampling=10):
    """
    Give a list of 'Fraction' objects that correspond to possible Mean Motion Resonances for a given period ratio.
    Optional argument are the uncertainty that will determine the range of perio ratios to test around the given value.
    The denominator and numerator limit above which all the resonances will be skipped

    Parameter :
    periodRatio : [float] the periodRatio between the two considered parameters

    Optional parameters :
    uncertainty=0.05 : the period range to test will be [p*(1-0.05) ; p * (1+0.05)]
    denominator_limit=12 : 13:12 resonance will be ok, but 14:13 resonance will be skipped
    numerator_limit=20 : 20:12 resonance will be ok, but 21:12 resonance will be skipped
    sampling=10 : The number of period ratio values for the range of values of period ratio to be tested around the nominal value
                  (given the uncertainty)

    Return :
    list of 'Fraction' objects, each one representing a Mean Motion Resonance to be tested

    NOTE: From a package developed by Autiwa <autiwa@gmail.com>
    """

    periodMin = periodRatio * (1 - uncertainty)
    periodMax = periodRatio * (1 + uncertainty)

    # We do not want period ratios less than 1 (this only happens for coorbitals I think)
    if (periodMin < 1.):
        periodMin = 1.

    periodWidth = periodMax - periodMin
    deltaPeriod = periodWidth/sampling


    periods = [periodMin + deltaPeriod * i for i in range(sampling)]

    resonances = []
    for period_i in periods:
        fraction = Fraction(period_i).limit_denominator(denominator_limit)
        resonances.append(fraction)

    # We exclude all values that appears several time to only keep one occurence of each value
    resonances = list(set(resonances))

    # We sort the resonances to get the more interesting first (3:2 before 32:27 for instance)
    tmp = [(res.numerator, res) for res in resonances]
    tmp.sort()
    resonances = [element[1] for element in tmp if element[1].numerator < numerator_limit]

    #~ print(uncertainty)
    #~ print(periodMin)
    #~ print(periodMax)
    #~ print(resonances)
    #~ exit()

    return resonances

def isResonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer, nb_points=50, angle_center_value=0, std_threshold=20.):
    """
    Given a resonance as a Fraction object, and g, n M for inner and
    outer planet, the function return if there is the resonance between
    the two planets

    Parameters :
    res : a Fraction object (for instance Fraction(3,2))
    g_inner, n_inner, M_inner : g, n, M for the inner planet (in degrees)
    g_outer, n_outer, M_outer : g, n, M for the outer planet (in degrees)

    Optional parameters :
    nb_points : [50] The number of points, at the end the each time span, that we will use to test resonance
    angler_center_value : [0.] (in degrees) we will overlap all the angle, through congruence, to have angle centered over this value.
                          For instance, with 0, angle will be between -180 and +180 degrees.
    std_threshold = [20.] in degrees, the value below which we will consider that a resonant angle prove the existence of a MMR.

    Return :
    True or False

    NOTE: From a package developed by Autiwa <autiwa@gmail.com>
    """
    outer_period_nb = res.denominator
    inner_period_nb = res.numerator

    angle_min = angle_center_value - 180.
    angle_max = angle_center_value + 180.

    # Resonances are usually displayed as (p+q):p where q is the order of
    # the resonance. We retreive thoses parameters
    p = outer_period_nb
    q = inner_period_nb - outer_period_nb

    # We calculate the resonant angles
    long_of_peri_inner = g_inner + n_inner
    mean_longitude_inner = M_inner + long_of_peri_inner

    long_of_peri_outer = g_outer + n_outer
    mean_longitude_outer = M_outer + long_of_peri_outer

    phi = np.empty((q+1, nb_points)) # we create an array, empty for the moment, that will contain all the resonant angles associated with the supposed resonance.

    temp_value = inner_period_nb * mean_longitude_outer - outer_period_nb * mean_longitude_inner

    for i in range(q+1):
        phi[i] = temp_value - i * long_of_peri_inner - (q - i) * long_of_peri_outer

    # We take modulo 2*pi of the resonant angle
    phi = phi%(360.)
    too_low = phi < angle_min
    too_high = phi > angle_max
    phi[too_low] = phi[too_low] + 360.
    phi[too_high] = phi[too_high] - 360.

    delta_longitude = long_of_peri_outer - long_of_peri_inner
    delta_longitude = delta_longitude%(360.)
    too_low = delta_longitude < angle_min
    too_high = delta_longitude > angle_max
    delta_longitude[too_low] = delta_longitude[too_low] + 360.
    delta_longitude[too_high] = delta_longitude[too_high] - 360.

    # If one of the std's is small (Typically, around 25, but
    # I had once a 80 that was not a resonance, so I think a threshold around 40 is a good one)
    standard_deviation = min(phi.std(1))

    if (standard_deviation < std_threshold):
        return True
    else:
        return False

def calculate_resonance_data(planets_keys, planets_data):
    # We initialize the arrays
    t = [] # time in years
    a = [] # demi-grand axe en ua
    e = [] # eccentricity
    g = [] # argument of pericentre or argument of perihelion (degrees)
    n = [] # longitude of ascending node (degrees)
    M = [] # Mean anomaly (degrees)
    q = [] # periastron or perihelion (AU)
    Q = [] # apastron or aphelion (AU)

    for key in planets_keys:
        planet_data = planets_data[key]
        t.append(planet_data['current_time']) # Years
        a.append(planet_data['semi-major_axis']) # AU
        e.append(planet_data['eccentricity'])
        # Change from longitude of perihelion to argument of perihelion
        longitude_of_perihelion_degrees = (planet_data['longitude_of_perihelion'] * RAD2DEG) % 360.
        longitude_of_ascending_node_degrees = (planet_data['longitude_of_ascending_node'] * RAD2DEG) % 360.
        g.append((longitude_of_perihelion_degrees - longitude_of_ascending_node_degrees) % 360.)
        n.append(longitude_of_ascending_node_degrees)
        M.append((planet_data['mean_anomaly'] * RAD2DEG) % 360.)
        #q.append(planet_data['perihelion_distance'])
        #Q.append(planet_data['semi-major_axis']*2 - planet_data['perihelion_distance'])
        q.append(planet_data['semi-major_axis'] * (1 - planet_data['eccentricity']))
        Q.append(planet_data['semi-major_axis'] * (1 + planet_data['eccentricity']))
    return t, a, e, g, n, M, q, Q

def get_subplot_shape(number_of_plots):
    """
    We chose the number of plot in x and y axis for the p.multi
    environment in order to plot ALL the resonant angles and have x and
    y numbers as close on from another as possible.

    HOW TO :
    fig = pl.figure(1)
    subplot_i = 0
    (nb_line, nb_row) = get_subplot_shape(5)
    subplot_i += 1
    plot_dof = fig.add_subplot(nb_line, nb_row, subplot_i)


    NOTE: From a package developed by Autiwa <autiwa@gmail.com>
    """
    nb_plots_x = 1
    nb_plots_y = 1
    while (nb_plots_x * nb_plots_y < number_of_plots):
        if (nb_plots_x == nb_plots_y):
            nb_plots_y += 1
        else:
            nb_plots_x += 1

    return (nb_plots_x, nb_plots_y)
