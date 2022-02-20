import datetime
import numpy as np
from posidonius.constants import *
from posidonius.particles.axes import Axes

def _calculate_eccentricity_vector(reduced_masse, position, velocity):

    x = position.x()
    y = position.y()
    z = position.z()
    u = velocity.x()
    v = velocity.y()
    w = velocity.z()

    hx = y * w  -  z * v
    hy = z * u  -  x * w
    hz = x * v  -  y * u

    v_vect_h_x = v * hz - w * hy
    v_vect_h_y = w * hx - u * hz
    v_vect_h_z = u * hy - v * hx

    r = np.sqrt(x*x + y*y + z*z)

    ex = (v_vect_h_x / reduced_masse) - (x / r)
    ey = (v_vect_h_y / reduced_masse) - (y / r)
    ez = (v_vect_h_z / reduced_masse) - (z / r)

    return (ex, ey, ez)

def _find_indices_around_target_value(target_x, x):
    positions = np.where(x >= np.asarray(target_x))[0]
    if len(positions) > 0:
        right = positions[0]
        if right > 0:
            left = right - 1
        else:
            left = right
    else:
        left = len(x) - 1
        right = left
    return left, right

def linear_interpolation(target_x, x, y):
    """Linear interpolation"""
    # Find index
    left, right = _find_indices_around_target_value(target_x, x)
    if left == right:
        # Target value out of range, use limit values
        target_y = y[left]
    else:
        ### Interpolate
        ## Linear
        #target_y = (y[left] * (x[right] - target_x) + y[right] * (target_x - x[left])) / (x[right] - x[left])
        ## Linear (alternative)
        target_percent = (target_x - x[left])/(x[right] - x[left]) # Transform target to percent as in transforming x[left]..x[right] to 0..1
        target_y = y[left] * (1. - target_percent) + y[right] * target_percent
    return target_y

def cosine_interpolation(target_x, x, y):
    """Cosine interpolation"""
    # Find index
    left, right = _find_indices_around_target_value(target_x, x)
    if left == right:
        # Target value out of range, use limit values
        target_y = y[left]
    else:
        ### Interpolate
        ## Cosine interpolate (http://paulbourke.net/miscellaneous/interpolation/)
        # - Smooth around the real data points (contrary to the linear interpolation)
        target_percent = (target_x - x[left])/(x[right] - x[left]) # Transform target to percent as in transforming x[left]..x[right] to 0..1
        target_percent = (1. - math.cos(target_percent*PI)) / 2. # Transform target percent so that it gets smoothed when close to 0 or 1 (i.e., closer to x[left] or x[right])
        target_y = y[left] * (1. - target_percent) + y[right] * target_percent
    return target_y


def calculate_center_of_mass(masses, positions, velocities):
    if not isinstance(masses, list):
        masses = [masses]
    if not isinstance(positions, list):
        positions = [positions]
    if not isinstance(velocities, list):
        velocities = [velocities]
    if len(masses) != len(positions):
        raise Exception("The number of masses needs to match the number of positions")
    elif len(masses) != len(velocities):
        raise Exception("The number of masses needs to match the number of velocities")
    elif len(masses) == 0:
        raise Exception("At least there should be one mass")

    if isinstance(masses[0], np.ndarray):
        # Case where the first dimension corresponds to different particles
        # and the second dimension to different moments in time of that object
        #
        # masses = [(1., 1.), (0.8, 0.8)]
        # positions = [Axes(...), (...)]
        # velocities = [Axes(...), (...)]
        n = len(masses[0])
        if n != len(positions[0].x()):
            raise Exception("The number of snapshots for masses needs to match the number of snapshots for positions")
        elif n != len(velocities[0].x()):
            raise Exception("The number of snapshots for masses needs to match the number of snapshots for velocities")
        center_of_mass_position = Axes(np.zeros(n), np.zeros(n), np.zeros(n))
        center_of_mass_velocity = Axes(np.zeros(n), np.zeros(n), np.zeros(n))
        center_of_mass_mass = np.zeros(n)
    else:
        center_of_mass_position = Axes(0., 0., 0.)
        center_of_mass_velocity = Axes(0., 0., 0.)
        center_of_mass_mass = 0.

    for particle_mass, particle_position, particle_velocity in zip(masses, positions, velocities):
        center_of_mass_position.set_x(center_of_mass_position.x()*center_of_mass_mass + particle_position.x()*particle_mass)
        center_of_mass_position.set_y(center_of_mass_position.y()*center_of_mass_mass + particle_position.y()*particle_mass)
        center_of_mass_position.set_z(center_of_mass_position.z()*center_of_mass_mass + particle_position.z()*particle_mass)
        center_of_mass_velocity.set_x(center_of_mass_velocity.x()*center_of_mass_mass + particle_velocity.x()*particle_mass)
        center_of_mass_velocity.set_y(center_of_mass_velocity.y()*center_of_mass_mass + particle_velocity.y()*particle_mass)
        center_of_mass_velocity.set_z(center_of_mass_velocity.z()*center_of_mass_mass + particle_velocity.z()*particle_mass)

        center_of_mass_mass = center_of_mass_mass + particle_mass
        # Assumption: no zero zero mass
        center_of_mass_position.set_x(center_of_mass_position.x() / center_of_mass_mass)
        center_of_mass_position.set_y(center_of_mass_position.y() / center_of_mass_mass)
        center_of_mass_position.set_z(center_of_mass_position.z() / center_of_mass_mass)
        center_of_mass_velocity.set_x(center_of_mass_velocity.x() / center_of_mass_mass)
        center_of_mass_velocity.set_y(center_of_mass_velocity.y() / center_of_mass_mass)
        center_of_mass_velocity.set_z(center_of_mass_velocity.z() / center_of_mass_mass)
    return center_of_mass_mass, center_of_mass_position, center_of_mass_velocity

def mass_radius_relation(planet_mass, planet_mass_type='factor', planet_percent_rock = 0.70):
    """
    For earth-like planets, calculate the planet radius factor given the planet
    mass and the percentage of rock.

    - planet_mass_type='factor': Planet mass is a factor of R_EARTH (recommended)
    - planet_mass_type='AU': Planet mass is in AU
    """
    # Planetary radius in AU (rearth in AU) Rocky planet
    #   Earth-like => mass-radius relationship from Fortney 2007
    if planet_mass_type == 'factor':
        planet_mass_factor = planet_mass
        planet_radius_factor = (0.0592*planet_percent_rock + 0.0975) * np.power(np.log10(planet_mass_factor), 2) \
                                + (0.2337*planet_percent_rock + 0.4938) * (np.log10(planet_mass_factor)) \
                                + 0.3102*planet_percent_rock + 0.7932
    elif planet_mass_type == 'AU':
        planet_radius_factor = (0.0592*planet_percent_rock+0.0975) * np.power(np.log10(planet_mass) + np.log10(M2EARTH), 2) \
                                 + (0.2337*planet_percent_rock+0.4938) * (np.log10(planet_mass) + np.log10(M2EARTH)) \
                                 + 0.3102*planet_percent_rock+0.7932
    else:
        raise Exception("Unknown planet mass type: {}".format(planet_mass_type))
    return planet_radius_factor # Factor to be multiplied by R_EARTH

def calculate_keplerian_orbital_elements(target_mass, target_position, target_velocity, masses=[], positions=[], velocities=[]):
    """
    Calculate the keplerian orbital elements

    - target_mass: mass (typically the target planet mass+central body mass)
    - target_position: position (Axes)
    - target_velocity: position (Axes)

    In case of instead of a central massive body there are two (i.e., circumbinary
    planet), then it is necessary to provide their masses, positions and velocities
    (as list of Axes objects) and target_mass must correspond only to the planet
    mass. Order is important, for instance:

    - target_mass: planet mass
    - masses: [star1_mass, star2_mass]
    - positions: [star1_position, star2_position]
    - velocities: [star1_velocity, star2_velocity]

    It returns:

    - a: semi-major axis (in AU)
    - q: perihelion distance
    - eccentricity: eccentricity
    - i: inclination
    - p: longitude of perihelion (NOT the argument of perihelion)
    - n: longitude of ascending node
    - l: mean anomaly (or mean longitude if eccentricity < 1.e-8)
    """

    if not isinstance(masses, list):
        masses = [masses]
    if not isinstance(positions, list):
        positions = [positions]
    if not isinstance(velocities, list):
        velocities = [velocities]
    if len(masses) != len(positions):
        raise Exception("The number of masses needs to match the number of positions")
    elif len(masses) != len(velocities):
        raise Exception("The number of masses needs to match the number of velocities")

    if len(positions) > 0:
        center_of_mass_mass, center_of_mass_position, center_of_mass_velocity = calculate_center_of_mass(masses, positions, velocities)
        target_position = Axes(target_position.x()-center_of_mass_position.x(), target_position.y()-center_of_mass_position.y(), target_position.z()-center_of_mass_position.z())
        target_velocity = Axes(target_velocity.x()-center_of_mass_velocity.x(), target_velocity.y()-center_of_mass_velocity.y(), target_velocity.z()-center_of_mass_velocity.z())

    gm = G*(sum(masses)+target_mass)
    a, q, eccentricity, i, p, n, l = _calculate_keplerian_orbital_elements(gm, target_position, target_velocity)

    if len(positions) > 0:
        delta_a = np.sqrt(center_of_mass_position.x()**2 + center_of_mass_position.y()**2 + center_of_mass_position.z()**2)
        delta_q = delta_a * (1.0 - eccentricity);                     # perihelion distance
        a += delta_a
        q += delta_q
    return (a, q, eccentricity, i, p, n, l)

def _calculate_keplerian_orbital_elements(gm, position, velocity):
    # ! Based on the implementation of Chambers in Mercury
    # Calculates Keplerian orbital elements given relative coordinates and
    # velocities, and GM = G times the sum of the masses.

    # Input
    x = position.x()
    y = position.y()
    z = position.z()
    u = velocity.x()
    v = velocity.y()
    w = velocity.z()
    gm = np.asarray((gm,)) if type(x) is float else gm
    x = np.asarray((x,)) if type(x) is float else x
    y = np.asarray((y,)) if type(y) is float else y
    z = np.asarray((z,)) if type(z) is float else z
    u = np.asarray((u,)) if type(u) is float else u
    v = np.asarray((v,)) if type(v) is float else v
    w = np.asarray((w,)) if type(w) is float else w
    # --- Output
    # a # semi-major axis (in AU)
    # q # perihelion distance
    # ecc # eccentricity
    # i = 0. # inclination
    # p # longitude of perihelion (NOT argument of perihelion!!)
    # n # longitude of ascending node
    # l # mean anomaly (or mean longitude if eccentricity < 1.e-8)
    # Local
    ex, ey, ez = _calculate_eccentricity_vector(gm, position, velocity)
    e_scal_r = ex*x + ey*y + ez*z
    hx = y * w  -  z * v
    hy = z * u  -  x * w
    hz = x * v  -  y * u
    h2 = np.power(hx, 2) + np.power(hy, 2) + np.power(hz, 2)
    v2 = u * u  +  v * v  +  w * w
    r_scal_v = x * u  +  y * v  +  z * w
    r = np.sqrt(x*x + y*y + z*z)
    h = np.sqrt(h2)
    s = h2 / gm

    ###-------------------------------------------------------------------------
    ### a: Semi-major axis
    a  = gm * r / (2.0 * gm  -  r * v2)


    # Inclination and node
    ci = hz / h
    ci = np.clip(ci, -1., 1.)

    ###-------------------------------------------------------------------------
    ### n: longitude_of_ascending_node
    ### i: inclination
    # Algorithm:
    #if abs(ci) < 1.:
        #i = np.arccos(ci)
        #n = np.arctan2(hx, -hy)
        #if n < 0.:
            #n = n + TWO_PI
    #else:
        #if ci > 0.:
            #i = 0.
        #if ci < 0.:
            #i = PI
        #n = 0.
    i = np.zeros(len(ci))
    n = np.zeros(len(ci))
    less_than_one = np.abs(ci) < 1.
    if np.any(less_than_one):
        i[less_than_one] = np.arccos(ci[less_than_one])
        #n[less_than_one] = np.arctan(hx[less_than_one] / -hy[less_than_one])
        n[less_than_one] = np.arctan2(hx[less_than_one], -hy[less_than_one])
        positive_hy = hy[less_than_one] > 0.
        if np.any(positive_hy):
            i[positive_hy] = TWO_PI - i[positive_hy]
    negative_n = np.logical_and(n < 0, less_than_one)
    if np.any(negative_n):
        n[negative_n] += TWO_PI

    rest = np.logical_not(less_than_one)
    if np.any(rest):
        positive_ci = np.logical_and(ci > 0., rest)
        negative_ci = np.logical_and(ci < 0., rest)
        if np.any(positive_ci):
            i[positive_ci] = 0.
        if np.any(negative_ci):
            i[negative_ci] = PI
        n[rest] = 0.
    no_inclination = i == 0.
    if np.any(no_inclination):
        n[no_inclination] = 0.

    ###-------------------------------------------------------------------------
    ### eccentricity
    ### q: perihelion_distance
    temp = 1.  +  s * (v2 / gm  -  2. / r)

    eccentricity = np.zeros(len(ci))
    ## Algorithm:
    #if temp <= 0.:
        #eccentricity = 0.
    #else:
        #eccentricity = np.sqrt(temp)
    positive_temp = temp > 0.
    if np.any(positive_temp):
        eccentricity[positive_temp] = np.sqrt(temp[positive_temp])
    q = s / (1. + eccentricity)
    small_eccentricity = eccentricity < 3.0e-8
    if np.any(small_eccentricity):
        eccentricity[small_eccentricity] = 0.

    # if the eccentricity is calculated below 3.0e-8
    # we take the 0 eccentricity approximation
    null_eccentricity = eccentricity == 0.
    big_eccentricity = np.logical_not(null_eccentricity)

    ###-------------------------------------------------------------------------
    ### p: longitude_of_perihelion
    ### l: mean_anomaly
    # Algorithm:
    #if eccentricity < 3.0e-8:
        #p = 0.
    #else:
        ## Longitude of perihelion
        #cf = (s - r) / (eccentricity*r)
        #if abs(cf) > 1.:
            ##cf = cf.signum()
            #cf = np.sign(cf)

        #f = np.arccos(cf)
        #if r_scal_v < 0.:
            #f = TWO_PI - f

        #p = true_longitude - f
        #p = modulus((p + TWO_PI + TWO_PI), TWO_PI)
    arg_pericentre = np.zeros(len(ci))
    # The eccentricity vector through the rotation matrix R3
    ex_rot = np.zeros(len(ci))
    ey_rot = np.zeros(len(ci))
    ez_rot = np.zeros(len(ci))
    # The eccentricity vector through the rotation matrix R3 then the rotation matrix R
    ex_rot_2 = np.zeros(len(ci))
    ey_rot_2 = np.zeros(len(ci))
    ez_rot_2 = np.zeros(len(ci))
    if np.any(null_eccentricity):
        arg_pericentre[null_eccentricity] = 0.
    if np.any(big_eccentricity):
        # The argument of pericentre is calculated with theangle made by the eccentricity vector
        # with respect to the referential rotated first by the argument of ascending node with respect to the z axis
        # and then by the inclination woth repect to the mew x axis

        # Calculate the eccentricity vector through the rotation matrix R3
        # with respect to the longitude of ascending node
        ex_rot[big_eccentricity] = ex[big_eccentricity]*np.cos(n[big_eccentricity]) + ey[big_eccentricity]*np.sin(n[big_eccentricity])
        ey_rot[big_eccentricity] = -ex[big_eccentricity]*np.sin(n[big_eccentricity]) + ey[big_eccentricity]*np.cos(n[big_eccentricity])
        ez_rot[big_eccentricity] = ez[big_eccentricity]

        # Calculate the eccentricity vector (rotate with ) through the rotation matrix R1
        # with respect to the inclination
        ex_rot_2[big_eccentricity] = ex_rot[big_eccentricity]
        ey_rot_2[big_eccentricity] = ey_rot[big_eccentricity]*np.cos(i[big_eccentricity]) + ez_rot[big_eccentricity]*np.sin(i[big_eccentricity])
        ez_rot_2[big_eccentricity] = -ey_rot[big_eccentricity]*np.sin(i[big_eccentricity]) + ez_rot[big_eccentricity]*np.cos(i[big_eccentricity])

        arg_pericentre[big_eccentricity] = np.arctan2(ey_rot_2[big_eccentricity],ex_rot_2[big_eccentricity])

    p = np.zeros(len(eccentricity))
    p = arg_pericentre + n

    ###-------------------------------------------------------------------------
    ### True anomaly
    # Algorithm:
    #if hy != 0.:
        #to = -hx/hy
        #temp = (1. - ci) * to
        #tmp2 = to * to
        #true_longitude = np.arctan2(y*(1.+tmp2*ci)-x*temp, (x*(tmp2+ci)-y*temp))
    #else:
        #true_longitude = np.arctan2(y*ci, x)

    #if ci < 0.:
        #true_longitude = true_longitude + PI
    f = np.zeros(len(ci))
    cosf_cal = np.zeros(len(ci))
    if np.any(null_eccentricity):
        x_rot = np.zeros(len(ci))
        x_rot[null_eccentricity] = x[null_eccentricity]*np.cos(n[null_eccentricity]) + y[null_eccentricity]*np.sin(n[null_eccentricity])
        # y_rot = -x[null_eccentricity]*np.sin(n[null_eccentricity]) + y[null_eccentricity]*np.cos(n[null_eccentricity])
        # z_rot = z
        cosf_cal[null_eccentricity] = x_rot[null_eccentricity] / r[null_eccentricity]
        cosf_cal[null_eccentricity] = np.clip(cosf_cal[null_eccentricity], -1., 1.)
        f[null_eccentricity] = np.arccos(cosf_cal[null_eccentricity])
    if np.any(big_eccentricity):
        cosf_cal[big_eccentricity] = e_scal_r[big_eccentricity] / (eccentricity[big_eccentricity]*r[big_eccentricity])
        cosf_cal[big_eccentricity] = np.clip(cosf_cal[big_eccentricity], -1., 1.)
        ellipse = np.logical_and(big_eccentricity, eccentricity < 1.)

        f[big_eccentricity] = np.arccos(cosf_cal[big_eccentricity])
        negative_r_scal_v = r_scal_v < 0.
        if np.any(negative_r_scal_v):
            f[negative_r_scal_v] = TWO_PI -f[negative_r_scal_v]

        # print("True ano cal {} rad and {} deg".format( f, f/posidonius.constants.DEG2RAD))
        # cosf_cal = (1/e_cal)*( (a_cal*(1-e_cal*e_cal)/r_cal) -1.)
        # f = np.arccos(cosf_cal)
        # print("True ano cal {} rad and {} deg".format( f, f/posidonius.constants.DEG2RAD))

    ###-------------------------------------------------------------------------
    ### l: mean_anomaly
    # Algorithm:
    #if eccentricity < 3.0e-8:
        #l = true_longitude
    #else:
        #ce = (v2*r - gm) / (eccentricity*gm)

        ## Mean anomaly for ellipse
        #if eccentricity < 1.:
            #if abs(ce) > 1.:
                ##ce = ce.signum()
                #ce = np.sign(ce)

            #bige = np.arccos(ce)
            #if r_scal_v < 0.:
                #bige = TWO_PI - bige

            #l = bige - eccentricity*np.sin(bige)
        #else:
            ## Mean anomaly for hyperbola
            #if ce < 1.:
                #ce = 1.

            #bige = np.log(ce + np.sqrt(ce*ce-1.))
            #if r_scal_v < 0.:
                #bige = - bige

            #l = eccentricity * np.sinh(bige) - bige
    l = np.zeros(len(ci))
    if np.any(null_eccentricity): #mean_anomaly = true_longitude = true_anomaly + long_pericentre
        l[null_eccentricity] = f[null_eccentricity]

    if np.any(big_eccentricity):
        # ce = np.zeros(len(true_longitude))
        # ce[big_eccentricity] = (v2[big_eccentricity]*r[big_eccentricity] - gm[big_eccentricity]) / (eccentricity[big_eccentricity]*gm[big_eccentricity])
        ce = np.zeros(len(ci))
        ce[big_eccentricity] = (1./eccentricity[big_eccentricity])*(1.-(r[big_eccentricity]/a[big_eccentricity]))
        ellipse = np.logical_and(big_eccentricity, eccentricity < 1.)
        if np.any(ellipse):
            change_sign = np.logical_and(eccentricity, np.abs(ce) > 1.)
            if np.any(change_sign):
                ce[change_sign] = np.sign(ce[change_sign])
            bige = np.zeros(len(ce))
            ce[ellipse] = np.clip(ce[ellipse], -1., 1.)
            bige[ellipse] = np.arccos(ce[ellipse])
            negative_rv = np.logical_and(ellipse, r_scal_v < 0.)
            if np.any(negative_rv):
                bige[negative_rv] = TWO_PI - bige[negative_rv]
            l[ellipse] = bige[ellipse] - eccentricity[ellipse]*np.sin(bige[ellipse])

        hyperbola = np.logical_and(big_eccentricity, eccentricity >= 1.)
        if np.any(hyperbola):
            negative_ce = np.logical_and(hyperbola, ce < 1.)
            if np.any(negative_ce):
                ce[negative_ce] = 1.
            bige = np.zeros(len(ce))
            bige[hyperbola] = np.log(ce[hyperbola] + np.sqrt(ce[hyperbola]*ce[hyperbola]-1.))
            negative_rv = np.logical_and(hyperbola, r_scal_v < 0.)
            if np.any(negative_rv):
                bige[negative_rv] = - bige[negative_rv]
            l[hyperbola] = eccentricity[hyperbola] * np.sinh(bige[hyperbola]) - bige[hyperbola]

    # Algorithm:
    #if l < 0.:
        #l = l + TWO_PI

    #if l > TWO_PI:
        #l = modulus(l, TWO_PI)
    negative_l = l < 0.
    if np.any(negative_l):
        l[negative_l] += TWO_PI
    bigger_than_two_pi = l > TWO_PI
    if np.any(bigger_than_two_pi):
        l[bigger_than_two_pi] = modulus(l[bigger_than_two_pi], TWO_PI)

    ###-------------------------------------------------------------------------

    if len(a) == 1:
        return (a[0], q[0], eccentricity[0], i[0], p[0], n[0], l[0])
    else:
        return (a, q, eccentricity, i, p, n, l)

def calculate_cartesian_coordinates(target_mass, q, e, i0, p, n0, l, masses=[], positions=[], velocities=[]):
    """
    Calculate the heliocentric position and velocity

    - target_mass: mass (typically the target planet mass+central body mass)
    - q: perihelion distance
    - e: eccentricity
    - i0: inclination (degrees)
    - p: longitude of perihelion
    - n0: longitude of the ascending node (degrees)
    - l: mean anomaly

    In case of instead of a central massive body there are two (i.e., circumbinary
    planet), then it is necessary to provide their masses, positions and velocities
    (as list of Axes objects) and target_mass must correspond only to the planet
    mass. Order is important, for instance:

    - target_mass: planet mass
    - masses: [star1_mass, star2_mass]
    - positions: [star1_position, star2_position]
    - velocities: [star1_velocity, star2_velocity]
    """
    if not isinstance(masses, list):
        masses = [masses]
    if not isinstance(positions, list):
        positions = [positions]
    if not isinstance(velocities, list):
        velocities = [velocities]
    if len(masses) != len(positions):
        raise Exception("The number of masses needs to match the number of positions")
    elif len(masses) != len(velocities):
        raise Exception("The number of masses needs to match the number of velocities")

    if len(positions) > 0:
        center_of_mass_mass, center_of_mass_position, center_of_mass_velocity = calculate_center_of_mass(masses, positions, velocities)
        delta_a = np.sqrt(center_of_mass_position.x()**2 + center_of_mass_position.y()**2 + center_of_mass_position.z()**2)
        delta_q = delta_a * (1.0 - e);                     # perihelion distance
        q -= delta_q

    gm = G*(sum(masses)+target_mass)
    x, y, z, vx, vy, vz = _calculate_cartesian_coordinates(gm, q, e, i0, p, n0, l)

    if len(positions) > 0:
        x += center_of_mass_position.x()
        y += center_of_mass_position.y()
        z += center_of_mass_position.z()
        vx += center_of_mass_velocity.x()
        vy += center_of_mass_velocity.y()
        vz += center_of_mass_velocity.z()

    planet_position = Axes(x, y, z)
    planet_velocity = Axes(vx, vy, vz)
    return planet_position, planet_velocity


def _calculate_cartesian_coordinates(gm, q, e, i0, p, n0, l):
    # Calculates Cartesian coordinates and velocities given Keplerian orbital
    # elements (for elliptical, parabolic or hyperbolic orbits).
    # ! Based on the implementation of Chambers in Mercury, which is
    #   based on a routine from Levison and Duncan's SWIFT integrator.
    # WARNING: It gives NaN velocities when eccentricity == 1. (also in the original implementation in mercury code)

    # Input
    # gm  = grav const * (central + secondary mass)
    # q  = perihelion distance
    # e  = eccentricity
    # i  = inclination (degrees)
    # p  = longitude of perihelion !!!
    # n  = longitude of the ascending node (degrees)
    # l  = mean anomaly

    ## Output
    #x # Cartesian positions  ( units the same as a )
    #y # Cartesian positions  ( units the same as a )
    #z # Cartesian positions  ( units the same as a )
    #u # Cartesian velocities ( units the same as sqrt(gm/a) )
    #v # Cartesian velocities ( units the same as sqrt(gm/a) )
    #w # Cartesian velocities ( units the same as sqrt(gm/a) )

    # Change from longitude of perihelion to argument of perihelion
    g0 = p - n0

    # Rotation factors
    si = np.sin(i0)
    ci = np.cos(i0)
    sg = np.sin(g0)
    cg = np.cos(g0)
    sn = np.sin(n0)
    cn = np.cos(n0)
    #(i, si, ci) = mco_sine(i0)
    #(g, sg, cg) = mco_sine(g0)
    #(n, sn, cn) = mco_sine(n0)
    z1 = cg * cn
    z2 = cg * sn
    z3 = sg * cn
    z4 = sg * sn
    d11 =  z1 - z4*ci
    d12 =  z2 + z3*ci
    d13 = sg * si
    d21 = -z3 - z2*ci
    d22 = -z4 + z1*ci
    d23 = cg * si

    # Semi-major axis
    a = q / (1. - e)

    #se
    #ce

    # Ellipse
    if e < 1.:
       romes = np.sqrt(1. - e*e)
       temp0 = kepler_solution_for_eccentrities_smaller_than_one(e,l)
       se = np.sin(temp0)
       ce = np.cos(temp0)
       z1 = a * (ce - e)
       z2 = a * romes * se
       temp = np.sqrt(gm/a) / (1. - e*ce)
       z3 = -se * temp
       z4 = romes * ce * temp
    else:
       # Parabola
       if e == 1.:
          tmp = kepler_solution_for_a_parabola(l)
          ce = tmp[1]
          z1 = q * (1. - ce*ce)
          z2 = 2. * q * ce
          z4 = np.sqrt(2.*gm/q) / (1. + ce*ce)
          z3 = -ce * z4
       else:
          # Hyperbola
          romes = np.sqrt(e*e - 1.)
          temp = kepler_solution_for_a_hyperbola(e,l)
          se = np.sinh(temp)
          ce = np.cosh(temp)
          z1 = a * (ce - e)
          z2 = -a * romes * se
          temp = np.sqrt(gm/abs(a)) / (e*ce - 1.)
          z3 = -se * temp
          z4 = romes * ce * temp

    x = d11 * z1  +  d21 * z2
    y = d12 * z1  +  d22 * z2
    z = d13 * z1  +  d23 * z2
    u = d11 * z3  +  d21 * z4
    v = d12 * z3  +  d22 * z4
    w = d13 * z3  +  d23 * z4

    return (x, y, z, u, v, w)



def kepler_solution_for_eccentrities_smaller_than_one(e, oldl):
    # ! Based on the implementation of Chambers in Mercury (mco_kep)
    # Solves Kepler's equation for eccentricities less than one.
    # Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
    #
    # e = eccentricity\n
    # l = mean anomaly      (radians)\n
    # u = eccentric anomaly (   "   )\n

    piby2 = 0.5 * PI

    # Reduce mean anomaly to lie in the range 0 < l < pi
    if oldl >= 0.:
       l = modulus(oldl, TWO_PI)
    else:
       l = (modulus(oldl, TWO_PI)) + TWO_PI

    sign = 1.
    if l > PI:
       l = TWO_PI - l
       sign = -1.

    ome = 1. - e

    if (l >= 0.45) or (e < 0.55):
        # Regions A,B or C in Nijenhuis
        # -----------------------------

        # Rough starting value for eccentric anomaly
        if l < ome:
           u1 = ome
        else:
            if l > (PI-1.-e):
                u1 = (l+e*PI)/(1.+e)
            else:
                u1 = l + e

        # Improved value using Halley's method
        flag = u1 > piby2
        if flag:
           x = PI - u1
        else:
           x = u1

        x2 = x*x
        sn = x*(1. + x2*(-0.16605 + x2*0.00761) )
        dsn = 1. + x2*(-0.49815 + x2*0.03805)
        if flag:
            dsn = -dsn

        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1. - e*dsn
        u2 = u1 - f0/(f1 - 0.5*f0*f2/f1)
    else:
        # Region D in Nijenhuis
        # ---------------------

        # Rough starting value for eccentric anomaly
        z1 = 4.*e + 0.5
        p = ome / z1
        q = 0.5 * l / z1
        p2 = p*p
        z2 = np.exp(np.log(np.sqrt(p2*p + q*q ) + q) / 1.5)
        u1 = 2.*q / ( z2 + p + p2/z2 )

        # Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - 0.075*u1*z3 / (ome + z1*z2 + 0.375*z3)
        u2 = l + e*u2*( 3. - 4.*u2*u2 )


    # Accurate value using 3rd-order version of Newton's method
    # N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy

    # First get accurate values for u2 - sin(u2) and 1 - cos(u2)
    bigg = u2 > piby2
    if bigg:
        z3 = PI - u2
    else:
        z3 = u2


    big = z3 > (0.5*piby2)
    if big:
       x = piby2 - z3
    else:
       x = z3


    x2 = x*x

    ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
    cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - x2/306.))))))))

    if big:
       z1 = cc + z3 - 1.
       z2 = ss + z3 + 1. - piby2
    else:
       z1 = ss
       z2 = cc


    if bigg:
       z1 = 2.*u2 + z1 - PI
       z2 = 2. - z2


    f0 = l - u2*ome - e*z1
    f1 = ome + e*z2
    f2 = 0.5*e*(u2-z1)
    f3 = e/6.*(1.-z2)
    z1 = f0/f1
    z2 = f0/(f2*z1+f1)
    return sign*( u2 + f0/((f3*z1+f2)*z2+f1) )



def kepler_solution_for_a_parabola(q0):
    # ! Based on the implementation of Duncan in Mercury (orbel_zget)
    # Solves the equivalent of Kepler's eqn. for a parabola
    # given Q (Fitz. notation.)
    #
    # ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
    #
    # @remarks For a parabola we can solve analytically.

    # Input:
    # q0 = parabola mean anomaly

    # Output
    #eccentric_anomaly # eccentric anomaly

    # We copy the input parameter to be able to modify it.
    q = q0

    #----
    if q < 0.:
       iflag = 1
       q = -q
    else:
        iflag = 0


    if q < 1.0e-3:
       eccentric_anomaly = q*(1. - (q*q/3.)*(1. -q*q))
    else:
       x = 0.5*(3.*q + (9.* np.power(q, 2) + 4.).sqrt())
       tmp = np.power(x, 1./3.)
       eccentric_anomaly = tmp - 1./tmp


    if iflag == 1:
       eccentric_anomaly = -eccentric_anomaly
       q = -q

    return (q, eccentric_anomaly)


def kepler_solution_for_a_hyperbola(e, n):
    # ! Based on the implementation of Duncan in Mercury (orbel_fhybrid)
    # Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
    # ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON
    #            For larger N, uses FGET

    # Inputs
    # e = eccentircity
    # n = hyperbola mean anomaly

    # Output
    #orbel_fhybrid # eccentric anomaly

    #...  Internals:
    #abn

    #----
    #...  Executable code

    abn = n
    if n <0.:
        abn = -abn

    if abn < 0.636*e -0.6:
       tmp = kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(e,n)
       orbel_fhybrid = tmp[1]
    else:
       orbel_fhybrid = kepler_solution_for_a_hyperbola_hybrid_approach(e,n)


    return orbel_fhybrid



def kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(e, capn0):
    # ! Based on the implementation of Duncan in Mercury (orbel_flon)
    # Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
    # ALGORITHM: Uses power series for N in terms of F and Newton,s method
    # REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)

    # Inputs:
    # e = eccentricty
    # capn = hyperbola mean anomaly

    # Output
    #orbel_flon

    # copy of the input capn0 that is not modified
    capn = capn0

    imax = 10
    tiny = 4.0e-15  # A small number
    a3 = 1037836800.
    a5 = 51891840.
    a7 = 1235520.
    a9 = 17160.
    a11 = 156.

    b3 = 3.*a3
    b5 = 5.*a5
    b7 = 7.*a7
    b9 = 9.*a9
    b11 = 11.*a11

    # Function to solve "Kepler's eqn" for F (here called
    # x) for given e and CAPN. Only good for smallish CAPN

    iflag = 0
    if capn < 0.:
       iflag = 1
       capn = -capn


    a1 = 6227020800. * (1. - 1./e)
    a0 = -6227020800.*capn/e
    b1 = a1

    #  Set iflag nonzero if capn < 0., in which case solve for -capn
    #  and change the sign of the final answer for F.
    #  Begin with a reasonable guess based on solving the cubic for small F
    a = 6.*(e-1.)/e
    b = -6.*capn/e
    sq = np.sqrt(0.25*b*b +a*a*a/27.)
    biga = np.power(-0.5*b + sq, 0.3333333333333333)
    bigb = -np.power(0.5*b + sq, 0.3333333333333333)
    x = biga + bigb
    #  write(6,*) 'cubic = ',x**3 +a*x +b
    orbel_flon = x
    # If capn is tiny (or zero) no need to go further than cubic even for
    # e =1.
    if capn >= tiny:
        converge = False
        for ignore in np.arange(imax)+1:
           x2 = x*x
           f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
           fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.*x2)))))
           dx = -f/fp
           orbel_flon = x + dx
           #   If we have converged here there's no point in going on
           if dx.abs() < tiny:
                converge = True
                break

           x = orbel_flon

        if not converge:
            # Abnormal return here - we've gone thru the loop
            # IMAX times without convergence
            if iflag == 1:
               orbel_flon = -orbel_flon
               capn = -capn

            print("[WARNING {} UTC] FLON : RETURNING WITHOUT COMPLETE CONVERGENCE".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            diff = e * np.sinh(orbel_flon)  - orbel_flon - capn
            print("N, F, ecc * F.sinh() - F - N : ")
            print("{}, {}, {}".format(capn,orbel_flon,diff))
            return (orbel_flon, capn)

    #  Normal return here, but check if capn was originally negative
    if iflag == 1:
       orbel_flon = -orbel_flon
       capn = -capn

    return (orbel_flon, capn)



def kepler_solution_for_a_hyperbola_hybrid_approach(e, capn):
    # ! Based on the implementation of Duncan in Mercury (orbel_fget)
    # Solves Kepler's eqn. for hyperbola using hybrid approach.
    # ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
    #              Cel. Mech. ".  Quartic convergence from Danby's book.

    # Inputs:
    # e = eccentricty
    # capn = hyperbola mean anomaly

    # Output
    #orbel_fget # eccentric anomaly

    #...  Internals:
    imax = 10
    tiny = 4.0e-15  # A small number

    #----
    # Function to solve "Kepler's eqn" for F (here called
    # x) for given e and CAPN.

    #  begin with a guess proposed by Danby
    if capn < 0.:
       tmp = -2.*capn/e + 1.8
       x = - np.log(tmp)
    else:
       tmp = 2.*capn/e + 1.8
       x = np.log(tmp)


    orbel_fget = x

    for _ in range(imax):
       shx = np.sinh(x)
       chx = np.cosh(x)
       esh = e*shx
       ech = e*chx
       f = esh - x - capn
       fp = ech - 1.
       fpp = esh
       fppp = ech
       dx = -f/fp
       dx = -f/(fp + dx*fpp/2.)
       dx = -f/(fp + dx*fpp/2. + dx*dx*fppp/6.)
       orbel_fget = x + dx
       #   If we have converged here there's no point in going on
       if np.abs(dx) <= tiny:
        return orbel_fget

       x = orbel_fget


    print("[WARNING {} UTC] FGET : RETURNING WITHOUT COMPLETE CONVERGENCE".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
    return orbel_fget




def modulus(a, b):
    #println!("Modulus::", -21.0f64 % 4.0f64)         # -1 because -21 divided by 4 gives -5 with a remainder of -1.
    #println!("Modulus::", 21.0f64 % 4.0f64)          #  1
    #println!("Modulus::", modulus(-21.0f64, 4.0f64)) #  3 because -21 + 4 x 6 is 3.
    #println!("Modulus::", modulus(21.0f64, 4.0f64))  #  1
    return a - np.floor(a / b) * b


def interpolate_b_spline(tdata, ydata, tval):
    ## OPTIMIZATION: Return last left to reduce future searches
    # Based on:
    #
    #    SPLINE_B_VAL evaluates a cubic B spline approximant.
    #
    #  Discussion:
    #
    #    The cubic B spline will approximate the data, but is not
    #    designed to interpolate it.
    #
    #    In effect, two "phantom" data values are appended to the data,
    #    so that the spline will interpolate the first and last data values.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    14 August 2005
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Carl de Boor,
    #    A Practical Guide to Splines,
    #    Springer Verlag, 1978.

    yval  = 0.

    # Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    (left, right) = find_indices_around_target_value(tdata, tval)

    if left == right:
        # Target value out of range, use limit values
        yval = ydata[left]
    else:
        # Evaluate the 5 nonzero B spline basis functions in the interval,
        # weighted by their corresponding data values.
        u = (tval - tdata[left]) / (tdata[right] - tdata[left])

        # B function associated with node LEFT - 1, (or "phantom node"),
        # evaluated in its 4th interval.
        bval = ( 1.0 - 3.0 * u + 3.0 * np.power(u, 2) - np.power(u, 3) ) / 6.0
        if left > 0:
           yval = yval + ydata[left - 1] * bval
        else:
           yval = yval + ( 2.0 * ydata[0] - ydata[1] ) * bval


        # B function associated with node LEFT,
        # evaluated in its third interval.
        bval = ( 4.0 - 6.0 * np.power(u, 2) + 3.0 * np.power(u, 3) ) / 6.0
        yval = yval + ydata[left] * bval

        # B function associated with node RIGHT,
        # evaluated in its second interval.
        bval = ( 1.0 + 3.0 * u + 3.0 * np.power(u, 2) - 3.0 * np.power(u, 3) ) / 6.0
        yval = yval + ydata[right] * bval

        # B function associated with node RIGHT+1, (or "phantom node"),
        # evaluated in its first interval.
        bval = np.power(u, 3) / 6.0
        ndata = len(ydata)
        if right + 1 < ndata:
           yval = yval + ydata[right + 1] * bval
        elif ndata >= 2:
           yval = yval + ( 2.0 * ydata[ndata - 1] - ydata[ndata - 2] ) * bval

    return (yval, left)


def find_indices_around_target_value(data, target_value):
    # Find the nearest interval [ x(LEFT), x(RIGHT) ] to XVAL.
    data = np.asarray(data)
    right_idx = data.searchsorted(target_value)
    if right_idx == 0:
        return (0, 0)
    elif right_idx == len(data):
        return (right_idx-1, right_idx-1)
    else:
        return (right_idx-1, right_idx)



def calculate_pseudo_synchronization_period(semi_major_axis, eccentricity, star_mass, planet_mass):
    alpha = (1.+15./2.*np.power(eccentricity, 2)+45./8.*np.power(eccentricity, 4)+5./16.
                 * np.power(eccentricity, 6))*1./(1.+3.*np.power(eccentricity, 2)+3./8.
                 * np.power(eccentricity, 4))*1./np.power(1.-np.power(eccentricity, 2), 1.5)
    pseudo_rot = alpha * np.sqrt(G_SI*M_SUN*(star_mass+planet_mass)) # L^(3/2).T^(-1)
    angular_frequency = pseudo_rot * np.power(semi_major_axis*AU, -3./2.) * HOUR * 24. # days^-1
    pseudo_synchronization_period = TWO_PI/(angular_frequency) # days

    return pseudo_synchronization_period # days



def calculate_spin(angular_frequency, inclination, obliquity):
    x = 0.
    y = -1.*angular_frequency * np.sin(obliquity+inclination)
    z = angular_frequency * np.cos(obliquity+inclination)
    return Axes(x, y, z)
