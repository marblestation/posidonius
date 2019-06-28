import os
import pandas as pd
import numpy as np
import struct

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def read(filename):
    f = open(filename, "rb")
    # (np.floor(np.log10(np.max((100., 10.)))) - 2.)*10.

    if not os.path.exists(filename):
        raise Exception("File does not exists!")

    fields = ('current_time', 'time_step', 'particle', 'position_x', 'position_y', 'position_z', 'spin_x', 'spin_y', 'spin_z', 'velocity_x', 'velocity_y', 'velocity_z', 'semi-major_axis', 'perihelion_distance', 'eccentricity', 'inclination', 'longitude_of_perihelion', 'longitude_of_ascending_node', 'mean_anomaly', 'orbital_angular_momentum_x', 'orbital_angular_momentum_y', 'orbital_angular_momentum_z', 'orbital_angular_momentum', 'denergy_dt', 'migration_timescale', 'total_energy', 'total_angular_momentum', 'mass', 'radius', 'radius_of_gyration_2', 'scaled_dissipation_factor', 'love_number', 'lag_angle')

    data = []
    while True:
        try:
            row = f.read(8+8+4+8*(len(fields)-3))
            #vrow = struct.unpack('> d d i' + ' d'*(len(fields)-3), row)
            vrow = struct.unpack('< d d i' + ' d'*(len(fields)-3), row)
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

def classify(n_particles, data, discard_first_hundred_years=False):
    # Ignore first 100 years
    data['current_time'] /= 365.25 # From days to years
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

