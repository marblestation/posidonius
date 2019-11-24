import os
import pandas as pd
import numpy as np
from numpy.lib.recfunctions import append_fields
import struct
import posidonius.tools
import posidonius.constants
from posidonius.particles.axes import Axes

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def read(filename):
    f = open(filename, "rb")
    # (np.floor(np.log10(np.max((100., 10.)))) - 2.)*10.

    if not os.path.exists(filename):
        raise Exception("File does not exists!")

    fields = ('current_time', 'time_step', 'particle', 'position_x', 'position_y', 'position_z', 'spin_x', 'spin_y', 'spin_z', 'velocity_x', 'velocity_y', 'velocity_z', 'mass', 'radius', 'radius_of_gyration_2', 'love_number', 'scaled_dissipation_factor', 'lag_angle', 'denergy_dt', 'migration_timescale', )

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

def classify(n_particles, data, reference_particle_index=0, discard_first_hundred_years=False):
    # Ignore first 100 years
    data['current_time'] /= 365.25 # From days to years
    if discard_first_hundred_years:
        data = data[data['current_time'] >= 100.]

    star = data['particle'] == reference_particle_index
    star_data = data[star]

    planets_data = {}
    planets_keys = [] # To ensure the order
    for i in range(n_particles-1):
        planets_data["{}".format(i+1)] = data[data['particle'] == i+1]
        planets_keys.append("{}".format(i+1))

    # From inertial to heliocentric coordinates
    for key in planets_keys:
        planets_data[key]['position_x'] -= star_data['position_x']
        planets_data[key]['position_y'] -= star_data['position_y']
        planets_data[key]['position_z'] -= star_data['position_z']
        planets_data[key]['velocity_x'] -= star_data['velocity_x']
        planets_data[key]['velocity_y'] -= star_data['velocity_y']
        planets_data[key]['velocity_z'] -= star_data['velocity_z']
        semimajor_axis = []
        eccentricity = []
        inclination = []
        gm = posidonius.constants.G * (star_data['mass'] + planets_data[key]['mass'])
        position = Axes(planets_data[key]['position_x'], planets_data[key]['position_y'], planets_data[key]['position_z'])
        velocity = Axes(planets_data[key]['velocity_x'], planets_data[key]['velocity_y'], planets_data[key]['velocity_z'])
        a, q, e, i, p, n, l = posidonius.tools.calculate_keplerian_orbital_elements(gm, position, velocity)
        planets_data[key] = append_fields(planets_data[key], ('semi-major_axis', 'eccentricity', 'inclination'), (a, e, i), usemask=False)
    star_data['position_x'] = 0.
    star_data['position_y'] = 0.
    star_data['position_z'] = 0.
    star_data['velocity_x'] = 0.
    star_data['velocity_y'] = 0.
    star_data['velocity_z'] = 0.
    return star_data, planets_data, planets_keys

