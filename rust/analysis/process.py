import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii

import struct
filename = "../target/universe_history.bin"

f = open(filename, "rb")
# (np.floor(np.log10(np.max((100., 10.)))) - 2.)*10.

fields = ('current_time', 'time_step', 'particle', 'position_x', 'position_y', 'position_z', 'spin_x', 'spin_y', 'spin_z', 'velocity_x', 'velocity_y', 'velocity_z', 'semi-major_axis', 'perihelion_distance', 'eccentricity', 'inclination', 'longitude_of_perihelion', 'longitude_of_ascending_node', 'mean_anomaly', 'orbital_angular_momentum_x', 'orbital_angular_momentum_y', 'orbital_angular_momentum_z', 'orbital_angular_momentum', 'denergy_dt', 'total_energy', 'total_angular_momentum', 'mass', 'radius', 'radius_of_gyration_2', 'scaled_dissipation_factor', 'love_number')

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
if len(data) % 2:
    # Force to have two lines per snapshot because there are 2 particles (sun+planet)
    data = data[:-1]
data = data.to_records()

#-------------------------------------------------------------------------------
# Constants
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

#-------------------------------------------------------------------------------
# Functions
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
# Main
#-------------------------------------------------------------------------------
# Ignore first 100 years
data['current_time'] /= 362.25 # From days to years
data = data[data['current_time'] >= 100.]

star = data['particle'] == 0
star_data = data[star]
star_mass = star_data['mass'][0]

planet = data['particle'] == 1
#planet = data['particle'] == 2
planet_data = data[planet]

# Same number of points
shared_idx = np.min((len(star_data), len(planet_data)))
planet_data = planet_data[:shared_idx]
star_data = star_data[:shared_idx]


star_norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
planet_norm_spin = np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
planet_rotation_period = 2*np.pi / planet_norm_spin
star_rotation_period = 2*np.pi / star_norm_spin

## Planet obliquity
relative_orbital_angular_momentum_x = planet_data['orbital_angular_momentum_x']/planet_data['orbital_angular_momentum']
relative_orbital_angular_momentum_y = planet_data['orbital_angular_momentum_y']/planet_data['orbital_angular_momentum']
relative_orbital_angular_momentum_z = planet_data['orbital_angular_momentum_z']/planet_data['orbital_angular_momentum']
numerator = relative_orbital_angular_momentum_x * planet_data['spin_x'] + \
                    relative_orbital_angular_momentum_y * planet_data['spin_y'] + \
                    relative_orbital_angular_momentum_z * planet_data['spin_z']
denominator= np.sqrt(np.power(relative_orbital_angular_momentum_x, 2) + \
                        np.power(relative_orbital_angular_momentum_y, 2) + \
                        np.power(relative_orbital_angular_momentum_z, 2)) * \
                        np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
planet_obliquity = numerator / denominator
ofilter = planet_obliquity <= 1.
planet_obliquity[ofilter] = np.arccos(planet_obliquity[ofilter])*180./np.pi
planet_obliquity[np.logical_not(ofilter)] = 1.e-6

## Star obliquity
# https://en.wikipedia.org/wiki/Axial_tilt
# Angle between the spin axis of the star and the planet's orbital plane (in other words, it's the inclination of the planet)
# or, equivalently, the angle between its equatorial plane and orbital plane
# At an obliquity of zero, the two axes point in the same direction; i.e., the rotational axis is perpendicular to the orbital plane.
numerator = relative_orbital_angular_momentum_x * star_data['spin_x'] + \
                    relative_orbital_angular_momentum_y * star_data['spin_y'] + \
                    relative_orbital_angular_momentum_z * star_data['spin_z']
denominator= np.sqrt(np.power(relative_orbital_angular_momentum_x, 2) + \
                        np.power(relative_orbital_angular_momentum_y, 2) + \
                        np.power(relative_orbital_angular_momentum_z, 2)) * \
                        np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
star_obliquity = numerator / denominator
ofilter = star_obliquity <= 1.
star_obliquity[ofilter] = np.arccos(star_obliquity[ofilter])*180./np.pi
star_obliquity[np.logical_not(ofilter)] = 1.e-6

## Planet precession angle
# https://en.wikipedia.org/wiki/Apsidal_precession
# Angle defined between an arbitrary direction and the semi-major axis direction
numerator = 1./np.sin(planet_obliquity) \
            *(planet_data['position_x']*planet_data['spin_x'] \
                + planet_data['position_y']*planet_data['spin_y'] \
                + planet_data['position_z']*planet_data['spin_z'])

denominator = np.sqrt(np.power(planet_data['position_x'], 2) + np.power(planet_data['position_y'], 2) + np.power(planet_data['position_z'], 2)) \
                * planet_norm_spin
planet_precession_angle = numerator / denominator
ofilter = planet_precession_angle <= 1.
planet_precession_angle[ofilter] = np.arccos(planet_precession_angle[ofilter])*180./np.pi
planet_precession_angle[np.logical_not(ofilter)] = 1.e-6


### Calculation of energydot and tidal flux, in W/m2
# Gravitationl energy lost of the system due to dissipation
# Masses in kg
k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
gravitational_energy_lost = energydot(planet_data['semi-major_axis']*AU, \
                                        planet_data['eccentricity'], \
                                        planet_norm_spin / day, \
                                        planet_obliquity * np.pi/180.0, \
                                        G, \
                                        planet_data['mass'] * Msun, \
                                        star_mass * Msun, \
                                        planet_data['radius'] * AU, \
                                        k2pdelta * day)

# The tidal heat flux depends on the eccentricity and on the obliquity of the planet.
# If the planet has no obliquity, no eccentricity and if its rotation is synchronized, the tidal heat flux is zero.
mean_tidal_flux = gravitational_energy_lost / (4 * np.pi * np.power(planet_data['radius'] * AU, 2))


denergy_dt = planet_data['denergy_dt'] * 6.90125e37 # conversation from Msun.AU^2.day^-3 to W
inst_tidal_flux = denergy_dt / (4 * np.pi * np.power(planet_data['radius'] * AU, 2))


star_angular_momentum = star_data['radius_of_gyration_2'] * (star_mass*Msun) * np.power(star_data['radius']*AU, 2) * (star_norm_spin/day)
planet_angular_momentum = planet_data['radius_of_gyration_2'] * (planet_data['mass']*Msun) * np.power(planet_data['radius']*AU, 2) * (planet_norm_spin/day)
total_planets_angular_momentum = planet_angular_momentum # If more than one planet is present, all of them should be added

# Sum on number of planets to have total orbital momentum
planet_orbital_angular_momentum = (star_mass*Msun) * (planet_data['mass']*Msun) / ((star_mass*Msun) + (planet_data['mass']*Msun)) * planet_data['orbital_angular_momentum'] * ((AU*AU)/day) # kg.m^2.s-1
total_planets_orbital_angular_momentum = planet_orbital_angular_momentum # If more than one planet is present, all of them should be added



# \Delta L / L
initial_total_angular_momentum = total_planets_orbital_angular_momentum[0] + total_planets_angular_momentum[0] + star_angular_momentum[0]
conservation_of_angular_momentum = np.abs(((total_planets_orbital_angular_momentum + total_planets_angular_momentum + star_angular_momentum) - initial_total_angular_momentum) / initial_total_angular_momentum)
conservation_of_angular_momentum[0] = conservation_of_angular_momentum[1]
#conservation_of_angular_momentum = np.abs(star_data['total_angular_momentum'] - star_data['total_angular_momentum'][0]) / star_data['total_angular_momentum'][0]


planet_mass = planet_data['mass'][0]
norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
corrotation_radius = ((G*Msun*(star_mass+planet_mass))**(1/3.)) * ((norm_spin/day)**(-2./3.))/AU

e = planet_data['eccentricity']
alpha = (1.+15./2.*e**2+45./8.*e**4+5./16.*e**6)*1./(1.+3.*e**2+3./8.*e**4)*1./(1.-e**2)**1.5
pseudo_rot = alpha * np.sqrt(G*Msun*(star_mass+planet_mass))
pseudo_synchronization_period  = 2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*AU)**(-3./2.) * hr)

fig = plt.figure(figsize=(16, 10))
ax = fig.add_subplot(4,3,1)
field = 'semi-major_axis'
ax.plot(planet_data['current_time'], planet_data[field])
ax.plot(planet_data['current_time'], corrotation_radius, color="red")
ax.set_ylabel(field+" (AU)")
ax.set_ylim([0.005, 0.028])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,2, sharex=ax)
field = 'planet obliquity (deg)'
ax.plot(planet_data['current_time'], planet_obliquity)
ax.set_ylabel(field)
ax.set_ylim([0.0001, 100.0])
ax.set_xscale('log')
ax.set_yscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,3, sharex=ax)
field = 'eccentricity'
ax.plot(planet_data['current_time'], planet_data[field])
ax.set_ylabel(field)
ax.set_ylim([0.001, 1.000])
ax.set_xscale('log')
ax.set_yscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)


ax = fig.add_subplot(4,3,4, sharex=ax)
field = 'inclination'
ax.plot(planet_data['current_time'], planet_data[field] * (180 / np.pi)) # From rad to degrees
ax.set_ylabel('planet '+field+ " (deg)")
ax.set_ylim([2.5, 5.5])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,5, sharex=ax)
# Energy loss dE/dt due to tides per planet surface
field = 'Energy lost\ndue to tides (W/m^2)'
ax.plot(planet_data['current_time'], inst_tidal_flux) # Instantaneous energy loss
ax.plot(planet_data['current_time'], mean_tidal_flux, color="red") # Mean energy loss
ax.set_ylabel(field)
ax.set_ylim([1e-2, 1e5])
ax.set_xscale('log')
ax.set_yscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,6, sharex=ax)
field = 'planet_rotation_period\n(hr)'
ax.plot(planet_data['current_time'], planet_rotation_period*24.)
ax.plot(planet_data['current_time'], pseudo_synchronization_period, color="red")
ax.set_ylabel(field)
ax.set_ylim([40, 160.0])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,7, sharex=ax)
field = '$\Delta L/L_{0}$'
ax.plot(planet_data['current_time'], conservation_of_angular_momentum)
ax.set_ylabel(field)
ax.set_ylim([0., 0.000007])
ax.set_xscale('log')
#ax.set_yscale('symlog')

ax = fig.add_subplot(4,3,8, sharex=ax)
# Energy loss dE/dt due to tides
field = 'Energy lost\ndue to tides (W)'
#ax.plot(planet_data['current_time'], planet_data[field])
ax.plot(planet_data['current_time'], denergy_dt) # Instantaneous energy loss
ax.plot(planet_data['current_time'], gravitational_energy_lost, color="red") # Mean energy loss
ax.set_ylabel(field)
ax.set_ylim([1e12, 1e19])
ax.set_xscale('log')
ax.set_yscale('symlog')

ax = fig.add_subplot(4,3,9, sharex=ax)
field = 'star_rotation_period\n(days)'
ax.plot(planet_data['current_time'], star_rotation_period)
ax.set_ylabel(field)
ax.set_ylim([2.915, 2.92])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,10, sharex=ax)
field = 'star_obliquity (deg)'
ax.plot(planet_data['current_time'], star_obliquity)
ax.set_ylabel(field)
ax.set_ylim([2.5, 5.5])
ax.set_xscale('log')
#ax.set_yscale('symlog')

ax = fig.add_subplot(4,3,11, sharex=ax)
field = 'planet_precession_angle\n(deg)'
ax.plot(planet_data['current_time'], planet_precession_angle)
ax.set_ylabel(field)
ax.set_ylim([80., 100.])
ax.set_xscale('log')
#ax.set_yscale('symlog')

# How to compute the total energy if it is not provided:
#star_e_kin = 0.5 * star_data['mass'] * (np.power(star_data['velocity_x'], 2) + \
                                        #np.power(star_data['velocity_y'], 2) + \
                                        #np.power(star_data['velocity_z'], 2))
#planet_e_kin = 0.5 * planet_data['mass'] * (np.power(planet_data['velocity_x'], 2) + \
                                        #np.power(planet_data['velocity_y'], 2) + \
                                        #np.power(planet_data['velocity_z'], 2))
#e_kin = star_e_kin + planet_e_kin
#dx = planet_data['position_x'] - star_data['position_x']
#dy = planet_data['position_y'] - star_data['position_y']
#dz = planet_data['position_z'] - star_data['position_z']
#K2 = 0.01720209895**2
#e_pot = (-1. * K2 * (planet_data['mass']) * (star_data['mass']))  / np.sqrt(np.power(dx, 2) + np.power(dy, 2) + np.power(dz, 2))
#total_energy = e_kin + e_pot
## conservation of energy (kinetic+potential)
#relative_energy_error = (total_energy - total_energy[0]) / total_energy[0]

# conservation of energy (kinetic+potential)
relative_energy_error = (star_data['total_energy'] - star_data['total_energy'][0]) / star_data['total_energy'][0]

ax = fig.add_subplot(4,3,12, sharex=ax)
field = '$\Delta E/E_{0}$'
ax.plot(planet_data['current_time'], relative_energy_error)
ax.set_ylabel(field)
ax.set_ylim([-0.35, 0.05])
ax.set_xscale('log')
#ax.set_yscale('symlog')

ax.set_xlim([100.0, 1.0e8])

plt.tight_layout()
#plt.savefig("../target/output.png")
#plt.savefig("output.png")
#plt.savefig(os.path.dirname(filename) + "/" + os.path.splitext(os.path.basename(filename))[0] + ".png")
plt.savefig(os.path.splitext(os.path.basename(filename))[0] + ".png")
#plt.show()


data = pd.DataFrame(planet_data['current_time'], columns=['current_time'])
data['semi-major_axis_AU'] = planet_data['semi-major_axis']
data['corrotation_radius_AU'] = corrotation_radius
data['planet_obliquity_deg'] = planet_obliquity
data['eccentricity'] = planet_data['eccentricity']
data['inclination_deg'] = planet_data['inclination'] * (180 / np.pi)
data['energy_lost_due_to_tides_W_per_m2'] = inst_tidal_flux
data['mean_energy_lost_due_to_tides_W_per_m2'] = mean_tidal_flux
data['planet_rotation_period_hours'] = planet_rotation_period*24
data['planet_pseudo_synchronization_period'] = pseudo_synchronization_period
data['conservation_of_angular_momentum'] = conservation_of_angular_momentum
data['energy_lost_due_to_tides_W'] = denergy_dt
data['mean_energy_lost_due_to_tides_W'] = gravitational_energy_lost
data['star_rotation_period_days'] = star_rotation_period
data['star_obliquity_deg'] = star_obliquity
data['planet_precession_angle_deg'] = planet_precession_angle
data['conservation_of_energy'] = relative_energy_error
data.to_csv(os.path.splitext(os.path.basename(filename))[0] + ".txt", sep="\t")

import pudb
pudb.set_trace()

