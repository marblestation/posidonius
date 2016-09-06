import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#from astropy.io import ascii

raw = np.loadtxt("../aeipnl.out")
ref = pd.DataFrame(raw, columns=("current_time", "semi-major_axis", "eccentricity", "inclination", "longitude_of_perihelion", "longitude_of_ascending_node", "mean_anomaly"))
ref = ref[['current_time']]
ref['ref'] = True
ref = ref.set_index('current_time')

raw = np.loadtxt("../spins.out")
star1 = pd.DataFrame(raw, columns=("current_time", "spin_x", "spin_y", "spin_z", "radius", "radius_of_gyration_2", "love_number", "dissipation_factor"))
star1 = star1.set_index('current_time')
star1['select'] = ref['ref']np


raw = np.loadtxt("../horb1.out")
star2 = pd.DataFrame(raw, columns=("current_time", "orbital_angular_momentum_x", "orbital_angular_momentum_y", "orbital_angular_momentum_z"))
del star2['current_time']

star_data = star1.join(star2)
star_data['mass'] = 0.08
star_data['position_x'] = 0.
star_data['position_y'] = 0.
star_data['position_z'] = 0.


raw = np.loadtxt("../spinp1.out")
planet1 = pd.DataFrame(raw, columns=("current_time", "spin_x", "spin_y", "spin_z", "radius", "radius_of_gyration_2"))

raw = np.loadtxt("../horb1.out")
planet2 = pd.DataFrame(raw, columns=("current_time", "orbital_angular_momentum_x", "orbital_angular_momentum_y", "orbital_angular_momentum_z"))
del planet2['current_time']

raw = np.loadtxt("../aeipnl.out")
planet3 = pd.DataFrame(raw, columns=("current_time", "semi-major_axis", "eccentricity", "inclination", "longitude_of_perihelion", "longitude_of_ascending_node", "mean_anomaly"))
del planet3['current_time']

planet4 = planet1.join(planet2)
planet_data = planet4.join(planet3)
planet_data['mass'] = 3.e-6
planet_data['position_x'] = 0.
planet_data['position_y'] = 0.
planet_data['position_z'] = 0.
planet_data['denergy_dt'] = 0.

#import struct
#filename = "output.bin"

#f = open(filename, "rb")
## (np.floor(np.log10(np.max((100., 10.)))) - 2.)*10.

#fields = ('current_time', 'time_step', 'particle', 'position_x', 'position_y', 'position_z', 'spin_x', 'spin_y', 'spin_z', 'velocity_x', 'velocity_y', 'velocity_z', 'acceleration_x', 'acceleration_y', 'acceleration_z', 'dspin_dt_x', 'dspin_dt_y', 'dspin_dt_z', 'torque_x', 'torque_y', 'torque_z', 'orthogonal_component_of_the_tidal_force_due_to_stellar_tide', 'orthogonal_component_of_the_tidal_force_due_to_planetary_tide', 'radial_component_of_the_tidal_force', 'radial_component_of_the_tidal_force_conservative_part', 'radial_component_of_the_tidal_force_dissipative_part', 'tidal_acceleration_x', 'tidal_acceleration_y', 'tidal_acceleration_z', 'radial_velocity', 'norm_velocity_vector', 'distance', 'semi-major_axis', 'perihelion_distance', 'eccentricity', 'inclination', 'longitude_of_perihelion', 'longitude_of_ascending_node', 'mean_anomaly', 'orbital_angular_momentum_x', 'orbital_angular_momentum_y', 'orbital_angular_momentum_z', 'orbital_angular_momentum', 'denergy_dt', 'mass', 'radius', 'radius_of_gyration_2', 'dissipation_factor', 'love_number')

#data = []
##i = 0
#while True:
    #try:
        #row = f.read(8+8+4+8*46)
        #vrow = struct.unpack('> d d i' + ' d'*46, row)
    #except:
        #break
    #else:
        #data.append(vrow)
        ##if data is None:
            ##data = pd.DataFrame(dict(zip(fields, vrow)), columns=fields, index=(i,))
        ##else:
            ##data = pd.concat((data, pd.DataFrame(dict(zip(fields, vrow)), columns=fields, index=(i,))))
        ##i += 1
#data = pd.DataFrame(data, columns=fields, index=np.arange(len(data)))
#if len(data) % 2:
    ## Force to have two lines per snapshot because there are 2 particles (sun+planet)
    #data = data[:-1]
#data = data.to_records()

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

#data = pd.read_csv("output.txt", sep="\t")
#data = ascii.read("../posidonius.ias15/output.txt", delimiter="\t")
#data = ascii.read("../target/output.txt", delimiter="\t")
#data = ascii.read("output.txt", delimiter="\t")
#data = data.as_array()
#data = data[data['current_time'] >= 100.]

#star = data['particle'] == 0
#star_data = data[star]
star_mass = star_data['mass'][0]

#planet = data['particle'] == 1
#planet_data = data[planet]


star_norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
planet_norm_spin = np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
planet_rotation_period = 2*np.pi / planet_norm_spin
star_rotation_period = 2*np.pi / star_norm_spin

## Planet obliquity
numerator = planet_data['orbital_angular_momentum_x'] * planet_data['spin_x'] + \
                    planet_data['orbital_angular_momentum_y'] * planet_data['spin_y'] + \
                    planet_data['orbital_angular_momentum_z'] * planet_data['spin_z']
denominator= np.sqrt(np.power(planet_data['orbital_angular_momentum_x'], 2) + np.power(planet_data['orbital_angular_momentum_y'], 2) + np.power(planet_data['orbital_angular_momentum_z'], 2)) * \
                        np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
planet_obliquity = numerator / denominator
ofilter = planet_obliquity <= 1.
planet_obliquity[ofilter] = np.arccos(planet_obliquity[ofilter])*180./np.pi
planet_obliquity[np.logical_not(ofilter)] = 1.e-6

## Star obliquity
numerator = planet_data['orbital_angular_momentum_x'] * star_data['spin_x'] + \
                    planet_data['orbital_angular_momentum_y'] * star_data['spin_y'] + \
                    planet_data['orbital_angular_momentum_z'] * star_data['spin_z']
denominator= np.sqrt(np.power(planet_data['orbital_angular_momentum_x'], 2) + np.power(planet_data['orbital_angular_momentum_y'], 2) + np.power(planet_data['orbital_angular_momentum_z'], 2)) * \
                        np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
star_obliquity = numerator / denominator
ofilter = star_obliquity <= 1.
star_obliquity[ofilter] = np.arccos(star_obliquity[ofilter])*180./np.pi
star_obliquity[np.logical_not(ofilter)] = 1.e-6

## Planet precession angle
numerator = 1./np.sin(planet_obliquity) \
            *(planet_data['position_x']*planet_data['spin_x'] \
                + planet_data['position_y']*planet_data['spin_y'] \
                + planet_data['position_z']*planet_data['spin_z'])

denominator = np.sqrt(np.power(star_data['position_x'], 2) + np.power(star_data['position_y'], 2) + np.power(star_data['position_z'], 2)) \
                * planet_norm_spin
planet_precession_angle = numerator / denominator
ofilter = planet_precession_angle <= 1.
planet_precession_angle[ofilter] = np.arccos(planet_precession_angle[ofilter])*180./np.pi
planet_precession_angle[np.logical_not(ofilter)] = 1.e-6




k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)

### Calculation of energydot and tidal flux, in W/m2
# Gravitationl energy lost of the system due to dissipation
# Masses in kg
gravitational_energy_lost = energydot(planet_data['semi-major_axis']*AU, \
                                        planet_data['eccentricity'], \
                                        planet_norm_spin / day, \
                                        planet_obliquity * np.pi/180.0, \
                                        G, \
                                        planet_data['mass'] * Msun, \
                                        star_mass * Msun, \
                                        planet_data['radius'] * Rsun, \
                                        k2pdelta * day)

# The tidal heat flux depends on the eccentricity and on the obliquity of the planet.
# If the planet has no obliquity, no eccentricity and if its rotation is synchronized, the tidal heat flux is zero.
tidal_flux = gravitational_energy_lost / (4 * np.pi * np.power(planet_data['radius'] * Rsun, 2))

denergy_dt = planet_data['denergy_dt'] * 6.90125e37 # conversation from Msun.AU^2.day^-3 to W
inst_tidal_flux = denergy_dt / (4 * np.pi * np.power(planet_data['radius'] * Rsun, 2))


star_angular_momentum = star_data['radius_of_gyration_2'] * (star_mass*Msun) * np.power(star_data['radius']*Rsun, 2) * (star_norm_spin/day)
planet_angular_momentum = planet_data['radius_of_gyration_2'] * (planet_data['mass']*Msun) * np.power(planet_data['radius']*Rsun, 2) * (planet_norm_spin/day)
total_planets_angular_momentum = planet_angular_momentum # If more than one planet is present, all of them should be added

# Sum on number of planets to have total orbital momentum
#planet_orbital_angular_momentum = (star_mass*Msun) * (planet_data['mass']*Msun) / ((star_mass*Msun) + (planet_data['mass']*Msun)) * planet_data['orbital_angular_momentum'] * ((AU*AU)/day) # kg.m^2.s-1
planet_orbital_angular_momentum = (star_mass*Msun) * (planet_data['mass']*Msun) / ((star_mass*Msun) + (planet_data['mass']*Msun)) \
                            * np.sqrt(np.power(planet_data['orbital_angular_momentum_x'], 2) + np.power(planet_data['orbital_angular_momentum_y'], 2) + np.power(planet_data['orbital_angular_momentum_z'], 2)) \
                            * ((AU*AU)/day) # kg.m^2.s-1
total_planets_orbital_angular_momentum = planet_orbital_angular_momentum # If more than one planet is present, all of them should be added



# \Delta L / L
initial_total_angular_momentum = total_planets_orbital_angular_momentum[0] + total_planets_angular_momentum[0] + star_angular_momentum[0]
conservation_of_angular_momentum = np.abs(((total_planets_orbital_angular_momentum + total_planets_angular_momentum + star_angular_momentum) - initial_total_angular_momentum) / initial_total_angular_momentum)
conservation_of_angular_momentum[0] = conservation_of_angular_momentum[1]

#conservation_of_angular_momentum = np.abs(total_planets_orbital_angular_momentum - total_planets_orbital_angular_momentum[0]) / total_planets_orbital_angular_momentum[0]
#conservation_of_angular_momentum = np.abs(total_planets_angular_momentum + star_angular_momentum - (total_planets_angular_momentum[0] + star_angular_momentum[0])) / (total_planets_angular_momentum[0] + star_angular_momentum[0])
#m = total_planets_angular_momentum + star_angular_momentum
#current = m[1:]
#previous = m[:-1]
#first = m[0]
#conservation_of_angular_momentum = (current-previous)/previous
#conservation_of_angular_momentum = np.hstack(([conservation_of_angular_momentum[0]], conservation_of_angular_momentum))

planet_mass = planet_data['mass'][0]
norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
corrotation_radius = ((G*Msun*(star_mass+planet_mass))**(1/3.)) * ((norm_spin/86400.)**(-2./3.))/AU

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
ax.set_ylim([0.001, 0.110])
ax.set_xscale('log')
ax.set_yscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)


ax = fig.add_subplot(4,3,4, sharex=ax)
field = 'inclination'
ax.plot(planet_data['current_time'], planet_data[field] * (180 / np.pi)) # From rad to degrees
ax.set_ylabel(field+ " (deg)")
ax.set_ylim([2.5, 5.5])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,5, sharex=ax)
field = 'tidal_flux (W/m^2)'
ax.plot(planet_data['current_time'], tidal_flux)
#ax.plot(planet_data['current_time'], inst_tidal_flux)
ax.set_ylabel(field)
#ax.set_ylim([0.001, 10000.0])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,6, sharex=ax)
field = 'planet_rotation_period (hr)'
ax.plot(planet_data['current_time'], planet_rotation_period*24.)
ax.set_ylabel(field)
#ax.set_ylim([40, 150.0])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,7, sharex=ax)
field = '$\Delta L/L_{0}$'
ax.plot(planet_data['current_time'], conservation_of_angular_momentum)
ax.set_ylabel(field)
#ax.set_ylim([40, 150.0])
ax.set_xscale('log')
#ax.set_yscale('symlog')

ax = fig.add_subplot(4,3,8, sharex=ax)
field = 'denergy_dt (W)'
#ax.plot(planet_data['current_time'], planet_data[field])
ax.plot(planet_data['current_time'], denergy_dt)
ax.set_ylabel(field)
#ax.set_ylim([2.5, 5.5])
ax.set_xscale('log')
ax.set_yscale('symlog')

ax = fig.add_subplot(4,3,9, sharex=ax)
field = 'star_rotation_period (days)'
ax.plot(planet_data['current_time'], star_rotation_period)
ax.set_ylabel(field)
#ax.set_ylim([40, 150.0])
ax.set_xscale('log')
#plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(4,3,10, sharex=ax)
field = 'star_obliquity (deg)'
ax.plot(planet_data['current_time'], star_obliquity)
ax.set_ylabel(field)
#ax.set_ylim([40, 150.0])
ax.set_xscale('log')
#ax.set_yscale('symlog')

ax = fig.add_subplot(4,3,11, sharex=ax)
field = 'planet_precession_angle (deg)'
ax.plot(planet_data['current_time'], planet_precession_angle)
ax.set_ylabel(field)
#ax.set_ylim([2.5, 5.5])
ax.set_xscale('log')
#ax.set_yscale('symlog')

ax.set_xlim([100.0, 1.0e7])

plt.tight_layout()
#plt.savefig("../target/output.png")
plt.savefig("output.png")
#plt.show()

import pudb
pudb.set_trace()


#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#field = 'eccentricity'
#ax.plot(np.log10(planet_data['current_time']), planet_data[field])
#ax.set_ylim([np.min(planet_data[field]), np.max(planet_data[field]) ])
#plt.show()



#tidefluxi(i,j) = enerdot(ai(i,j)*AU,ei(i,j),rotpi(i,j),oblpi(i,j)*!Pi/180.d0,G $
    #,mb(i,0)*Msun,Ms,Rpi(i,j)*Rearth,k2pDeltap(i)*day)/(4.d0*!Pi*(Rpi(i,j)*Rearth)^2)
#Ip(i) = rg2p(i,0)*mb(i,0)*Msun*(Rp(i,0)*rsun)^2
