import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import struct
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')

    args = parser.parse_args()
    filename = args.historic_snapshot_filename


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

    planets_data = {}
    planets_keys = [] # To ensure the order
    for i in xrange(n_particles-1):
        planets_data["{}".format(i+1)] = data[data['particle'] == i+1]
        planets_keys.append("{}".format(i+1))

    ## Same number of points
    #shared_idx = np.min((len(star_data), len(planet_data)))
    #planet_data = planet_data[:shared_idx]
    #star_data = star_data[:shared_idx]


    ################################################################################
    ## Star
    ################################################################################
    star_norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
    star_rotation_period = 2*np.pi / star_norm_spin


    ################################################################################



    ################################################################################
    ## Planets
    ################################################################################
    total_planets_angular_momentum = 0.
    total_planets_orbital_angular_momentum = 0.
    planets_computed_data = {}
    for key in planets_keys:
        planet_data = planets_data[key]
        planet_computed_data = {}

        planet_norm_spin = np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
        planet_rotation_period = 2*np.pi / planet_norm_spin
        planet_computed_data['planet_rotation_period'] = planet_rotation_period

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
        planet_computed_data['planet_obliquity'] = planet_obliquity

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
        planet_computed_data['planet_precession_angle'] = planet_precession_angle

        planet_mass = planet_data['mass'][0]
        norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
        corrotation_radius = ((G*Msun*(star_mass+planet_mass))**(1/3.)) * ((norm_spin/day)**(-2./3.))/AU
        planet_computed_data['corrotation_radius'] = corrotation_radius

        e = planet_data['eccentricity']
        alpha = (1.+15./2.*e**2+45./8.*e**4+5./16.*e**6)*1./(1.+3.*e**2+3./8.*e**4)*1./(1.-e**2)**1.5
        pseudo_rot = alpha * np.sqrt(G*Msun*(star_mass+planet_mass))
        pseudo_synchronization_period  = 2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*AU)**(-3./2.) * hr)
        planet_computed_data['pseudo_synchronization_period'] = pseudo_synchronization_period


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
        planet_computed_data['gravitational_energy_lost'] = gravitational_energy_lost

        # The tidal heat flux depends on the eccentricity and on the obliquity of the planet.
        # If the planet has no obliquity, no eccentricity and if its rotation is synchronized, the tidal heat flux is zero.
        mean_tidal_flux = gravitational_energy_lost / (4 * np.pi * np.power(planet_data['radius'] * AU, 2))
        planet_computed_data['mean_tidal_flux'] = mean_tidal_flux


        denergy_dt = planet_data['denergy_dt'] * 6.90125e37 # conversation from Msun.AU^2.day^-3 to W
        planet_computed_data['denergy_dt'] = denergy_dt

        inst_tidal_flux = denergy_dt / (4 * np.pi * np.power(planet_data['radius'] * AU, 2))
        planet_computed_data['inst_tidal_flux'] = inst_tidal_flux

        ################################################################################
        ## Star obliquity (with respect to a planet)
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
        planet_computed_data['star_obliquity'] = star_obliquity
        ################################################################################


        ################################################################################
        ## Sum over all the planets values:
        planet_angular_momentum = planet_data['radius_of_gyration_2'] * (planet_data['mass']*Msun) * np.power(planet_data['radius']*AU, 2) * (planet_norm_spin/day)
        total_planets_angular_momentum += planet_angular_momentum # If more than one planet is present, all of them should be added

        # Sum on number of planets to have total orbital momentum
        planet_orbital_angular_momentum = (star_mass*Msun) * (planet_data['mass']*Msun) / ((star_mass*Msun) + (planet_data['mass']*Msun)) * planet_data['orbital_angular_momentum'] * ((AU*AU)/day) # kg.m^2.s-1
        total_planets_orbital_angular_momentum += planet_orbital_angular_momentum # If more than one planet is present, all of them should be added
        ################################################################################

        ################################################################################
        # Save computed data
        planets_computed_data[key] = planet_computed_data
        ################################################################################



    # \Delta L / L
    star_angular_momentum = star_data['radius_of_gyration_2'] * (star_mass*Msun) * np.power(star_data['radius']*AU, 2) * (star_norm_spin/day)
    initial_total_angular_momentum = total_planets_orbital_angular_momentum[0] + total_planets_angular_momentum[0] + star_angular_momentum[0]
    conservation_of_angular_momentum = np.abs(((total_planets_orbital_angular_momentum + total_planets_angular_momentum + star_angular_momentum) - initial_total_angular_momentum) / initial_total_angular_momentum)
    conservation_of_angular_momentum[0] = conservation_of_angular_momentum[1]
    #conservation_of_angular_momentum = np.abs(star_data['total_angular_momentum'] - star_data['total_angular_momentum'][0]) / star_data['total_angular_momentum'][0]



    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(5,3,1)
    field = 'semi-major_axis'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planet_data[field], label=key)
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['corrotation_radius'], label=key+" corrotation")
    ax.set_ylabel(field+" (AU)")
    #ax.set_ylim([0.005, 0.028])
    ax.set_xscale('log')
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(5,3,2, sharex=ax)
    field = 'planet obliquity (deg)'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['planet_obliquity'], label=key)
    ax.set_ylabel(field)
    ax.set_ylim([0.0001, 100.0])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(5,3,3, sharex=ax)
    field = 'eccentricity'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planet_data[field], label=key)
    ax.set_ylabel(field)
    ax.set_ylim([0.001, 1.000])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)


    ax = fig.add_subplot(5,3,4, sharex=ax)
    field = 'inclination'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planet_data[field] * (180 / np.pi), label=key) # From rad to degrees
    ax.set_ylabel('planet '+field+ " (deg)")
    #ax.set_ylim([0.5, 5.5])
    ax.set_xscale('log')
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(5,3,5, sharex=ax)
    # Energy loss dE/dt due to tides per planet surface
    field = 'Energy lost\ndue to tides (W/m^2)'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['inst_tidal_flux'], label=key+" Instantaneous") # Instantaneous energy loss
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['mean_tidal_flux'], label=key+" Mean") # Mean energy loss
    ax.set_ylabel(field)
    ax.set_ylim([1e-2, 1e5])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(5,3,6, sharex=ax)
    field = 'planet_rotation_period\n(hr)'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['planet_rotation_period']*24., label=key)
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['pseudo_synchronization_period'], label=key+" pseudo-sync")
    ax.set_ylabel(field)
    #ax.set_ylim([40, 160.0])
    ax.set_xscale('log')
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(5,3,7, sharex=ax)
    field = '$\Delta L/L_{0}$'
    ax.plot(planet_data['current_time'], conservation_of_angular_momentum)
    ax.set_ylabel(field)
    #ax.set_ylim([0., 0.000007])
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax = fig.add_subplot(5,3,8, sharex=ax)
    # Energy loss dE/dt due to tides
    field = 'Energy lost\ndue to tides (W)'
    #ax.plot(planet_data['current_time'], planet_data[field])
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['denergy_dt'], label=key+" Instantaneous") # Instantaneous energy loss
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['gravitational_energy_lost'], label=key+" Mean") # Mean energy loss
    ax.set_ylabel(field)
    ax.set_ylim([1e12, 1e19])
    ax.set_xscale('log')
    ax.set_yscale('symlog')
    ax.legend(loc=0, prop={'size':8})

    ax = fig.add_subplot(5,3,9, sharex=ax)
    field = 'star_rotation_period\n(days)'
    ax.plot(planet_data['current_time'], star_rotation_period)
    ax.set_ylabel(field)
    #ax.set_ylim([2.915, 2.92])
    ax.set_xscale('log')
    #plt.setp(ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(5,3,10, sharex=ax)
    field = 'star_obliquity (deg)'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['star_obliquity'], label=key)
    ax.set_ylabel(field)
    #ax.set_ylim([0.5, 5.5])
    ax.set_xscale('log')
    #ax.set_yscale('symlog')
    ax.legend(loc=0, prop={'size':8})

    ax = fig.add_subplot(5,3,11, sharex=ax)
    field = 'planet_precession_angle\n(deg)'
    for key in planets_keys:
        planet_data = planets_data[key]
        ax.plot(planet_data['current_time'], planets_computed_data[key]['planet_precession_angle'], label=key)
    ax.set_ylabel(field)
    #ax.set_ylim([80., 100.])
    ax.set_xscale('log')
    #ax.set_yscale('symlog')
    ax.legend(loc=0, prop={'size':8})

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

    ax = fig.add_subplot(5,3,12, sharex=ax)
    field = '$\Delta E/E_{0}$'
    ax.plot(planet_data['current_time'], relative_energy_error)
    ax.set_ylabel(field)
    ax.set_ylim([-0.35, 0.05])
    ax.set_xscale('log')
    #ax.set_yscale('symlog')

    ax = fig.add_subplot(5,3,13, sharex=ax)
    field = 'Star radius (AU)'
    ax.plot(star_data['current_time'], star_data['radius'])
    ax.set_ylabel(field)
    ax.set_xscale('log')

    ax = fig.add_subplot(5,3,14, sharex=ax)
    field = 'Star radius of\ngyration 2'
    ax.plot(star_data['current_time'], star_data['radius_of_gyration_2'])
    ax.set_ylabel(field)
    ax.set_xscale('log')

    ax = fig.add_subplot(5,3,15, sharex=ax)
    field = 'Star love number'
    ax.plot(star_data['current_time'], star_data['love_number'])
    ax.set_ylabel(field)
    ax.set_xscale('log')

    ax.set_xlim([100.0, 1.0e8])
    plt.tight_layout()

    output_figure_filename = os.path.dirname(filename) + "/" + os.path.splitext(os.path.basename(filename))[0] + ".png"
    #output_figure_filename = os.path.splitext(os.path.basename(filename))[0] + ".png"
    plt.savefig(output_figure_filename)
    #plt.show()
    print("Output figure file written to: {}".format(output_figure_filename))


    all_data = None
    for key in planets_keys:
        planet_data = planets_data[key]
        data = pd.DataFrame(planet_data['current_time'], columns=['current_time'])
        data['planet'] = key
        data['semi-major_axis_AU'] = planet_data['semi-major_axis']
        data['corrotation_radius_AU'] = planets_computed_data[key]['corrotation_radius']
        data['planet_obliquity_deg'] = planets_computed_data[key]['corrotation_radius']
        data['eccentricity'] = planet_data['eccentricity']
        data['inclination_deg'] = planet_data['inclination'] * (180 / np.pi)
        data['energy_lost_due_to_tides_W_per_m2'] = planets_computed_data[key]['inst_tidal_flux']
        data['mean_energy_lost_due_to_tides_W_per_m2'] = planets_computed_data[key]['mean_tidal_flux']
        data['planet_rotation_period_hours'] = planets_computed_data[key]['planet_rotation_period']*24
        data['planet_pseudo_synchronization_period'] = planets_computed_data[key]['pseudo_synchronization_period']
        data['energy_lost_due_to_tides_W'] = planets_computed_data[key]['denergy_dt']
        data['mean_energy_lost_due_to_tides_W'] = planets_computed_data[key]['gravitational_energy_lost']
        data['star_obliquity_deg'] = planets_computed_data[key]['star_obliquity']
        data['planet_precession_angle_deg'] = planets_computed_data[key]['planet_precession_angle']
        data['conservation_of_angular_momentum'] = conservation_of_angular_momentum
        data['star_rotation_period_days'] = star_rotation_period
        data['conservation_of_energy'] = relative_energy_error
        if all_data is None:
            all_data = data
        else:
            all_data = pd.concat((all_data, data))
    output_text_filename = os.path.dirname(filename) + "/" + os.path.splitext(os.path.basename(filename))[0] + ".txt"
    #output_text_filename = os.path.splitext(os.path.basename(filename))[0] + ".txt"
    all_data.to_csv(output_text_filename, sep="\t")

    print("Output plain text file written to: {}".format(output_text_filename))


