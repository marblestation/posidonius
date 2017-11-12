"""
Based on a script developed by Dr. Christophe Cossou
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
import struct
import argparse
import posidonius

################
## Parameters ##
################

## The running time of the script is very sensitive to this parameter that determine the number of time, for each resonance,
# that we search for a corresping fraction to the orbital period ratio
NUMBER_OF_VALUES = 10 # sampling for period ratio around the given value
DENOMINATOR_LIMIT = 12 # Maximum value allowed of the denominator when we want to get a fraction from a decimal value
NUMERATOR_LIMIT = 20 # maximum value allowed for the numerator
UNCERTAINTY = 5 # In percentage
NB_LAST_POINTS = 15 # Number of points we want to test the libration of angles.

NB_MEASUREMENTS = 500 # The number of times we test the resonances between planets (because the total number of output can vary from one simulation to another)

# the threshold of the standard deviation of a resonant angle,
# below which we consider there is libration and thus, a resonance.
STD_THRESHOLD = 70

# Extreme angles for angle folding. We will apply a mod 2pi to all angles,
# but theses values determines between wich values we will constrains the angles.
ANGLE_CENTER_VALUE = 90.

# We get arguments from the script
isLog = False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')

    args = parser.parse_args()
    filename = args.historic_snapshot_filename
    n_particles, data = posidonius.analysis.history.read(filename)
    star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, discard_first_hundred_years=False)

    output_figure_dirname = os.path.dirname(filename)
    if len(output_figure_dirname) > 0:
        output_figure_dirname += "/"
    output_figure_filename = output_figure_dirname + os.path.splitext(os.path.basename(filename))[0] + "_timed_resonances.png"

    planet_names = planets_keys
    nb_planets = len(planet_names)

    # t: time in years
    # a: demi-grand axe en ua
    # e: eccentricity
    # g: argument of pericentre or argument of perihelion (degrees)
    # n: longitude of ascending node (degrees)
    # M: Mean anomaly (degrees)
    # q: periastron or perihelion (AU)
    # Q: apastron or aphelion (AU)
    t, a, e, g, n, M, q, Q = posidonius.analysis.resonances.calculate_resonance_data(planets_keys, planets_data)
    max_lengths = len(t[0])

    # If the user require more measurement than timestep available, we force to have nb_measurements equal the number of timestep
    if (NB_MEASUREMENTS < max_lengths):
        time_delay = max_lengths / NB_MEASUREMENTS
    else:
        time_delay = 1

    # We want at least to test every single point once.
    # Thus, we make sure that "NB_LAST_POINTS" is at least equal to the space between two tests of the resonances
    if (time_delay > NB_LAST_POINTS):
        print("Warning: The interval between two checks of resonance is greater than the number of points used to test it.")
        print("         'NB_LAST_POINTS' extended to 'time_delay'=%d" % time_delay)
        NB_LAST_POINTS = time_delay

    ## We calculate q and Q
    #q = [ai * (1 - ei) for (ai, ei) in zip(a,e)]
    #Q = [ai * (1 + ei) for (ai, ei) in zip(a,e)]

    ####################
    # We declare the arrays needed for the resonances
    ####################
    # to store the successive resonances for each planet. If a resonance is written at the index i, that mean that the planet i is the inner planet of the resonance given in the array.
    resonance_index_range = [[] for planet in range(nb_planets)] # For each planet, each tuple is a range of time during which the planet is in resonance with the planet just after it.
    resonance_type = [[] for planet in range(nb_planets)] # For each planet, store the type of resonance 3:2, 4:3 and so on. The index i of the sublist (of a given planet) correspond to the range in time of the array resonance_time

    dynamic_order = [[] for planet in range(nb_planets)] # a list, for each planet, of its position through time. If 2, this means that the plant is the 2nd closest from its host star.
    time_order = [[] for planet in range(nb_planets)] # the corresponding time in years
    resonance_with = [[] for planet in range(nb_planets)] # the index of the outer planet for the resonance (the inner planet index being the index of the sublist considered)
    resonance_inner_rank = [[] for planet in range(nb_planets)] # The rank, in distance, of the inner planet in the resonance (useful to place the resonance in the second plot)


    ### We create the first point for some arrays
    # We sort the planet in the order of distance from the host star
    distance = [ai[0] for ai in a]
    planet_index_sorted_by_distance = np.argsort(distance) # the i-th closest planet in distance is the planet planet_index_sorted_by_distance[i]
    ordering_planets = 1 + np.argsort(planet_index_sorted_by_distance) # the i-th planet is the ordering_planets[i] in distance starting at 1

    # we append the ordering of the current planets in order to display resonances later
    for (planet_idx, order) in enumerate(ordering_planets):
        dynamic_order[planet_idx].append(order)
        time_order[planet_idx].append(t[planet_idx][0])
    ###
    for instant_index in range(NB_LAST_POINTS,max_lengths,time_delay):
        # We display a progress bar of the computation
        # The extra spaces are to make sure that no old character from the previous line will appear
        sys.stdout.write("Progression %6.2f %% : %i / %i                   \r" % ((instant_index * 100. / float(max_lengths)), instant_index, max_lengths))
        sys.stdout.flush()

        range_start = instant_index - (NB_LAST_POINTS-1)
        range_stop  = instant_index

        # We make sublist only with the planets that still exist at the present time
        still_here_planets = []
        a_temp = []
        g_temp = []
        n_temp = []
        M_temp = []

        for planet in range(nb_planets):
            if (t[planet].size > instant_index):
                still_here_planets.append(planet)
                a_temp.append(a[planet][range_start:range_stop+1])
                g_temp.append(g[planet][range_start:range_stop+1])
                n_temp.append(n[planet][range_start:range_stop+1])
                M_temp.append(M[planet][range_start:range_stop+1])

        # We sort the planet in the order of distance from the host star
        distance_end = [ai[-1] for ai in a_temp] # The distance of each planet from the host star, in AU, at the end of the current sub-range (of a_temp)
        distance_begin = [ai[0] for ai in a_temp] # The distance of each planet from the host star, in AU, at the beginning of the current sub-range (of a_temp)
        planet_index_sorted_by_distance = np.argsort(distance_end) # the i-th closest planet in distance is the planet planet_index_sorted_by_distance[i]
        ordering_planets = 1 + np.argsort(planet_index_sorted_by_distance) # the i-th planet is the ordering_planets[i] in distance (starting at 1)

        # we append the ordering of the current planets in order to display resonances later
        for (planet_idx, order) in zip(still_here_planets, ordering_planets):
            dynamic_order[planet_idx].append(order)
            time_order[planet_idx].append(t[planet_idx][instant_index])

        for (inner, outer) in zip(planet_index_sorted_by_distance[0:-1], planet_index_sorted_by_distance[1:]):
            g_inner = g_temp[inner]
            n_inner = n_temp[inner]
            M_inner = M_temp[inner]

            g_outer = g_temp[outer]
            n_outer = n_temp[outer]
            M_outer = M_temp[outer]

            # We get the various possible resonance between the two planets
            periodRatio_begin = (distance_begin[outer] / distance_begin[inner])**1.5
            periodRatio_end = (distance_end[outer] / distance_end[inner])**1.5

            # If the period ratio is too different between the beginning and the end of the range, we do not calculate possible resonances to gain time.
            if (abs(periodRatio_begin - periodRatio_end) > 0.02):
                continue

            resonances = posidonius.analysis.resonances.get_possible_resonances(periodRatio_end, uncertainty=0.01 * float(UNCERTAINTY),
                        denominator_limit=DENOMINATOR_LIMIT, numerator_limit=NUMERATOR_LIMIT, sampling=NUMBER_OF_VALUES)

            # For each resonance we check if this one exist between the two considered planets
            index = 0
            while (index < len(resonances)):
                res = resonances[index]
                is_resonance = posidonius.analysis.resonances.is_resonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer,
                              nb_points=NB_LAST_POINTS, angle_center_value=ANGLE_CENTER_VALUE, std_threshold=STD_THRESHOLD)
                if (is_resonance):
                    isExtend = False # boolean that say if the current resonance is the extension of the last resonance listed for the inner planet
                    if (len(resonance_type[still_here_planets[inner]]) != 0):
                        last_type = resonance_type[still_here_planets[inner]][-1]
                        last_index_range = resonance_index_range[still_here_planets[inner]][-1]

                        # if the two index ranges overlap
                        if ((last_type == res) and (last_index_range[1] >= range_start-1)):
                            isExtend = True

                    # if the current resonance already existed before, we only extend the index range of validity for the last resonance of the inner planet index
                    if isExtend:
                        resonance_index_range[still_here_planets[inner]][-1][1] = instant_index
                    else:
                        # We test here the previous resonance. If she was only at one instant, we delete it
                        if (len(resonance_index_range[still_here_planets[inner]]) != 0):
                            tmp = resonance_index_range[still_here_planets[inner]][-1]
                            if (tmp[0] == tmp[1]):
                                del(resonance_type[still_here_planets[inner]][-1])
                                del(resonance_index_range[still_here_planets[inner]][-1])
                                del(resonance_with[still_here_planets[inner]][-1])
                                del(resonance_inner_rank[still_here_planets[inner]][-1])

                        resonance_type[still_here_planets[inner]].append(res)
                        # To avoid overlap, we define the resonance at the position of the range_stop, without associating, by default, any length.
                        resonance_index_range[still_here_planets[inner]].append([range_stop, range_stop])
                        resonance_with[still_here_planets[inner]].append(still_here_planets[outer])
                        resonance_inner_rank[still_here_planets[inner]].append(ordering_planets[inner])

                    # If we find a resonance, we do not search for another one
                    break
                index += 1

    ####################
    # We now want to display in a fashion way the resonances
    ####################
    sys.stdout.write("Generating graphics                          \r")
    sys.stdout.flush()

    # We generate a list of colors
    #colors = [ '#'+li for li in autiwa.colorList(nb_planets)]
    import matplotlib.cm as cm
    colors = cm.rainbow(np.linspace(0, 1, nb_planets))

    fig = plt.figure(figsize=(7,14))
    plt.clf()
    fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)
    # We creat subplots. For subplot(311), it means that it has 2 lines & 3 columns and the current subplot is the first one (2*3 = 6 plots in total)
    plot_a = fig.add_subplot(311)
    for planet in range(nb_planets):
        sys.stdout.write("Generating graphics  %5.1f %%                          \r" % ((planet+1) * 25. / float(nb_planets)))
        sys.stdout.flush()
        if isLog:
            plot_a.semilogx(t[planet], a[planet], color=colors[planet], label=planet_names[planet])
            plot_a.semilogx(t[planet], q[planet], color=colors[planet])
            plot_a.semilogx(t[planet], Q[planet], color=colors[planet])
        else:
            plot_a.plot(t[planet], a[planet], color=colors[planet], label=planet_names[planet])
            plot_a.plot(t[planet], q[planet], color=colors[planet])
            plot_a.plot(t[planet], Q[planet], color=colors[planet])

    ylims = list(plt.ylim())
    for planet in range(nb_planets):
        sys.stdout.write("Generating graphics  %5.1f %%                          \r" % (20.+(planet+1) * 20. / float(nb_planets)))
        sys.stdout.flush()
        for (outer, (id_begin, id_end)) in zip(resonance_with[planet], resonance_index_range[planet]):
            plot_a.plot([t[planet][id_begin], t[planet][id_begin]], [a[planet][id_begin], a[outer][id_begin]], 'k--')
            plot_a.plot([t[planet][id_end], t[planet][id_end]], [a[planet][id_end], a[outer][id_end]], 'k--')

    # For a huge number of planet, the legend will be horrible, so we skip this line in that case
    if (nb_planets<=10):
        plot_a.legend()

    plot_a.set_xlabel("time [years]")
    plot_a.set_ylabel("a [AU]")
    plot_a.grid(True)

    plot_res = fig.add_subplot(312, sharex=plot_a)
    for planet in range(nb_planets):
        sys.stdout.write("Generating graphics  %5.1f %%                          \r" % (40.+(planet+1) * 20. / float(nb_planets)))
        sys.stdout.flush()
        for (res, inner_rank, (id_begin, id_end)) in zip(resonance_type[planet], resonance_inner_rank[planet], resonance_index_range[planet]):
            x_position = (t[planet][id_begin] + t[planet][id_end]) / 2.
            y_position = inner_rank + .5

            x_left = t[planet][id_begin]
            x_right = t[planet][id_end]
            y_bottom = inner_rank + 0.1
            y_top = inner_rank + 0.9

            text = "%i:%i" % (res.numerator, res.denominator)
            plot_res.text(x_position, y_position, text, horizontalalignment='center', verticalalignment='center', rotation='vertical')
            plot_res.fill([x_left, x_right, x_right, x_left],[y_bottom, y_bottom, y_top, y_top], fill=True, color='#cccccc')
            plot_res.plot([x_left, x_left], [y_bottom, y_top], 'k-')
            plot_res.plot([x_right, x_right], [y_bottom, y_top], 'k-')

    for planet in range(nb_planets):
        sys.stdout.write("Generating graphics  %5.1f %%                          \r" % (60.+(planet+1) * 20. / float(nb_planets)))
        sys.stdout.flush()
        if isLog:
            plot_res.semilogx(time_order[planet], dynamic_order[planet], color=colors[planet], label=planet_names[planet])
        else:
            plot_res.plot(time_order[planet], dynamic_order[planet], color=colors[planet], label=planet_names[planet])
    plot_res.set_xlabel("time [years]")
    plot_res.set_ylabel("planet order")
    plot_res.set_ylim((0, nb_planets+1))
    plot_res.yaxis.set_ticks(range(1, nb_planets+1))
    plot_res.yaxis.set_ticklabels(range(1, nb_planets+1))
    plot_res.grid(True)

    plot_e = fig.add_subplot(313, sharex=plot_a)
    for planet in range(nb_planets):
        sys.stdout.write("Generating graphics  %5.1f %%                          \r" % ((planet+1) * 25. / float(nb_planets)))
        sys.stdout.flush()
        if isLog:
            plot_e.loglog(t[planet], e[planet], color=colors[planet], label=planet_names[planet])
        else:
            plot_e.semilogy(t[planet], e[planet], color=colors[planet], label=planet_names[planet])
    plot_e.set_xlabel("time [years]")
    plot_e.set_ylabel("eccentricity")
    plot_e.grid(True)

    myxfmt = ScalarFormatter(useOffset=True)
    #myxfmt._set_offset(1e5)
    myxfmt.set_scientific(True)
    myxfmt.set_powerlimits((-3, 3))
    plot_a.xaxis.set_major_formatter(myxfmt)
    plot_res.yaxis.set_major_formatter(FormatStrFormatter('%i'))

    sys.stdout.write("Saving graphics                          \r")
    sys.stdout.flush()
    plt.savefig(output_figure_filename)
    print("Output figure file written to: {}".format(output_figure_filename))

    #sys.stdout.write("Displaying graphics                          \n")
    #sys.stdout.flush()
    #plt.show()

    ## TODO
    # Store the list of fraction for a given period ratio. If a period ratio is less than
    #  delta ratio/2 of an existing one, we do not calculate again the fraction, but instead retrieve the existing list

    ## Tricks
    # One thing to understand is the fact that when checking resonances at t=ti, we actually search for a resonance between t=ti-dt and t=ti.
    # If, for a particular reason, the current resonance can be extended, the resonance validity is extended to ti-dt to ti'=ti+step
    # where step is the separation in time between each check of resonances.
    # Though, if step is smaller than dt, then a new resonance can be declare in a range that overlap the previous one.
