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


ANGLE_MIN = - 90.
ANGLE_MAX = 270.

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')
    parser.add_argument('inner_planet', action='store', type=int, help='Inner planet number')
    parser.add_argument('outer_planet', action='store', type=int, help='Outer planet number')
    parser.add_argument('inner_period', action='store', type=int, help='Inner period')
    parser.add_argument('outer_period', action='store', type=int, help='Outer period')

    args = parser.parse_args()

    if args.inner_planet > args.outer_planet:
        raise Exception("Inner planet number should be always smaller than the outer planet number")

    if args.outer_period > args.inner_period:
        raise Exception("Outer period cannot be greater than the inner period!")


    filename = args.historic_snapshot_filename
    n_particles, data = posidonius.analysis.history.read(filename)
    star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, discard_first_hundred_years=False)

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

    inner_planet = args.inner_planet
    outer_planet = args.outer_planet
    inner_period_nb = args.inner_period
    outer_period_nb = args.outer_period

    output_figure_dirname = os.path.dirname(filename)
    output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_resonance_{}_{}_between_{}_and_{}.png".format(inner_period_nb, outer_period_nb, inner_planet, outer_planet))

    t = t[inner_planet] # time in years (== t[outer_planet])
    a_inner = a[inner_planet - 1] # demi-grand axe en ua
    g_inner = g[inner_planet - 1] # argument of pericentre (degrees)
    n_inner = n[inner_planet - 1] # longitude of ascending node (degrees)
    M_inner = M[inner_planet - 1] # Mean anomaly (degrees)

    a_outer = a[outer_planet - 1] # demi-grand axe en ua
    g_outer = g[outer_planet - 1] # argument of pericentre (degrees)
    n_outer = n[outer_planet - 1] # longitude of ascending node (degrees)
    M_outer = M[outer_planet - 1] # Mean anomaly (degrees)

    sys.stdout.write("Calculating angles %5.1f %%                          \r" % (50.))
    sys.stdout.flush()

    nb_times = len(t)

    # Resonances are usually displayed as (p+q):p where q is the order of
    # the resonance. We retreive thoses parameters
    p = outer_period_nb
    q = inner_period_nb - outer_period_nb

    # We calculate the resonant angles
    long_of_peri_inner = g_inner + n_inner
    mean_longitude_inner = M_inner + long_of_peri_inner

    long_of_peri_outer = g_outer + n_outer
    mean_longitude_outer = M_outer + long_of_peri_outer

    phi = np.empty((q+1, nb_times)) # we create an array, empty for the moment, that will contain all the resonant angles associated with the supposed resonance.

    temp_value = inner_period_nb * mean_longitude_outer - outer_period_nb * mean_longitude_inner

    for i in range(q+1):
        phi[i] = temp_value - i * long_of_peri_inner - (q - i) * long_of_peri_outer

    # We take modulo 2*pi of the resonant angle
    phi = phi%(360.)
    too_low = phi < ANGLE_MIN
    too_high = phi > ANGLE_MAX
    phi[too_low] = phi[too_low] + 360.
    phi[too_high] = phi[too_high] - 360.

    delta_longitude = long_of_peri_outer - long_of_peri_inner
    delta_longitude = delta_longitude%(360.)
    too_low = delta_longitude < ANGLE_MIN
    too_high = delta_longitude > ANGLE_MAX
    delta_longitude[too_low] = delta_longitude[too_low] + 360.
    delta_longitude[too_high] = delta_longitude[too_high] - 360.


    # We chose the number of plot in x and y axis for the p.multi
    # environment in order to plot ALL the resonant angles and have x and
    # y numbers as close on from another as possible. There are q+1
    # resonant angles, the period ratio and w1 - w2, i.e q+3 plots
    (nb_lines, nb_rows) = posidonius.analysis.resonances.get_subplot_shape(q+3)

    # on trace les plots
    myxfmt = ScalarFormatter(useOffset=True)
    #myxfmt._set_offset(1e5)
    myxfmt.set_scientific(True)
    myxfmt.set_powerlimits((-3, 3))

    sys.stdout.write("Displaying %5.1f %%                          \r" % (75.))
    sys.stdout.flush()

    fig = plt.figure(figsize=(12,12))
    plt.clf()
    fig.subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.95, wspace=0.26, hspace=0.26)

    # We create a variable to store the index of the subplot
    subplot_index = 0
    fig.suptitle("resonance "+str(inner_period_nb)+":"+str(outer_period_nb)+" between "+str(inner_planet)+" and "+str(outer_planet))

    subplot_index += 1
    plot_period = fig.add_subplot(nb_lines, nb_rows, subplot_index)
    plot_period.plot(t, (a_outer / a_inner)**1.5)
    plot_period.set_xlabel("time [years]")
    plot_period.set_ylabel("period ratio")
    #~ plt.legend()
    plot_period.grid(True)

    # For each resonant angle, it is important to show only the points and not lines,
    # so that we do not mask interesting features by meaningless lines.
    subplot_index += 1
    plot_dl = fig.add_subplot(nb_lines, nb_rows, subplot_index, sharex=plot_period)
    plot_dl.plot(t, delta_longitude, '.')
    plot_dl.set_xlabel("time [years]")
    plot_dl.set_ylabel("w2 - w1")
    plot_dl.grid(True)

    for i in range(q+1):
      sys.stdout.write("Displaying %5.1f %%                          \r" % (75. + float((i+3) / (q+3)) * 25.))
      sys.stdout.flush()
      subplot_index += 1
      plot_phi = fig.add_subplot(nb_lines, nb_rows, subplot_index, sharex=plot_period)
      plot_phi.plot(t, phi[i], '.')
      plot_phi.set_xlabel("time [years]")
      plot_phi.set_ylabel(unicode("phi %i" % i, 'utf8'))
      plot_phi.grid(True)

    plot_period.xaxis.set_major_formatter(myxfmt)

    fig.savefig(output_figure_filename)
    print("Output figure file written to: {}".format(output_figure_filename))

    #plt.show() # We display the plot
