import os
import numpy as np
import pandas as pd
import argparse
import posidonius

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')

    args = parser.parse_args()

    filename = args.historic_snapshot_filename
    n_particles, data = posidonius.analysis.history.read(filename)

    output_text_dirname = os.path.dirname(filename)
    for particle in np.unique(data['particle']):
        subdata = data[data['particle'] == particle].copy()
        subdata = pd.DataFrame(subdata)
        del subdata['index']
        del subdata['particle']
        output_text_filename = os.path.join(output_text_dirname, os.path.splitext(os.path.basename(filename))[0] + "_{}.txt".format(particle))
        subdata.to_csv(output_text_filename, sep="\t", index=False)
        print("> Raw output for particle '{}' written to plain text file: {}".format(particle, output_text_filename))


