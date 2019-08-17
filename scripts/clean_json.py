import fnmatch
import os
import json
from itertools import chain

paths = ('tests/data', 'posidonius/tests/data',)
for root, dirnames, filenames in chain.from_iterable(os.walk(path) for path in paths):
    for filename in fnmatch.filter(filenames, 'case.json'):
        full_path = os.path.join(root, filename)
        with open(full_path) as input_json_file:
            data = json.load(input_json_file)
        with open(full_path, "w") as output_json_file:
            json.dump(data, output_json_file, indent=2, sort_keys=True)
        print("Transformed: {}".format(full_path))
