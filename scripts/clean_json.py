import fnmatch
import os
import json

for root, dirnames, filenames in os.walk('tests/data'):
    for filename in fnmatch.filter(filenames, 'case.json'):
        full_path = os.path.join(root, filename)
        with open(full_path) as input_json_file:
            data = json.load(input_json_file)
        with open(full_path, "w") as output_json_file:
            json.dump(data, output_json_file, indent=2, sort_keys=True)
        print("Transformed: {}".format(full_path))
