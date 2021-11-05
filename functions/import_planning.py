from inputs.inputs import *

with open("./inputs/input_planning.json", "r") as read_file:
    input_planning = json.load(read_file)

# Assign the value of each Json's key to itself.
for key in input_planning:
    if isinstance(input_planning[key], list):
        globals()[key] = np.asarray(input_planning[key])
    else:
        globals()[key] = input_planning[key]