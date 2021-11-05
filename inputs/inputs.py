# region Libraries
import numpy as np
import scipy as cp
import pandas as pd
import pickle
from datetime import datetime
import math
import sys
import docplex.mp
from docplex.mp.model import Model

# endregion

# region json
import json
from json import JSONEncoder


class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyArrayEncoder, self).default(obj)


# endregion

# region Import Inputs
with open("./inputs/inputs.json", "r") as read_file:
    inputs = json.load(read_file)

# Assign the value of each Json's key to itself.
for key in inputs:
    if isinstance(inputs[key], list):
        globals()[key] = np.asarray(inputs[key])
    else:
        globals()[key] = inputs[key]


# endregion
