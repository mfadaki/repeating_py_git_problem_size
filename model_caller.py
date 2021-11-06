# region Inputs and Libraries
from functions.aug_nij import *
from functions.model_singleP import *
from functions.model import *
import os
import glob
from functions.create_nij import *
from inputs.inputs import *

# endregion

# region Planning Horizon

from functions.planning import *

daily_period, weekly_period, remainder_days, phl, ph = planning(60, 14)
# the second arg is comp time like 14

input_planning = {
    "daily_period": daily_period,
    "weekly_period": weekly_period,
    "remainder_days": remainder_days,
    "phl": phl,
    "ph": ph,
}

with open("./inputs/input_planning.json", "w") as write_file:
    json.dump(input_planning, write_file, cls=NumpyArrayEncoder)

# endregion

# region Create n_ij
# Create the n_ij for each sub_multi-period and save it in a json file.
#
# - If weekly_period >= daily_period, from the start, we have a full daily period and at least one weekly period. For example, daily_period=5, days=100 --> weekly period = 19. So, we consider the multiperiods of 5*2=10 for a number of times until the period that there will be less than 5 weeks left to the end of horizon. So, for these terminating stages, the length of the multi-period reduces from 10, to 9, then 8, then 7 ....
# - On the other hand, if weekly_period < daily_period, from the start, we need to start with n_ij with dimension of less than t=2*daily_period.
# - We need to continue considering the 5d + 5w until, period 95. after that, size of n-ij is reduced by 1 for each iteration as well as the number of weeks to be averaged. 'cc' is avtivated when we get to the critical point of period 95 in this example.


create_nij(daily_period, weekly_period, remainder_days, phl, ph)
# endregion

# region Most recent Json
# Selct the most recent JSON file in the 'results' directory

list_of_files = glob.glob("./results/n_ij/*.json")
latest_json_file = max(list_of_files, key=os.path.getctime)
print(latest_json_file)

print(list_of_files)

print(list_of_files[0])

latest_json_file_filename_with_extention = os.path.basename(latest_json_file)
latest_json_file_filename_without_extention = os.path.splitext(
    latest_json_file_filename_with_extention
)[0]
no_n_ij = int(latest_json_file_filename_without_extention)
no_n_ij
# endregion

# region Loop

for c_main in range(no_n_ij + 1):
    # for c_main in range(3): #c_main is number of cycles
    with open("./results/n_ij/" + str(c_main) + ".json", "r") as read_file:
        n_ij_doc = json.load(read_file)
    n_ij = np.asarray(n_ij_doc["n_ij"])
    t = len(n_ij[0][0])
    print(t)

    if c_main != 0:
        for ci in range(i):
            for cj in range(j):
                n_ij[ci, cj, 0] = n_ij[ci, cj, 0] + last_unmet_demand[ci, cj]

    # exec(open("model.py").read())
    # sol_to_json(c_main, 1)
    model_fn(1, n_ij, c_main, t)

    with open("./results/solutions/" + str(c_main) + ".json", "r") as read_file:
        result_last_day = json.load(read_file)
    last_unmet_demand = np.asarray(result_last_day["unmet_demand"])


print(last_unmet_demand)

# endregion

# region (After Target Horizon) Single Period
# ### Periods after the target horizon (e.g. 60)

# gap_stop = 0.0001
# c_main = 8

with open("./results/solutions/" + str(c_main) + ".json", "r") as read_file:
    result_last_day = json.load(read_file)
last_unmet_demand = np.asarray(result_last_day["unmet_demand"])


compt_power = daily_period * 2
while np.sum(last_unmet_demand) > 0.1:
    c_main += 1
    sum_ud = np.sum(last_unmet_demand)
    avl_stock = k * s * compt_power

    if sum_ud > avl_stock:
        t_next = compt_power
    else:
        t_next_temp = np.sum(last_unmet_demand) / (k * s)
        if t_next_temp < 1 and t_next_temp > 0:
            t_next = 1
        else:
            t_next = int(t_next_temp) + 1

    n_ij = np.zeros((i, j, t_next))
    for ct in range(t_next):
        for ci in range(i):
            for cj in range(j):
                n_ij[ci, cj, ct] = int(last_unmet_demand[ci, cj] / t_next)

    nij_json = {"n_ij": n_ij}
    with open("./results/n_ij/" + str(c_main) + ".json", "w") as write_file:
        json.dump(nij_json, write_file, cls=NumpyArrayEncoder)

    t = len(n_ij[0][0])

    if t != 1:
        # exec(open("model.py").read())
        # sol_to_json(c_main, 0)
        model_fn(0, n_ij, c_main, t)
    else:
        n_ij = n_ij[:, :, 0]
        # exec(open("model_singleP.py").read())
        n_ij = aug_nij_fn(n_ij)
        model_singleP_fn(n_ij, c_main, t)

    with open("./results/solutions/" + str(c_main) + ".json", "r") as read_file:
        result_last_day = json.load(read_file)
    last_unmet_demand = np.asarray(result_last_day["unmet_demand"])

# endregion
