# Instructions:
# Change the path
# If number of solutions in the first function is not right, set the number of json solutions manually in line after 40
# change the parameter name in the last function. for n-ij, we need to select nd_ijt

# import numpy as np
# import scipy as cp
# import pandas as pd
# import math
# import sys
# import os
# import json
# from json import JSONEncoder

from functions.import_planning import *
path = "./results/results_t6/solutions/"


# Selct the most recent JSON file in the 'results' directory

def no_sol_fn():
    import glob
    import os

    list_of_files = glob.glob(path + "*.json")
    latest_json_file = max(list_of_files, key=os.path.getctime)
    print(latest_json_file)

    # No of Jsons
    latest_json_file_filename_with_extention = os.path.basename(latest_json_file)
    latest_json_file_filename_without_extention = os.path.splitext(
        latest_json_file_filename_with_extention
    )[0]
    no_solions = int(latest_json_file_filename_without_extention)
    return no_solions


no_solions = no_sol_fn()
print(no_solions)

no_solions = 25  # Use it if we want to manually set the no of solutions

# Get total time ...


def t_tot_fn(no_solions):
    t_total = 0
    with open(path + "0.json", "r") as read_file:
        sol_json_0 = json.load(read_file)
    daily_period = sol_json_0["daily_period"]
    for _ in range(no_solions):
        with open(path + str(_) + ".json", "r") as read_file:
            sol_json = json.load(read_file)
        print(_)
        print("is_within_target_period= ", sol_json["is_within_target_period"])
        print("single_period= ", sol_json["single_period"])
        print("t= ", sol_json["t"])

        if sol_json["t"] >= daily_period:
            t = daily_period
        else:
            t = sol_json["t"]
        t_total = t_total + t
    return t_total


t_total = t_tot_fn(no_solions)
print("t_total= ", t_total)


def risk_fn(t_total):
    risk_multi = np.zeros((t_total)).astype(int)
    c_period = 0
    for c in range(no_solions + 1):
        with open(path + str(c) + ".json", "r") as read_file:
            sol_json = json.load(read_file)

        t = sol_json["t"]

        if t >= daily_period:
            risk_multi[c_period: c_period + daily_period] = np.asarray(
                sol_json["risk"]
            )[:daily_period]
            c_period += daily_period
        else:
            if t == 1:
                risk_multi[c_period: c_period + t] = np.asarray(sol_json["risk"])
            else:
                risk_multi[c_period: c_period + t] = np.asarray(sol_json["risk"])[:t]
            c_period += t
    return risk_multi


risk_multi = risk_fn(t_total)

np.savetxt("risk_multi.csv", risk_multi, delimiter=",")

# region Any Param Summation


def param_sum_fn(t_total, param):
    param_multi = np.zeros((t_total)).astype(int)
    c_period = 0
    for _ in range(no_solions):
        with open(path + str(_) + ".json", "r") as read_file:
            sol_json = json.load(read_file)

        t = sol_json["t"]
        print(t)
        # print(np.sum(np.asarray(sol_json["n_ij"])))
        # print(np.sum(np.asarray(sol_json['nd_ij_sv'])))
        # print(np.sum(np.asarray(sol_json["x_ij_sv"])[:, :, :1]))
        # print("***")

        if t >= daily_period:
            for __ in range(daily_period):
                param_multi[c_period] = np.sum(np.asarray(sol_json[param])[:, :, __])
                c_period += 1
        else:
            if t == 1:
                if param == "nd_ij_sv":
                    param = "n_ij"
                param_multi[c_period] = np.sum(np.asarray(sol_json[param]))
                c_period += 1
            else:
                for __ in range(t):
                    param_multi[c_period] = np.sum(
                        np.asarray(sol_json[param])[:, :, __]
                    )
                    c_period += 1
    return param_multi


param_multi = param_sum_fn(t_total, "nd_ij_sv")

np.savetxt("param_multi.csv", param_multi, delimiter=",")

# endregion

######################################################
