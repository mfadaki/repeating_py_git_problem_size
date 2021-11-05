#
# - Split the target horizon (e.g. 60) based on the computational power that the ccomputer can solve a multi-period problem.
# - For example, if computational power would be 7 days, then for the first iteration, we consider 7 days + 7 weeks.
#     - The average of week 2 is considered for day=8, average of week 3 for day=9, ..., day=13 (starting from 0).
#     - The Integer of average is considered.
# - if remainder_days=0, the number of period changes ass there is no last period with partial number of days.

import numpy as np


def planning(days, compt_power):
    # compt_power = 14 #in terms of days
    daily_period = int(compt_power / 2)
    weekly_period = int((days - daily_period) / daily_period)
    remainder_days = days - weekly_period * daily_period - daily_period
    print("daily_period = ", daily_period)
    print("weekly_period = ", weekly_period)
    print("remainder_days = ", remainder_days)

    phl = 1 + weekly_period + 1  # planning Horizon Length

    if remainder_days != 0:
        ph = np.zeros((phl, 2))  # planning_horizon
    else:
        ph = np.zeros((phl - 1, 2))  # planning_horizon

    start = 0
    for c in range(phl - 1):
        end = start + daily_period
        ph[c, 0] = start
        ph[c, 1] = end
        start = end

    if remainder_days != 0:
        ph[c + 1, 0] = start
        ph[c + 1, 1] = start + remainder_days

    ph = ph.astype(int)
    print("phl = ", phl)
    print("ph = \n", ph)

    return daily_period, weekly_period, remainder_days, phl, ph
