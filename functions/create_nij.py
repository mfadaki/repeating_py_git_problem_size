from inputs.inputs import *

n_ij_original = n_ij.copy()


def save_json(n_ij, c):
    nij_json = {"n_ij": n_ij}

    with open("./results/n_ij/" + str(c) + ".json", "w") as write_file:
        json.dump(nij_json, write_file, cls=NumpyArrayEncoder)


def create_nij(daily_period, weekly_period, remainder_days, phl, ph):

    cc = 0

    if weekly_period >= daily_period:

        if remainder_days == 0:
            c_range = phl - daily_period
        else:
            c_range = phl - daily_period - 1

        # for c in range(c_range):
        for c in range(phl - 1):

            if remainder_days == 0:
                if c < c_range - 1:
                    cc = 0
                else:
                    cc = cc + 1
            else:
                if c < c_range:
                    cc = 0
                else:
                    cc = cc + 1

            print(c)

            n_ij = np.zeros((i, j, daily_period * 2 - cc))

            t = len(n_ij[0][0])
            print("t=", t)

            n_ij[:, :, 0:daily_period] = n_ij_original[:, :, ph[c, 0] : ph[c, 1]]

            cd = daily_period
            for cw in range(c + 1, c + daily_period + 1 - cc):
                print("cd=", cd)
                n_ij[:, :, cd] = np.average(
                    n_ij_original[:, :, ph[cw, 0] : ph[cw, 1]], axis=2
                ).astype(int)
                cd = cd + 1

            save_json(n_ij, c)

        if (
            remainder_days != 0
        ):  # we need to create another n_ij for the remaining days of the last period.
            # Can we add the remainer to the estimated remaining of the unvacinated people to run multi period?
            # Load the unmet demand from the privious runs. Then, estimate the remaining number of days to run.
            # If > computational power, rum multiperiod for computational power days, else, run just one multiperiod

            c = c + 1
            print(c)

            n_ij = np.zeros((i, j, remainder_days))
            t = len(n_ij[0][0])
            print("last_remaining_t=", t)

            n_ij[:, :, 0:remainder_days] = n_ij_original[:, :, ph[c, 0] : ph[c, 1]]

            save_json(n_ij, c)

    # %% if weekly_period < daily_period:
    if weekly_period < daily_period:

        c_range = phl - 1

        if remainder_days == 0:
            cw_range = c_range
        else:
            cw_range = c_range + 1

        for c in range(c_range):

            print(c)

            n_ij = np.zeros((i, j, daily_period + c_range - c))

            t = len(n_ij[0][0])
            print("t=", t)

            n_ij[:, :, 0:daily_period] = n_ij_original[:, :, ph[c, 0] : ph[c, 1]]

            cd = daily_period
            for cw in range(c + 1, cw_range):
                print("cd=", cd)
                n_ij[:, :, cd] = np.average(
                    n_ij_original[:, :, ph[cw, 0] : ph[cw, 1]], axis=2
                ).astype(int)
                cd = cd + 1

            save_json(n_ij, c)

        if (
            remainder_days != 0
        ):  # we need to create another n_ij for the remaining days of the last period.
            # Can we add the remainer to the estimated remaining of the unvacinated people to run multi period?
            # Load the unmet demand from the privious runs. Then, estimate the remaining number of days to run.
            # If > computational power, rum multiperiod for computational power days, else, run just one multiperiod

            c = c + 1
            print(c)

            n_ij = np.zeros((i, j, remainder_days))
            t = len(n_ij[0][0])
            print("last_remaining_t=", t)

            n_ij[:, :, 0:remainder_days] = n_ij_original[:, :, ph[c, 0] : ph[c, 1]]

            save_json(n_ij, c)
