from inputs.inputs import *

with open("./inputs/input_planning.json", "r") as read_file:
    input_planning = json.load(read_file)

# Assign the value of each Json's key to itself.
for key in input_planning:
    if isinstance(input_planning[key], list):
        globals()[key] = np.asarray(input_planning[key])
    else:
        globals()[key] = input_planning[key]


def sol_to_json(
    c_main, is_within_target_period, t_fn, m, x_ij, y_ii, z_ii, nd_ij, vrh_i, phi_i, s_t
):
    t = t_fn
    x_ij_sv = np.zeros((i, j, t))
    y_ii_sv = np.zeros((i, i, t))
    z_ii_sv = np.zeros((i, i, t))
    nd_ij_sv = np.zeros((i, j, t))
    vrh_i_sv = np.zeros((i, t))
    phi_i_sv = np.zeros((i, t))
    s_t_sv = np.zeros((t))
    for ci in range(i):
        for cii in range(i):
            for ct in range(t):
                y_ii_sv[ci, cii, ct] = y_ii[ci, cii, ct].solution_value
                z_ii_sv[ci, cii, ct] = z_ii[ci, cii, ct].solution_value

    for ci in range(i):
        for cj in range(j):
            for ct in range(t):
                x_ij_sv[ci, cj, ct] = x_ij[ci, cj, ct].solution_value
                nd_ij_sv[ci, cj, ct] = nd_ij[ci, cj, ct].solution_value
                phi_i_sv[ci, ct] = phi_i[ci, ct].solution_value
                vrh_i_sv[ci, ct] = vrh_i[ci, ct].solution_value
                s_t_sv[ct] = s_t[ct].solution_value

    of = m.solution.objective_value
    sol_time = m.solve_details.time
    gap_reached = m.solve_details.mip_relative_gap

    # Calculating alpha_used
    alpha_used = np.zeros((t))

    for ct in range(t):
        alpha_temp = 0
        for ci in range(i):
            for cii in range(i):
                if ci != cii:
                    alpha_temp = alpha_temp + z_ii_sv[ci, cii, ct] * t_ii[ci, cii]
        alpha_used[ct] = alpha_temp

    # CALC RISK
    cooef_risk = np.zeros((i, j, t))

    for ct in range(t):
        for ci in range(i):
            for cj in range(j):
                makhraj = 0
                for cii in range(i):
                    makhraj = makhraj + n_ij[cii, cj, ct]
                if makhraj == 0:
                    kasr = 0
                else:
                    kasr = n_ij[ci, cj, ct] / makhraj

                cooef_risk[ci, cj, ct] = (w_j[cj] + gamma * kasr) * (
                    nd_ij_sv[ci, cj, ct] - x_ij_sv[ci, cj, ct]
                )

    risk = np.sum(cooef_risk, axis=(0, 1))

    # Calcs Unmet Demand
    unmet_demand = np.zeros((i, j))
    # unmet_demand = nd_ij_sv[:,:,t-1] - x_ij_sv[:,:,t-1]
    if t >= daily_period and is_within_target_period == 1:
        unmet_demand = (
            nd_ij_sv[:, :, int(daily_period - 1)] - x_ij_sv[:, :, int(daily_period - 1)]
        )
    else:
        unmet_demand = nd_ij_sv[:, :, int(t - 1)] - x_ij_sv[:, :, int(t - 1)]
    # for ct in range(t):
    # unmet_demand += nd_ij_sv[:,:,ct] - x_ij_sv[:,:,ct]

    # Save to JSON
    output_json = {
        "i": i,
        "j": j,
        "t": t,
        "daily_period": daily_period,
        "weekly_period": weekly_period,
        "remainder_days": remainder_days,
        "phl": phl,
        "ph": ph,
        "select": select,
        "days": days,
        "s": s,
        "k": k,
        "sol_time": sol_time,
        "alpha": alpha,
        "theta": theta,
        "epsilon": epsilon,
        "gamma": gamma,
        "of": of,
        "risk": risk,
        "gap_reached": gap_reached,
        "alpha_used": alpha_used,
        "x_ij_sv": x_ij_sv,
        "y_ii_sv": y_ii_sv,
        "z_ii_sv": z_ii_sv,
        "nd_ij_sv": nd_ij_sv,
        "vrh_i_sv": vrh_i_sv,
        "phi_i_sv": phi_i_sv,
        "s_t_sv": s_t_sv,
        "n_ij": n_ij,
        "nu_ii": nu_ii,
        "unmet_demand": unmet_demand,
        "is_within_target_period": is_within_target_period,
        "single_period": 0,
    }

    # Writing json to file
    path_file = "./results/solutions/" + str(c_main) + ".json"
    with open(path_file, "w") as write_file:
        json.dump(output_json, write_file, cls=NumpyArrayEncoder)

    # return x_ij_sv, y_ii_sv, z_ii_sv, nd_ij_sv, phi_i_sv, vrh_i_sv, s_t_sv, alpha_used, risk, unmet_demand
