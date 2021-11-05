from inputs.inputs import *

with open("./inputs/input_planning.json", "r") as read_file:
    input_planning = json.load(read_file)

# Assign the value of each Json's key to itself.
for key in input_planning:
    if isinstance(input_planning[key], list):
        globals()[key] = np.asarray(input_planning[key])
    else:
        globals()[key] = input_planning[key]


def model_singleP_fn(n_ij_fn, c_main, t_fn):

    n_ij = n_ij_fn.copy()
    t = t_fn

    m = Model(name="vaccine_singleP")

    phi_i = m.integer_var_list(
        keys=i, lb=0, name="phi"
    )  # number of vaccine packages allocated to demand point i

    x_ij = m.integer_var_matrix(
        keys1=i, keys2=j, lb=0, name="x"
    )  # he number of vaccine units administrated

    y_ii = m.continuous_var_matrix(
        keys1=i, keys2=i, lb=0, name="y"
    )  # The number of vaccine units transshipped demand pint i ∈ D to demand point j ̸= i ∈ D.

    z_ii = m.binary_var_matrix(
        keys1=i, keys2=i, name="z"
    )  # The binary decision variable to define if a transshipment form demand pint i ∈ D to demand point j ̸= i ∈ D is occurred.

    delta = m.continuous_var(name="delta")

    p_i = m.binary_var_list(keys=i, name="p")
    p_prime_i = m.binary_var_list(keys=i, name="p_prime_i")

    u = m.binary_var(name="u")

    # C2
    for c1 in range(i):
        m.add_constraint(
            (m.sum(x_ij[c1, c2] for c2 in range(j)))
            <= (
                k * phi_i[c1]
                - m.sum(y_ii[c1, c3] for c3 in range(i) if c3 != c1)
                + m.sum(y_ii[c4, c1] for c4 in range(i) if c4 != c1)
            )
        )

    # C3
    m.add_constraint(
        m.sum(x_ij[c1, c2] for c1 in range(i) for c2 in range(j))
        == m.sum(phi_i[c3] for c3 in range(i)) * k - delta
    )

    # C4
    m.add_constraint((m.sum(phi_i[c1] for c1 in range(i))) <= s)

    # C5
    for c1 in range(i):
        m.add_constraint((m.sum(x_ij[c1, c2] for c2 in range(j))) <= c_i[c1])

    # C6
    m.add_constraint(
        (
            m.sum(
                z_ii[c1, c2] * t_ii[c1, c2]
                for c1 in range(i)
                for c2 in range(i)
                if c2 != c1
            )
        )
        <= alpha
    )

    # C7,8,9
    for c1 in range(i):
        m.add_constraint((m.sum(z_ii[c1, c2] for c2 in range(i))) <= p_i[c1] * M)
        m.add_constraint((m.sum(z_ii[c2, c1] for c2 in range(i))) <= p_prime_i[c1] * M)
        m.add_constraint(p_i[c1] + p_prime_i[c1] <= 1)

    # C11
    for c1 in range(i):
        for c2 in range(j):
            m.add_constraint(x_ij[c1, c2] <= n_ij[c1, c2])

    # C12
    for c1 in range(i):
        for c2 in range(i):
            if c2 != c1:
                m.add_constraint(y_ii[c1, c2] <= M * z_ii[c1, c2])

    # C13
    # m.add_constraint(delta <= k-1)

    print()

    # Creating nu_ii for C10
    nu_ii = np.zeros(shape=(i, i))  # Create matrix nu_ij

    n_i = np.sum(n_ij, axis=1)  # Total demand of each hospital sumation over j

    for c in range(i):
        if n_i[c] != 0:
            # print(c)
            hos_with_demand_index = np.where(n_i > 0)[0]
            hos_with_demand_index = np.delete(
                hos_with_demand_index, np.where(hos_with_demand_index == c)
            )
            t_i_all_demand = t_ii[c, hos_with_demand_index]

            short_dist_index = np.argsort(t_i_all_demand)[:select]
            short_dist = t_i_all_demand[
                short_dist_index
            ]  # values of distance from the 5 hospitals from the main hospital c

            hos_js = [
                np.where(t_ii[c] == c1)[0][0] for c1 in t_i_all_demand[short_dist_index]
            ]
            for c2 in hos_js:
                nu_ii[c, c2] = 1
                nu_ii[c2, c] = 1

    # C10
    # z_ii = 0 if hospital has no demand or it is not within the 5 closest hospitals.
    for c1 in range(i):
        for c2 in range(i):
            if nu_ii[c1, c2] == 0:
                m.add_constraint(z_ii[c1, c2] == 0)
                m.add_constraint(y_ii[c1, c2] == 0)

    for c1 in range(i):
        for c2 in range(i):
            if c2 != c1:
                m.add_constraint(z_ii[c1, c2] <= nu_ii[c1, c2])

    theta_prime = np.max(
        [t_ii[c1, c2] * nu_ii[c1, c2] for c1 in range(i) for c2 in range(i)]
    )

    # Cut by delta
    m.add_constraint(
        delta
        <= (
            len(np.nonzero(np.sum(n_ij, axis=1))[0])
            / min(select, len(np.nonzero(np.sum(n_ij, axis=1))[0]))
        )
        * (k - 1)
        + M * (1 - u)
    )

    m.add_constraint(
        alpha * u
        >= alpha
        - (m.sum(t_ii[c1, c2] * y_ii[c1, c2] for c1 in range(i) for c2 in range(i)))
        - theta_prime
    )

    # KPI
    alpha_used = m.sum(
        z_ii[c1, c2] * t_ii[c1, c2] for c1 in range(i) for c2 in range(i) if c2 != c1
    )

    m.add_kpi(alpha_used, "alpha_used")

    # KPI
    risk = m.sum(
        (
            w_j[c5]
            + (
                gamma
                * (
                    n_ij[c1, c5] / (sum(n_ij[c3, c5] for c3 in range(i)))
                    if (sum(n_ij[c3, c5] for c3 in range(i))) != 0
                    else 0
                )
            )
        )
        * (n_ij[c1, c5] - x_ij[c1, c5])
        for c1 in range(i)
        for c5 in range(j)
    )

    m.add_kpi(risk, "risk")

    # OF
    m.minimize(
        1
        * (
            m.sum(
                (
                    w_j[c5]
                    + (
                        gamma
                        * (
                            n_ij[c1, c5] / (sum(n_ij[c3, c5] for c3 in range(i)))
                            if (sum(n_ij[c3, c5] for c3 in range(i))) != 0
                            else 0
                        )
                    )
                )
                * (n_ij[c1, c5] - x_ij[c1, c5])
                for c1 in range(i)
                for c5 in range(j)
            )
        )
        + epsilon
        * (
            m.sum(
                y_ii[c6, c7] * t_ii[c6, c7]
                for c6 in range(i)
                for c7 in range(i)
                if c7 != c6
            )
        )
        + epsilon * delta * theta
    )

    m.print_information()

    # m.parameters.mip.tolerances.mipgap = 1e-5
    m.parameters.mip.tolerances.mipgap = gap_stop  # To stop when the gap gets to 1%

    # m.parameters.timelimit(120) #120 seconds is the allowed runtime of the cplex
    m.parameters.timelimit(run_time)

    msol = m.solve(url=url, key=key, log_output=True)
    assert msol, "solve failed"
    m.report()

    x_ij_sv = np.zeros((i, j, t))
    y_ii_sv = np.zeros((i, i, t))
    z_ii_sv = np.zeros((i, i, t))
    phi_i_sv = np.zeros((i, t))
    for ci in range(i):
        for cii in range(i):
            for ct in range(t):
                y_ii_sv[ci, cii, ct] = y_ii[ci, cii].solution_value
                z_ii_sv[ci, cii, ct] = z_ii[ci, cii].solution_value

    for ci in range(i):
        for cj in range(j):
            for ct in range(t):
                x_ij_sv[ci, cj, ct] = x_ij[ci, cj].solution_value
                phi_i_sv[ci, ct] = phi_i[ci].solution_value

    of = m.solution.objective_value
    sol_time = m.solve_details.time
    gap_reached = m.solve_details.mip_relative_gap

    unmet_demand = n_ij - x_ij_sv[:, :, 0]

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
        "risk": risk.solution_value,
        "gap_reached": gap_reached,
        "delta": delta.solution_value,
        "alpha_used": alpha_used.solution_value,
        "x_ij_sv": x_ij_sv,
        "y_ii_sv": y_ii_sv,
        "z_ii_sv": z_ii_sv,
        "phi_i_sv": phi_i_sv,
        "n_ij": n_ij,
        "nu_ii": nu_ii,
        "unmet_demand": unmet_demand,
        "is_within_target_period": 0,
        "single_period": 1,
    }

    # Writing json to file
    path_file = "./results/solutions/" + str(c_main) + ".json"
    with open(path_file, "w") as write_file:
        json.dump(output_json, write_file, cls=NumpyArrayEncoder)
