from inputs.inputs import *
from functions.sol_to_json import *


def model_fn(is_within_target_period, n_ij_fn, c_main, t_fn):

    n_ij = n_ij_fn.copy()
    t = t_fn
    # gap_stop = 0.99  # endregion

    m = Model(name="vaccine")

    # region DVs
    phi_i = m.integer_var_matrix(
        keys1=i, keys2=t, lb=0, name="phi"
    )  # number of vaccine packages allocated to demand point i

    x_ij = m.integer_var_cube(
        keys1=i, keys2=j, keys3=t, lb=0, name="x"
    )  # he number of vaccine units administrated

    y_ii = m.continuous_var_cube(
        keys1=i, keys2=i, keys3=t, lb=0, name="y"
    )  # The number of vaccine units transshipped demand pint i ∈ D to demand point j ̸= i ∈ D.

    z_ii = m.binary_var_cube(
        keys1=i, keys2=i, keys3=t, name="z"
    )  # The binary decision variable to define if a transshipment form demand pint i ∈ D to demand point j ̸= i ∈ D is occurred.

    nd_ij = m.integer_var_cube(keys1=i, keys2=j, keys3=t, lb=0, name="nd")

    vrh_i = m.integer_var_matrix(keys1=i, keys2=t, lb=0, name="vrh")

    s_t = m.integer_var_list(keys=t, name="st")

    p_i = m.binary_var_matrix(keys1=i, keys2=t, name="p")  # T_in in the model
    p_prime_i = m.binary_var_matrix(keys1=i, keys2=t, name="p_prime_i")

    # u = m.binary_var(name='u') ????????????????

    # endregion

    # region Constraints
    # C2
    for ct in range(t):
        for c1 in range(i):
            if ct == 0:
                m.add_constraint(
                    (m.sum(x_ij[c1, c2, ct] for c2 in range(j)))
                    == (
                        k * phi_i[c1, ct]
                        - vrh_i[c1, ct]
                        + m.sum(y_ii[c4, c1, ct] for c4 in range(i) if c4 != c1)
                        - m.sum(y_ii[c1, c3, ct] for c3 in range(i) if c3 != c1)
                    )
                )
            else:
                m.add_constraint(
                    (m.sum(x_ij[c1, c2, ct] for c2 in range(j)))
                    == (
                        k * phi_i[c1, ct]
                        + (vrh_i[c1, ct - 1] - vrh_i[c1, ct])
                        + m.sum(y_ii[c4, c1, ct] for c4 in range(i) if c4 != c1)
                        - m.sum(y_ii[c1, c3, ct] for c3 in range(i) if c3 != c1)
                    )
                )

    # C3
    for ct in range(t):
        for ci in range(i):
            for cj in range(j):
                if ct == 0:
                    m.add_constraint(nd_ij[ci, cj, ct] == n_ij[ci, cj, ct])
                else:
                    m.add_constraint(
                        nd_ij[ci, cj, ct]
                        == n_ij[ci, cj, ct]
                        + m.sum(
                            n_ij[ci, cj, ct2] - x_ij[ci, cj, ct2] for ct2 in range(ct)
                        )
                    )

    # C4
    for ct in range(t):
        if ct == 0:
            m.add_constraint(s_t[ct] == s)
        else:
            # m.add_constraint(s_t[ct] == s - m.sum(vrh_i[ci,ct] for ci in range(i)))
            m.add_constraint(s_t[ct] == s + m.sum(vrh_i[ci, ct - 1] for ci in range(i)))

    # C5
    for ct in range(t):
        m.add_constraint((m.sum(phi_i[ci, ct] for ci in range(i))) <= s_t[ct])

    # C6
    for ct in range(t):
        for ci in range(i):
            m.add_constraint((m.sum(x_ij[ci, cj, ct] for cj in range(j))) <= c_i[ci])

    # C7
    for ct in range(t):
        m.add_constraint(
            (
                m.sum(
                    z_ii[ci, cii, ct] * t_ii[ci, cii]
                    for ci in range(i)
                    for cii in range(i)
                    if ci != cii
                )
            )
            <= alpha
        )

    # C8,9,10
    for ct in range(t):
        for ci in range(i):
            m.add_constraint(
                (m.sum(z_ii[ci, cii, ct] for cii in range(i))) <= p_i[ci, ct] * M
            )
            m.add_constraint(
                (m.sum(z_ii[cii, ci, ct] for cii in range(i))) <= p_prime_i[ci, ct] * M
            )
            m.add_constraint(p_i[ci, ct] + p_prime_i[ci, ct] <= 1)

    # C11
    # z_ii = 0 if hospital has no demand or it is not within the 5 closest hospitals.
    for ct in range(t):
        for ci in range(i):
            for cii in range(i):
                if nu_ii[ci, cii] == 0:
                    m.add_constraint(z_ii[ci, cii, ct] == 0)
                    m.add_constraint(y_ii[ci, cii, ct] == 0)

    for ct in range(t):
        for ci in range(i):
            for cii in range(i):
                if cii != ci:
                    m.add_constraint(z_ii[ci, cii, ct] <= nu_ii[ci, cii])

    # C12
    for ct in range(t):
        for ci in range(i):
            for cj in range(j):
                m.add_constraint(x_ij[ci, cj, ct] <= nd_ij[ci, cj, ct])

    # C13
    for ct in range(t):
        for ci in range(i):
            for cii in range(i):
                if cii != ci:
                    m.add_constraint(y_ii[ci, cii, ct] <= M * z_ii[ci, cii, ct])

    print()

    # nu
    theta_prime = np.max(
        [t_ii[c1, c2] * nu_ii[c1, c2] for c1 in range(i) for c2 in range(i)]
    )

    # endregion

    # region OF
    m.minimize(
        1
        * (
            m.sum(
                (
                    w_j[c5]
                    + (
                        gamma
                        * (
                            n_ij[c1, c5, ct]
                            / (sum(n_ij[c3, c5, ct] for c3 in range(i)))
                            if (sum(n_ij[c3, c5, ct] for c3 in range(i))) != 0
                            else 0
                        )
                    )
                )
                * (nd_ij[c1, c5, ct] - x_ij[c1, c5, ct])
                for c1 in range(i)
                for c5 in range(j)
                for ct in range(t)
            )
        )
        + epsilon
        * (
            m.sum(
                y_ii[c6, c7, ct2] * t_ii[c6, c7]
                for c6 in range(i)
                for c7 in range(i)
                for ct2 in range(t)
                if c7 != c6
            )
        )
        + epsilon * theta * m.sum(vrh_i[c8, ct3] for c8 in range(i) for ct3 in range(t))
    )

    # endregion

    # region Model Settings and Solve
    m.print_information()

    # m.parameters.mip.tolerances.mipgap = 1e-5
    m.parameters.mip.tolerances.mipgap = gap_stop  # To stop when the gap gets to 1%

    # m.parameters.timelimit(120) #120 seconds is the allowed runtime of the cplex
    m.parameters.timelimit(run_time)

    msol = m.solve(url=url, key=key, log_output=True)
    assert msol, "solve failed"
    m.report()

    sol_to_json(
        c_main,
        is_within_target_period,
        t,
        m,
        x_ij,
        y_ii,
        z_ii,
        nd_ij,
        vrh_i,
        phi_i,
        s_t,
    )

    # endregion
