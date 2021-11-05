from inputs.inputs import *

# with open("./results/solutions/" + str(25) + ".json", "r") as read_file:
#     result_last_day = json.load(read_file)
# n_ij = np.asarray(result_last_day["n_ij"])
# before = np.sum(n_ij)
# max_before = np.max(n_ij)


def aug_nij_fn(n_ij):
    d_mc = np.sum(n_ij, axis=1)  # Demand of each medical center
    ur_mc = d_mc / c_i  # Utilization rate of medical centers

    for ci in range(i):
        if d_mc[ci] > c_i[ci]:
            exss_d_i = d_mc[ci] - c_i[ci]
            rem_n_ij = np.zeros((j))

            rem_c_i = c_i[ci]
            for cj in range(j - 1, -1, -1):
                if n_ij[ci, cj] < rem_c_i:
                    rem_n_ij[cj] = n_ij[ci, cj]
                    rem_c_i = rem_c_i - n_ij[ci, cj]
                else:
                    rem_n_ij[cj] = rem_c_i
                    rem_c_i = 0

            # Array to get the excess demand of medical center fro each risk category
            exss_d_ij = np.zeros((j))
            exss_d_ij = n_ij[ci, :] - rem_n_ij
            # n_ij[ci, :] = n_ij[ci, :] - exss_d_ij

            # First, we try to redistribute this amount to the adjacent MCs when they use less than 50% of their capacity
            adj_mc = np.where(nu_ii[ci] == 1)[0]
            adj_mc_less50 = [xmc for xmc in adj_mc if ur_mc[xmc] < 0.5]
            ur_adj_mc = [ur_mc[xmc] for xmc in adj_mc_less50]
            # Sort based on the utilization rate
            adj_mc_less50_sorted = [
                xmc for _, xmc in sorted(zip(ur_adj_mc, adj_mc_less50))
            ]

            prob_select = -1 * np.sort(-1 * np.subtract(1, ur_adj_mc))
            prob_select = prob_select / np.sum(prob_select)

            while np.sum(exss_d_ij) != 0 and len(adj_mc_less50_sorted) != 0:
                nbr_rcv = np.random.choice(adj_mc_less50_sorted, 1, p=prob_select)[0]

                while (
                    np.sum(n_ij[nbr_rcv, :]) < c_i[nbr_rcv] and np.sum(exss_d_ij) != 0
                ):
                    for cjj in range(j - 1, -1, -1):
                        if exss_d_ij[cjj] != 0:
                            exss_d_ij[cjj] = exss_d_ij[cjj] - 1
                            n_ij[ci, cjj] = n_ij[ci, cjj] - 1
                            n_ij[nbr_rcv, cjj] = n_ij[nbr_rcv, cjj] + 1
                index_del = np.where(adj_mc_less50_sorted == nbr_rcv)
                adj_mc_less50_sorted = np.delete(adj_mc_less50_sorted, index_del)
                prob_select = np.delete(prob_select, index_del)
                prob_select = prob_select / np.sum(prob_select)
    return n_ij

    # n_ij[ci, :] = n_ij[ci, :] + exss_d_ij

    # print('________________')


# n_ij = aug_nij_fn(n_ij)
# print(before)
# print(max_before)
# print(np.sum(n_ij))
# print(np.max(n_ij))
