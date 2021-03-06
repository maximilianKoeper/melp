import numpy as np
from scipy.optimize import minimize


# ---------------------------------------
def loop_correction_phi(detector, dt_phi_rel: dict, station: int):
    print("*Simple correction phi (sum loop)")
    for z in range(len(detector.TileDetector.row_ids(0, station))):
        tile_ids = detector.TileDetector.column_ids(z, station)

        sum_dt_column = 0.
        for id_index in tile_ids:
            sum_dt_column += dt_phi_rel[id_index]

        sum_dt_column /= len(tile_ids)

        for id_index in tile_ids:
            dt_phi_rel[id_index] -= sum_dt_column

    return dt_phi_rel


# ---------------------------------------
# OLD VERSION
# def loop_correction_phi(detector, dt_phi_rel: dict, station: int):
#    for z in range(len(detector.TileDetector.row_ids(0, station))):
#        tile_ids = detector.TileDetector.column_ids(z, station)
#
#        sum_dt_column = 0.
#        number_unfilled_dt = 0
#        for id_index in tile_ids:
#            try:
#                sum_dt_column += dt_phi_rel[id_index]
#            except KeyError:
#                number_unfilled_dt += 1
#                continue
#
#        sum_dt_column /= (len(tile_ids) - number_unfilled_dt)
#
#        for id_index in tile_ids:
#            try:
#                dt_phi_rel[id_index] -= sum_dt_column
#            except KeyError:
#                continue
#
#    return dt_phi_rel

# ---------------------------------------
def single_phi_min_func(dt_corr: float, phi_pos: int, dt_phi: dict, dt_z: dict, detector, penalties: list) -> float:
    t1_tmp = 0
    if detector.TileDetector.getNeighbour(phi_pos, "l") is not False:
        id_l_tmp = detector.TileDetector.getNeighbour(phi_pos, "l")
        id_lu_tmp = detector.TileDetector.getNeighbour(id_l_tmp, "u")

        t1_tmp = dt_z[id_l_tmp]
        t1_tmp -= dt_phi[id_l_tmp]
        t1_tmp += dt_corr
        t1_tmp -= dt_z[id_lu_tmp]

    t2_tmp = 0
    if detector.TileDetector.getNeighbour(phi_pos, "r") is not False:
        id_r_tmp = detector.TileDetector.getNeighbour(phi_pos, "r")
        id_u_tmp = detector.TileDetector.getNeighbour(phi_pos, "u")

        t2_tmp = dt_z[phi_pos]
        t2_tmp += dt_phi[id_r_tmp]
        t2_tmp -= dt_z[id_u_tmp]
        t2_tmp -= dt_corr

    z = detector.TileDetector.tile[phi_pos].column()
    if detector.TileDetector.tile[phi_pos].id >= 300000:
        station = 300000
    elif detector.TileDetector.tile[phi_pos].id < 300000:
        station = 200000

    tile_ids = detector.TileDetector.column_ids(z, station)

    sum_dt_column = 0.
    for id_index in tile_ids:
        if id_index is not phi_pos:
            sum_dt_column += dt_phi[id_index]
        else:
            sum_dt_column += dt_corr

    dt_tmp = penalties[0] * (t1_tmp ** 2 + t2_tmp ** 2)
    dt_tmp += penalties[1] * sum_dt_column ** 2
    dt_tmp += penalties[2] * (dt_corr - dt_phi[phi_pos]) ** 2

    return dt_tmp * 1e5


# ---------------------------------------
def loop_adv_correction_phi(detector, dt_phi_rel: dict, dt_z_rel: dict, station: int, penalties: list):
    print("*optimizing relative phi offsets (single mode) | station: ", station)
    print("  Settings:\n   -> Penalty for measured offset:", penalties[2])
    print("   -> Penalty for column offset:", penalties[1])
    print("   -> Penalty for neighboring loop:", penalties[0])

    dt_phi_rel_truth = dt_phi_rel.copy()
    for z in range(len(detector.TileDetector.row_ids(0, station))):
        for phi_id in detector.TileDetector.column_ids(z, station):
            res = minimize(single_phi_min_func, x0=dt_phi_rel_truth[phi_id] + 0.0005,
                           args=(phi_id, dt_phi_rel_truth, dt_z_rel, detector, penalties),
                           tol=1e-8)
            dt_phi_rel[phi_id] = res.x[0]


# ---------------------------------------
def multi_phi_min_func_single(dt_corr: float, phi_pos: int, dt_phi: dict, dt_z: dict, detector,
                              penalties: list) -> float:
    t1_tmp = 0
    # phi_pos = detector.TileDetector.tile[id].row()
    if detector.TileDetector.getNeighbour(phi_pos, "l") is not False:
        id_l_tmp = detector.TileDetector.getNeighbour(phi_pos, "l")
        id_lu_tmp = detector.TileDetector.getNeighbour(id_l_tmp, "u")

        t1_tmp = dt_z[id_l_tmp]
        t1_tmp -= dt_phi[id_l_tmp]
        t1_tmp += dt_corr
        t1_tmp -= dt_z[id_lu_tmp]

    t2_tmp = 0
    if detector.TileDetector.getNeighbour(phi_pos, "r") is not False:
        id_r_tmp = detector.TileDetector.getNeighbour(phi_pos, "r")
        id_u_tmp = detector.TileDetector.getNeighbour(phi_pos, "u")

        t2_tmp = dt_z[phi_pos]
        t2_tmp += dt_phi[id_r_tmp]
        t2_tmp -= dt_z[id_u_tmp]
        t2_tmp -= dt_corr

    dt_tmp = penalties[0] * (t1_tmp ** 2 + t2_tmp ** 2)
    # dt_tmp += penalties[1] * sum_dt_column ** 2
    dt_tmp += penalties[2] * (dt_corr - dt_phi[phi_pos]) ** 2

    return dt_tmp * 1e3


def multi_phi_min_func(dt_corr: list, z: int, station: int, dt_phi: dict, dt_z: dict, detector,
                       penalties: list) -> float:
    sum_1 = 0
    sum_2 = 0

    for ids in detector.TileDetector.column_ids(z, station):
        phi_pos = detector.TileDetector.tile[ids].row()
        sum_1 += multi_phi_min_func_single(dt_corr[phi_pos], ids, dt_phi, dt_z, detector, penalties)
        sum_2 += dt_corr[phi_pos]

    return sum_1 ** 2 + penalties[1] * sum_2 ** 2


# ---------------------------------------
def loop_adv_full_correction_phi(detector, dt_phi_rel: dict, dt_z_rel: dict, station: int, penalties: list):
    print("*optimizing relative phi offsets (Column by column) | station: ", station)
    print("  Settings:\n   -> Penalty for measured offset:", penalties[2])
    print("   -> Penalty for column offset:", penalties[1])
    print("   -> Penalty for neighboring loop:", penalties[0])

    dt_phi_rel_start = dt_phi_rel.copy()
    for z in range(len(detector.TileDetector.row_ids(0, station))):
        dt_corr = []
        for phi in detector.TileDetector.column_ids(z, station):
            dt_corr.append(dt_phi_rel_start[phi])

        # print(multi_phi_min_func(dt_corr, z, station, dt_phi_rel, dt_z_rel, detector, penalties))
        res = minimize(multi_phi_min_func, x0=dt_corr,
                       args=(z, station, dt_phi_rel_start, dt_z_rel, detector, penalties), tol=1e-8)
        # print(res)
        for phi in detector.TileDetector.column_ids(z, station):
            dt_phi_rel[phi] = res.x[detector.TileDetector.tile[phi].row()]


# ---------------------------------------
def single_z_min_func(dt_corr: float, current_id: int, dt_phi: dict, dt_z: dict, detector, penalties: list) -> float:
    id_r_tmp = detector.TileDetector.getNeighbour(current_id, "r")
    id_u_tmp = detector.TileDetector.getNeighbour(current_id, "u")

    t1_tmp = dt_phi[current_id]
    t1_tmp += dt_z[id_u_tmp]
    t1_tmp -= dt_phi[id_r_tmp]
    t1_tmp -= dt_corr

    id_d_tmp = detector.TileDetector.getNeighbour(current_id, "d")
    id_dr_tmp = detector.TileDetector.getNeighbour(id_d_tmp, "r")

    t2_tmp = dt_corr
    t2_tmp -= dt_phi[id_dr_tmp]
    t2_tmp -= dt_z[id_d_tmp]
    t2_tmp += dt_phi[id_d_tmp]

    dt_tmp = penalties[0] * (t1_tmp ** 2 + t2_tmp ** 2)
    dt_tmp += penalties[1] * (dt_corr - dt_z[current_id]) ** 2

    return dt_tmp * 1e5


# ---------------------------------------
def loop_adv_correction_z(detector, dt_phi_rel: dict, dt_z_rel: dict, station: int, penalties: list):
    print("*optimizing relative z offsets")
    print("  Settings:\n   -> Penalty for measured offset:", penalties[1])
    print("   -> Penalty for neighboring loop:", penalties[0])

    dt_z_rel_truth = dt_z_rel.copy()
    for z in range(len(detector.TileDetector.row_ids(0, station)) - 1):
        for phi_id in detector.TileDetector.column_ids(z, station):
            res = minimize(single_z_min_func, x0=dt_z_rel_truth[phi_id] + 0.0001,
                           args=(phi_id, dt_phi_rel, dt_z_rel_truth, detector, penalties),
                           tol=1e-8)
            print(res.x[0] - dt_z_rel[phi_id])
            # print(res)
            dt_z_rel[phi_id] = res.x[0]
