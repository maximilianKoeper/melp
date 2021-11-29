import numpy as np
from scipy.optimize import minimize

from melp.taft.utils.misc_array_creation import create_array_from_dict


def fit_func(params: np.array, arr: np.array) -> float:
    result = 0.

    tmp_z = np.zeros((51, 56))
    for i in range(51):
        for j in range(56):
            tmp_z[i][j] = params[i + j]

    tmp_phi = np.zeros((52, 56))
    for i in range(51):
        for j in range(56):
            tmp_phi[i][j] = params[51 * 56 - 1 + i + j]

    DT_z_tmp = arr[1][:][:] + tmp_z
    DT_phi_tmp = arr[0][:][:]

    return result


def global_align(detector, dt_z: dict, dt_phi: dict, station: int) -> (dict, dict):
    dt_arr = create_array_from_dict(detector=detector, dt_dict_z=dt_z, dt_dict_phi=dt_phi, station=station)

    res = minimize(fit_func, x0=np.zeros(56 * 51 + 56 * 52), args=dt_arr, method="BFGS")

    print(res.success)
    print(res.x)
