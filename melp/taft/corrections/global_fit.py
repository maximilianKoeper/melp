import numpy as np


# def fourier_func(x: float or list, c: list, s: list, n: int) -> float or list:
def fourier_func(x: float or list, *args) -> float or np.array:
    result = np.zeros(x.shape)

    n = int(len(args) / 2)
    c = args[0:n + 1]
    s = args[n:2 * n]
    for i in range(1, n + 1):
        # result += c[i - 1] * np.cos(x * i) + s[i - 1] * np.sin(x * i)
        result = np.add(result, c[i - 1] * np.cos(x * i) + s[i - 1] * np.sin(x * i))

    return result


# def calibration_correction(x: tuple, c_phi: list, s_phi: list, c_z: list, s_z: list, n: int) -> float or list:
def calibration_correction(x: tuple, *args) -> float or np.array:
    # v_fourier_func = np.vectorize(fourier_func, excluded=['c', 's', 'n'])
    z, phi = x
    # print(phi, z)
    n = int(len(args) / 2)
    # result_tmp1 = fourier_func(2*np.pi*(np.array(phi)/56), *args[0:n+1])
    result_tmp2 = fourier_func(np.pi * (np.array(z) / 52), *args[n:2 * n])
    # result = np.multiply(result_tmp1, result_tmp2)
    # return result
    return result_tmp2


def fit_func(x: tuple, *args) -> float or np.array:
    if len(args) % 4 != 0:
        raise ValueError
    z1, phi1, z2, phi2 = x
    return calibration_correction((z1, phi1), *args) - calibration_correction((z2, phi2), *args)
