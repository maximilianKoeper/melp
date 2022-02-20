import numpy as np
import scipy.optimize as opt
from melp.taft.utils.cosmic import get_cosmic_data_from_file


def correct_z_two_event(detector, cosmic_station, **kwargs) -> np.array:
    print("*Two event correction")
    a, b = get_cosmic_data_from_file(kwargs.get("cosmic_file"), detector, cosmic_station, **kwargs)
    print("\n  -> Cosmic Muons: ", len(a))

    args = np.zeros(kwargs["cosmic_n_modes"]*2) + 0.0001
    popt, cov = opt.curve_fit(fit_func, b, a, p0=args, method="lm")

    station_offset = 200000
    if cosmic_station == 2:
        station_offset += 100000
    for row in range(len(detector.TileDetector.column_ids(0, station_offset))):
        for column in range(len(detector.TileDetector.row_ids(0, station_offset))):
            current_tile_id = detector.TileDetector.id_from_row_col(row=row, column=column, station_offset=station_offset)
            timing_correction = calibration_correction_z((column, row), *popt)
            timing_correction -= calibration_correction_z((0, 0), *popt)

            detector.TileDetector.tile[current_tile_id].update_calibration(timing_correction)

    return popt


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
def calibration_correction_both(x: tuple, *args) -> float or np.array:
    # v_fourier_func = np.vectorize(fourier_func, excluded=['c', 's', 'n'])
    z, phi = x
    # print(phi, z)
    n = int(len(args) / 2)
    result_tmp1 = fourier_func(2 * np.pi * (np.array(phi) / 56), *args[0:n + 1])
    result_tmp2 = fourier_func(np.pi * (np.array(z) / 52), *args[n:2 * n])
    result = np.multiply(result_tmp1, result_tmp2)
    # return result
    return result


def calibration_correction_z(x: tuple, *args) -> float or np.array:
    # v_fourier_func = np.vectorize(fourier_func, excluded=['c', 's', 'n'])
    z, phi = x
    # print(phi, z)
    result = fourier_func(np.pi * (np.array(z) / 52), *args)
    #result += args[0]*z
    return result


def fit_func(x: tuple, *args) -> float or np.array:
    #if len(args) % 4 != 0:
    #    raise ValueError
    # z1, phi1, z2, phi2 = x
    z2, phi2, z1, phi1 = x
    return calibration_correction_z((z1, phi1), *args) - calibration_correction_z((z2, phi2), *args)
