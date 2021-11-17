import numpy as np


# ---------------------------------------
# using angle distribution to correct for TOF
# Equation from Christian Graf
# --> 7.1 in Bachelor Thesis

def alpha_from_z(z: float) -> float:
    # TODO: This is only approximated!
    return (((20 / 328.8) * z) / 180) * np.pi


def tof_z(z: list) -> float:
    # TODO: check equation (seems a bit off)
    tile_length = 5. * (10 ** (-3))  # m
    c = 299792458  # m/s
    beta = (19 / 180) * np.pi  # deg
    alpha = alpha_from_z(z[2])

    # offset = 0.003
    offset = 0.

    time = (tile_length / (2 * c)) * (1 + np.tan(alpha) ** 2 + np.tan(beta) ** 2)
    return time * (10 ** 9) - offset

