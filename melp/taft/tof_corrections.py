import numpy as np
import warnings


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


# ---------------------------------------
def tof_correction_z(__detector__, dt_z_rel: dict, station_offset: int, tof_mode: str) -> dict:
    for phi in range(len(__detector__.TileDetector.column_ids(0, station_offset))):
        tile_ids = __detector__.TileDetector.row_ids(phi, station_offset)

        for id_index in tile_ids:
            try:
                dt_tmp = dt_z_rel[id_index]
            except KeyError:
                continue

            # TOF correction advanced
            # TODO: check tof_modes
            if tof_mode == "advanced":
                tof_time = tof_z(__detector__.TileDetector.tile[id_index].pos)
                if station_offset == 200000:
                    dt_tmp += tof_time
                else:
                    dt_tmp -= tof_time
            # TOF simple
            elif tof_mode == "simple":
                tof_time = 0.009
                if station_offset == 200000:
                    dt_tmp += tof_time
                else:
                    dt_tmp -= tof_time
            # No TOF correction
            else:
                warnings.warn("No TOF correction applied")

            dt_z_rel[id_index] = dt_tmp

    return dt_z_rel
