import warnings

import numpy as np

from melp.taft.utils.cosmic import cosmic_tof_correction

# ---------------------------------------
# using angle distribution to correct for TOF
# Equation from Christian Graf
# --> 7.1 in Bachelor Thesis

def alpha_from_z(z: float) -> float:
    # TODO: This is only approximated!
    return np.pi - (((25 / 328.8) * abs(z)) / 180) * np.pi
    # return np.pi + (((25 / 328.8) * z) / 180) * np.pi


def tof_z_graf(z: list) -> float:
    tile_length = 5. * (10 ** (-3))  # m
    c = 299792458  # m/s
    beta = (19 / 180) * np.pi  # deg
    alpha = alpha_from_z(z[2])

    time = (tile_length / (2 * c)) * (1 + np.tan(alpha) ** 2 + np.tan(beta) ** 2)
    return time * (10 ** 9)


def tof_z_new(xyz: list) -> float:
    tile_length = 5. * (10 ** (-3))  # m
    c = 299792458  # m/s
    beta = (19 / 180) * np.pi  # deg
    alpha = alpha_from_z(xyz[2])

    length_z = (tile_length / 2) * (np.cos(alpha))
    tof = - (length_z / c) * (10 ** 9)  # ns
    return tof


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
            if tof_mode == "advanced_graf":
                tof_time = tof_z_graf(__detector__.TileDetector.tile[id_index].pos)
                if station_offset == 200000:
                    dt_tmp += tof_time
                else:
                    dt_tmp -= tof_time
            elif tof_mode == "advanced_new":
                tof_time = tof_z_new(__detector__.TileDetector.tile[id_index].pos)
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


# ---------------------------------------------
def cosmic_correction_z(__detector__, **kwargs):
    print("analyzing cosmics file #1")
    kwargs["station"] = 1
    hist_z_1 = cosmic_tof_correction(kwargs["cosmic_file"], __detector__, **kwargs)
    correction_1 = np.median(hist_z_1)
    print("analyzing cosmics file #2")
    kwargs["station"] = 2
    hist_z_2 = cosmic_tof_correction(kwargs["cosmic_file"], __detector__, **kwargs)
    correction_2 = np.median(hist_z_2)
    print("done")
    phi = __detector__.TileDetector.column_ids(0, 200000)

    for p in range(len(phi)):
        corr = []
        row = __detector__.TileDetector.row_ids(p, 200000)
        for i in range(len(row)):
            corr.append(correction_1 * i)
        for i in range(len(row)):
            __detector__.TileDetector.tile[row[i]].dt_cal -= corr[i]

    for p in range(len(phi)):
        corr = []
        row = __detector__.TileDetector.row_ids(p, 300000)
        for i in range(len(row)):
            corr.append(correction_2 * i)
        for i in range(len(row)):
            __detector__.TileDetector.tile[row[i]].dt_cal -= corr[i]
