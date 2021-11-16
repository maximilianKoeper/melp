import ROOT

import numpy as np
import warnings

from melp.libs.misc import *
import melp

# ---------------------------------------------------------------------
#  Define global variables and functions to select Detector
# ---------------------------------------------------------------------

__detector__ = None


def select(selection: melp.Detector):
    global __detector__
    __detector__ = selection


def info():
    global __detector__
    print(__detector__)


# ---------------------------------------------------------------------
# calibrate function is the only function the user should call
# -> options for calibrating should be passed with **kwargs
# ---------------------------------------------------------------------

def calibrate(filename: str, **kwargs):
    # TODO: tidy up and streamline options
    global __detector__
    if __detector__ is None:
        print("ERROR: Detector not selected")
        return

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")

    # Get dict with histogramms with dt data
    histogram = fill_dt_histos(ttree_mu3e)

    # calculate residuals to truth
    # and dt between tiles (median of the histogram)
    dt_z, dt_phi, resid_z, resid_phi = get_median_from_hist(histogram, resid=True)

    # prealign in z direction
    if "tof" in kwargs.keys():
        cal_station_1_z, cal_station_2_z = pre_align_z(dt_z, tof_mode=kwargs["tof"])
    else:
        cal_station_1_z, cal_station_2_z = pre_align_z(dt_z)

    # prealign in phi direction
    cal_station_1_phi, cal_station_2_phi = pre_align_phi(dt_phi)

    # correction for phi loop (going around the loop should result in sum dt = 0)
    loop_correction_phi_simple(cal_station_1_phi)
    loop_correction_phi_simple(cal_station_2_phi)

    # TODO: combine phi and z information

    # return difference from truth
    cal = None
    if "station" in kwargs.keys():
        if kwargs["station"] == 1:
            cal = check_cal_z(cal_station_1_z, 1)
        elif kwargs["station"] == 2:
            cal = check_cal_z(cal_station_2_z, 2)

    cal_phi = None
    if "station" in kwargs.keys():
        if kwargs["station"] == 1:
            cal_phi = check_cal_phi(cal_station_1_phi, 1)
        elif kwargs["station"] == 2:
            cal_phi = check_cal_phi(cal_station_2_phi, 2)

    return resid_z, resid_phi, cal, cal_phi


# ---------------------------------------
def loop_correction_phi_simple(cal_station: dict) -> bool:
    for z_column in cal_station:
        time_diff = cal_station[z_column][-1] - cal_station[z_column][0]
        time_diff /= len(cal_station[z_column])

        tmp_arr = cal_station[z_column]
        for entry in range(len(tmp_arr)):
            tmp_arr[int(entry)] -= time_diff * entry
        cal_station[z_column] = tmp_arr

    return True


# ---------------------------------------
# using angle distribution to correct for TOF
# Equation from Christian Graf
# --> 7.1 in Bachelor Thesis

def alpha_from_z(z: float) -> float:
    # TODO: This is only approximated!
    return (((20 / 328.8) * z) / 180) * np.pi


def tof(z: list) -> float:
    # TODO: check equation (seems a bit off)
    l = 5. * (10 ** (-3))  # m
    c = 299792458  # m/s
    beta = (19 / 180) * np.pi  # deg
    alpha = alpha_from_z(z[2])

    # offset = 0.003
    offset = 0.

    time = (l / (2 * c)) * (1 + np.tan(alpha) ** 2 + np.tan(beta) ** 2)
    return time * (10 ** 9) - offset


# ---------------------------------------

def pre_align_phi(dt_phi):
    cal_station_1 = {}
    cal_station_2 = {}
    station_offset_arr = [200000, 300000]

    # calibration station 1 and 2
    for station_offset in station_offset_arr:
        for z_column in range(52):
            dt_tiles_cal = [0.]
            for tile in range(0, 56):
                try:
                    dt_tmp = dt_phi[(station_offset + z_column * 56 + tile)]
                except KeyError:
                    dt_tmp = 0
                    warnings.warn(f"Not enough data for phi calibration")

                dt_tiles_cal.append(float(dt_tiles_cal[-1] + dt_tmp))

            if station_offset == 200000:
                cal_station_1[z_column] = np.asarray(dt_tiles_cal)
            elif station_offset == 300000:
                cal_station_2[z_column] = np.asarray(dt_tiles_cal)

    return cal_station_1, cal_station_2


# ---------------------------------------

def check_cal_phi(cal_data, station):
    cal = {}

    if station == 1:
        station_offset = 200000
    elif station == 2:
        station_offset = 300000
    else:
        raise ValueError(f"expected station=1 or 2, got station = {station}")

    for z_column in range(52):
        dt_truth = [0]
        for tile in range(0, 56):
            tile_id = station_offset + z_column * 56 + tile
            dt_tmp = (__detector__.TileDetector.tile[tile_id].dt -
                      __detector__.TileDetector.tile[__detector__.TileDetector.getNeighbour(tile_id, "up")].dt)
            dt_truth.append(dt_truth[-1] + dt_tmp)

        cal[z_column] = np.array(cal_data[z_column]) + np.array(dt_truth)

    return cal


# ---------------------------------------
# - returns the absolute timing offset of the tile to the tile nearest to the target
#
# -> returned dictionaries:
# --> cal_station_1:
#           - keys: row of tiles (0 corresponds to row starting with tile 200000)
#                                (1 corresponds to row starting with tile 200001)
#           - entries: np.array with offsets (starting with nearest tile muon gun)
# --> cal_station_2:
#           - same as 1 but for second station
#

def pre_align_z(dt_z, tof_mode: str = "None"):
    cal_station_1 = {}
    cal_station_2 = {}
    station_offset_arr = [200000, 300000]

    # calibration station 1 and 2
    for station_offset in station_offset_arr:
        for phi_row in range(55):
            dt_tiles_cal = [0.]
            for tile in range(0, 51):
                try:
                    dt_tmp = dt_z[(station_offset + phi_row + tile * 56)]
                except KeyError:
                    dt_tmp = 0
                    warnings.warn(f"Not enough data for z calibration")

                # TOF correction advanced
                # TODO: check tof_modes
                if tof_mode == "advanced":
                    tof_time = tof(__detector__.TileDetector.tile[(station_offset + phi_row + tile * 56)].pos)
                    if station_offset == 200000:
                        dt_tmp += tof_time
                    else:
                        dt_tmp -= tof_time
                elif tof_mode == "simple":
                    tof_time = 0.009
                    if station_offset == 200000:
                        dt_tmp += tof_time
                    else:
                        dt_tmp -= tof_time
                else:
                    warnings.warn("No TOF correction applied")

                dt_tiles_cal.append(float(dt_tiles_cal[-1] + dt_tmp))

            if station_offset == 200000:
                cal_station_1[phi_row] = np.asarray(dt_tiles_cal)
            elif station_offset == 300000:
                cal_station_2[phi_row] = np.asarray(dt_tiles_cal)

    return cal_station_1, cal_station_2


# ---------------------------------------

def check_cal_z(cal_data, station):
    cal = {}

    if station == 1:
        station_offset = 200000
    elif station == 2:
        station_offset = 300000
    else:
        raise ValueError(f"expected station=1 or 2, got station = {station}")

    for phi_row in range(55):
        dt_truth = [0]
        for tile in range(0, 51):
            tile_id = station_offset + phi_row + tile * 56
            dt_tmp = (__detector__.TileDetector.tile[tile_id].dt -
                      __detector__.TileDetector.tile[__detector__.TileDetector.getNeighbour(tile_id, "right")].dt)
            dt_truth.append(dt_truth[-1] + dt_tmp)

        cal[phi_row] = np.array(cal_data[phi_row]) + np.array(dt_truth)

    return cal


# ---------------------------------------
# calculates the median from the histogram dict returned by fill_dt_histogram()
# returns residuals to true dt if resid=True

def get_median_from_hist(histogram: dict, resid=False):
    resid_z = []
    dt_z = {}
    for tile_entry in histogram:
        neighbour_z_id = __detector__.TileDetector.getNeighbour(tile_entry, "right")
        if neighbour_z_id is False:
            continue
        # setting prob to 0.5 for median
        prob = np.array([0.5])
        # q for result (has to be np.array because of type casting)
        q = np.array([0.])
        if histogram[tile_entry][0].ComputeIntegral() == 0:
            # warnings.warn(f"WARNING: INTEGRAL = 0, tile: {tile_entry}")
            continue
        histogram[tile_entry][0].GetQuantiles(1, q, prob)
        dt_z[tile_entry] = q
        if resid:
            resid_dt = q - (__detector__.TileDetector.tile[neighbour_z_id].dt - __detector__.TileDetector.tile[
                tile_entry].dt)
            resid_z.append(resid_dt)

    # phi dir:
    resid_phi = []
    dt_phi = {}
    for tile_entry in histogram:

        neighbour_phi_id = __detector__.TileDetector.getNeighbour(tile_entry, "up")

        # setting prob to 0.5 for median
        prob = np.array([0.5])
        # q for result (has to be np.array because of type casting)
        q = np.array([0.])
        if histogram[tile_entry][1].ComputeIntegral() == 0:
            # print("WARNING: INTEGRAL = 0", tile_entry)
            continue
        histogram[tile_entry][1].GetQuantiles(1, q, prob)
        dt_phi[tile_entry] = q
        if resid:
            resid_dt = q - (__detector__.TileDetector.tile[neighbour_phi_id].dt - __detector__.TileDetector.tile[
                tile_entry].dt)
            resid_phi.append(resid_dt)

    return dt_z, dt_phi, resid_z, resid_phi


# ---------------------------------------
#
# Generates dictionary with ROOT TH1D Histogramms
#   -> dict[tileid] = [hist_z, hist_pih]
#

def fill_dt_histos(ttree_mu3e) -> dict:
    cluster_counter = 0

    hist_dict = {}
    nbins, lo, hi = 10000, -64, 64
    # Generating empty histos:
    for tile in __detector__.TileDetector.tile:
        histo_name_z = str(tile) + "_z"
        histo_name_phi = str(tile) + "_phi"
        hist_dict[tile] = [ROOT.TH1D(histo_name_z, 'my hist', nbins, lo, hi),
                           ROOT.TH1D(histo_name_phi, 'my hist', nbins, lo, hi)]

    for frame in range(ttree_mu3e.GetEntries()):

        # Printing status info
        if frame % 10000 == 0:
            print("Searching clusters. Progress: ", np.round(frame / ttree_mu3e.GetEntries() * 100), " % , Found: ",
                  cluster_counter, end='\r')
        ttree_mu3e.GetEntry(frame)

        # TODO: index_finder cant handle multiple events on one tile in one frame!!!
        #       --> skipping frame (looses some data)
        # Analyzing frame
        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

            # -----------------------------
            # Look for clusters in z-dir
            neighbour_z_id = __detector__.TileDetector.getNeighbour(hit_tile, "right")
            if neighbour_z_id in ttree_mu3e.tilehit_tile and neighbour_z_id is not False:
                # find associated tile hit
                hit_tile_assoc = index_finder(list(ttree_mu3e.tilehit_tile), neighbour_z_id)

                # workaround for multiple hits in the same tile
                try:
                    hit_tile_assoc = int(*hit_tile_assoc)
                except (TypeError, ValueError):
                    continue

                # calculate dt
                # TODO: TOF maybe with edep ?
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[hit_tile].dt
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[neighbour_z_id].dt
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist_dict[hit_tile][0].Fill(dt)
                cluster_counter += 1

            # -----------------------------
            # Look for clusters in phi-dir
            neighbour_phi_id = __detector__.TileDetector.getNeighbour(hit_tile, "up")
            if neighbour_phi_id in ttree_mu3e.tilehit_tile and neighbour_phi_id is not False:
                hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]
                # find associated tile hit
                hit_tile_assoc = index_finder(list(ttree_mu3e.tilehit_tile), neighbour_phi_id)

                # workaround for multiple hits in the same tile
                try:
                    hit_tile_assoc = int(*hit_tile_assoc)
                except (TypeError, ValueError):
                    continue

                # calculate dt
                # TODO: TOF maybe with edep ?
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[hit_tile].dt
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[
                    neighbour_phi_id].dt
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist_dict[hit_tile][1].Fill(dt)
                cluster_counter += 1

    print("Searching clusters. Progress: ", 100, " % , Found: ", cluster_counter)
    return hist_dict
