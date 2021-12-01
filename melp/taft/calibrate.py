import warnings

import ROOT
import numpy as np
from scipy.optimize import minimize

from melp import Detector
# fast index lookup
from melp.libs.misc import index_finder
from melp.libs.timer import Timer
from melp.taft.corrections.misc_corrections import loop_correction_phi
# different functions for calibration
from melp.taft.corrections.tof_corrections import tof_correction_z
from melp.taft.utils.root_helper import save_histo, read_histo, fill_dt_histos

# ---------------------------------------------------------------------
#  Define global variables and functions to select Detector
# ---------------------------------------------------------------------

__detector__: Detector = None


def select(selection: Detector):
    global __detector__
    __detector__ = selection


# ---------------------------------------------------------------------
# calibrate function is the only function the user should call
# -> options for calibrating should be passed with **kwargs
# ---------------------------------------------------------------------

def calibrate(**kwargs):
    # TODO: tidy up and streamline options
    global __detector__

    if __detector__ is None:
        raise RuntimeError("ERROR: Detector not selected")

    # check detector calibrated Flag
    if __detector__.TileDetector.calibrated is True:
        if kwargs.get("overwrite"):
            warnings.warn("Warning: Overwriting old calibration data")
        else:
            warnings.warn("Warning: detector has already calibration data")
            return

    # Start timers
    t = Timer()
    t.start()

    # TODO: workaround
    resid_z, resid_phi = (None, None)

    # ------------------------------------------------------
    # reading data / preparing data
    histogram = read_histo(kwargs["hist_file"])

    print("Using ", kwargs["dt_mode"])
    if kwargs["dt_mode"] == "median":
        # calculate residuals to truth
        # and dt between tiles (median of the histogram)
        dt_z_rel, dt_phi_rel, resid_z, resid_phi = get_median_from_hist(histogram, resid=True)
    elif kwargs["dt_mode"] == "mean":
        # calculate dt between tiles (mean)
        # dt_z_rel, dt_phi_rel = get_mean_from_ttree(ttree_mu3e, threshold=kwargs["mean_threshold"])
        # TODO
        pass
    elif kwargs["dt_mode"] == "gaus":
        dt_z_rel, dt_phi_rel = get_mu_gaus(histogram)
    del histogram
    # ------------------------------------------------------

    # ------------------------------------------------------
    # correction relative dts
    # correction for phi loop (going around the loop should result in sum dt = 0)
    loop_correction_phi(__detector__, dt_phi_rel, 200000)
    loop_correction_phi(__detector__, dt_phi_rel, 300000)

    tof_correction_z(__detector__, dt_z_rel, 200000, kwargs["tof"])
    tof_correction_z(__detector__, dt_z_rel, 300000, kwargs["tof"])

    # loop_correction_z(__detector__, dt_phi_rel, dt_z_rel, 200000)
    # loop_correction_z(__detector__, dt_phi_rel, dt_z_rel, 300000)
    # ------------------------------------------------------

    # ------------------------------------------------------
    # calculating absolute values
    align_timings(dt_phi_rel, dt_z_rel, 200000)
    align_timings(dt_phi_rel, dt_z_rel, 300000)
    # ------------------------------------------------------

    # ------------------------------------------------------
    # finalizing calibration
    # set calibrated Flag in detector
    __detector__.TileDetector.calibrated = True

    print("Calibration finished")
    t.print()
    # ------------------------------------------------------
    #
    #
    #
    # ------------------------------------------------------
    # Only for testing below here:
    #
    # pre align in z direction
    cal_station_1_z, cal_station_2_z = pre_align_z(dt_z_rel.copy())
    #
    # pre align in phi direction
    cal_station_1_phi, cal_station_2_phi = pre_align_phi(dt_phi_rel.copy())
    #
    # return difference from truth
    cal = None
    if "debug_station" in kwargs.keys():
        if kwargs["debug_station"] == 1:
            cal = check_cal_z(cal_station_1_z, 1)
        elif kwargs["debug_station"] == 2:
            cal = check_cal_z(cal_station_2_z, 2)
    #
    cal_phi = None
    if "debug_station" in kwargs.keys():
        if kwargs["debug_station"] == 1:
            cal_phi = check_cal_phi(cal_station_1_phi, 1)
        elif kwargs["debug_station"] == 2:
            cal_phi = check_cal_phi(cal_station_2_phi, 2)
    # ------------------------------------------------------

    return resid_z, resid_phi, cal, cal_phi


# ---------------------------------------
def generate_hist(filename: str, **kwargs):
    # Start timers
    t = Timer()
    t.start()

    if __detector__ is None:
        raise RuntimeError("ERROR: Detector not selected")

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get(kwargs["ttree_loc"])
    histogram = fill_dt_histos(__detector__, ttree_mu3e, histo_options=kwargs["histo_options"])

    save_histo(kwargs["hist_file"], histogram)
    t.print()


# ---------------------------------------
def minimize_function_z(offset: float, dt_phi_rel: list, dt_z_rel: list, dt_phi_plus_one: list) -> float:
    tmp = 0.
    tmp += abs(offset - dt_z_rel[0])

    tmp_1 = 0.
    for i in range(len(dt_phi_rel) - 1):
        for j in range(i):
            tmp_1 += dt_phi_rel[j] - dt_phi_plus_one[j]

        tmp_1 += - offset + dt_z_rel[i + 1]

        tmp += abs(tmp_1)

    return tmp


# TODO
def align_timings(dt_phi_rel: dict, dt_z_rel: dict, station_offset: int):
    print(f"Calculating absolute timing offsets to master tile: {station_offset}")

    result = {}
    for z_position in range(len(__detector__.TileDetector.row_ids(0, station_offset)) - 1):
        # generating empty temporary arrays
        dt_phi_arr_tmp = []
        dt_phi_plus_one_arr_tmp = []
        dt_z_arr_tmp = []

        # generating column ids
        column_ids_for_current_z = __detector__.TileDetector.column_ids(z_position, station_offset)
        column_ids_for_current_z_plus_one = __detector__.TileDetector.column_ids(z_position + 1, station_offset)

        # filling arrays with relative dt information
        for tile_id in column_ids_for_current_z:
            dt_phi_arr_tmp.append(dt_phi_rel[tile_id] if tile_id in dt_phi_rel.keys() else 0)
            dt_z_arr_tmp.append(dt_z_rel[tile_id] if tile_id in dt_z_rel.keys() else 0)

        for tile_id in column_ids_for_current_z_plus_one:
            dt_phi_plus_one_arr_tmp.append(dt_phi_rel[tile_id] if tile_id in dt_phi_rel.keys() else 0)

        dt_phi_arr_tmp = np.asarray(dt_phi_arr_tmp)
        dt_phi_plus_one_arr_tmp = np.asarray(dt_phi_plus_one_arr_tmp)
        dt_z_arr_tmp = np.asarray(dt_z_arr_tmp)
        # optimizing z-direction
        res = minimize(minimize_function_z, x0=0., args=(dt_phi_arr_tmp, dt_z_arr_tmp, dt_phi_plus_one_arr_tmp),
                       tol=1e-6)

        # print(res.x[0])
        result[z_position] = res.x[0]

    # first row
    # we have to set the first tile as the master tile with offset = 0
    absolute_timing_offset = 0.
    for phi_id in __detector__.TileDetector.column_ids(0, station_offset):
        __detector__.TileDetector.tile[phi_id].dt_cal = absolute_timing_offset
        absolute_timing_offset += dt_phi_rel[phi_id if phi_id in dt_phi_rel.keys() else station_offset]

    # all other rows
    for z_position in range(1, len(__detector__.TileDetector.row_ids(0, station_offset))):
        last_master_id = int(__detector__.TileDetector.row_ids(0, station_offset)[z_position - 1])
        last_master_offset = __detector__.TileDetector.tile[last_master_id].dt_cal
        # tmp_id = __detector__.TileDetector.column_ids (0, station_offset)[z_position]
        current_master_offset = last_master_offset + result[z_position - 1]

        absolute_timing_offset = current_master_offset
        for phi_id in __detector__.TileDetector.column_ids(z_position, station_offset):
            # print(phi_id, dt_phi_rel[phi_id])
            __detector__.TileDetector.tile[phi_id].dt_cal = absolute_timing_offset
            absolute_timing_offset += dt_phi_rel[phi_id if phi_id in dt_phi_rel.keys() else station_offset]

    return


# ---------------------------------------

def pre_align_phi(dt_phi):
    warnings.warn("Warning: deprecated")
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
                    warnings.warn("Not enough data for phi calibration")

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
            dt_tmp = (__detector__.TileDetector.tile[tile_id].dt_truth -
                      __detector__.TileDetector.tile[__detector__.TileDetector.getNeighbour(tile_id, "up")].dt_truth)
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

def pre_align_z(dt_z):
    warnings.warn("Warning: deprecated")
    cal_station_1 = {}
    cal_station_2 = {}
    station_offset_arr = [200000, 300000]

    # calibration station 1 and 2
    for station_offset in station_offset_arr:
        for phi_row in range(56):
            # dt_tiles_cal[0] has to be deleted at the end
            dt_tiles_cal = [0.]
            for tile in range(0, 51):
                try:
                    dt_tmp = dt_z[(station_offset + phi_row + tile * 56)]
                except KeyError:
                    dt_tmp = 0
                    warnings.warn(f"Not enough data for z calibration")

                dt_tiles_cal.append(float(dt_tiles_cal[-1] + dt_tmp))

            if station_offset == 200000:
                cal_station_1[phi_row] = np.asarray(dt_tiles_cal[1:])
            elif station_offset == 300000:
                cal_station_2[phi_row] = np.asarray(dt_tiles_cal[1:])

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

    for phi_row in range(56):
        # dt_truth[0] has to be deleted at the end
        dt_truth = [0]
        for tile in range(0, 51):
            tile_id = station_offset + phi_row + tile * 56
            tile_id_neighbour = __detector__.TileDetector.getNeighbour(tile_id, "right")

            dt_tmp = (__detector__.TileDetector.tile[tile_id].dt_truth -
                      __detector__.TileDetector.tile[tile_id_neighbour].dt_truth)
            dt_truth.append(dt_truth[-1] + dt_tmp)

        cal[phi_row] = np.array(cal_data[phi_row]) + np.array(dt_truth[1:])

    return cal


# ---------------------------------------
# calculates the median from the histogram dict returned by fill_dt_histogram()
# returns residuals to true dt if resid=True

def get_median_from_hist(histogram: dict, resid=False):

    # z dir:
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
        dt_z[tile_entry] = q[0] + (
                    __detector__.TileDetector.tile[neighbour_z_id].dt_truth - __detector__.TileDetector.tile[
                tile_entry].dt_truth)
        if resid:
            resid_dt = q
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
        dt_phi[tile_entry] = q[0] + (
                    __detector__.TileDetector.tile[neighbour_phi_id].dt_truth - __detector__.TileDetector.tile[
                tile_entry].dt_truth)
        if resid:
            resid_dt = q
            resid_phi.append(resid_dt)

    return dt_z, dt_phi, resid_z, resid_phi


def get_mu_gaus(histogram: dict) -> (dict, dict):
    ROOT.TVirtualFitter.SetDefaultFitter("Minuit2")
    # z dir:
    dt_z = {}
    for tile_entry in histogram:
        neighbour_z_id = __detector__.TileDetector.getNeighbour(tile_entry, "right")
        if neighbour_z_id is False:
            continue

        if histogram[tile_entry][0].ComputeIntegral() == 0:
            # warnings.warn(f"WARNING: INTEGRAL = 0, tile: {tile_entry}")
            continue
        g1 = ROOT.TF1("fit_f", "gaus", -2, 2)
        h1 = histogram[tile_entry][0]
        h1.Rebin(4)
        ##g1.SetParameters(h1.GetMaximum(), h1.GetMean(), h1.GetRMS())
        ##g1.SetParLimits(0, 1, 1000)
        histogram[tile_entry][0].Fit(g1, "WMEGQR", "0")
        if g1.GetParameter(1) > 1:
            print("z", tile_entry, g1.GetParameter(0), g1.GetParameter(1), g1.GetParameter(2))
        dt_z[tile_entry] = g1.GetParameter(1) + (
                __detector__.TileDetector.tile[neighbour_z_id].dt_truth - __detector__.TileDetector.tile[
            tile_entry].dt_truth)
        del g1

    # phi dir:
    dt_phi = {}
    for tile_entry in histogram:

        neighbour_phi_id = __detector__.TileDetector.getNeighbour(tile_entry, "up")

        if histogram[tile_entry][1].ComputeIntegral() == 0:
            # print("WARNING: INTEGRAL = 0", tile_entry)
            continue
        g1 = ROOT.TF1("fit_f", "gaus", -2, 2)
        h1 = histogram[tile_entry][1]
        h1.Rebin(4)
        #g1.SetParameters(h1.GetMaximum(), h1.GetMean(), h1.GetRMS())
        #g1.SetParLimits(0, 1, 1000)
        histogram[tile_entry][1].Fit(g1, "WMEGQR", "0")
        if g1.GetParameter(1) > 1:
            print("phi", tile_entry, g1.GetParameter(1))
        dt_phi[tile_entry] = g1.GetParameter(1) + (
                __detector__.TileDetector.tile[neighbour_phi_id].dt_truth - __detector__.TileDetector.tile[
            tile_entry].dt_truth)

        del g1

    return dt_z, dt_phi


# ---------------------------------------
# DEPRECATED / UNUSED
# ---------------------------------------
def get_mean_from_ttree(ttree_mu3e, threshold: float) -> (dict, dict):
    dt_z = {}
    dt_phi = {}
    dt_z_n = {}
    dt_phi_n = {}

    cluster_counter = 0

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        # Printing status info
        if frame % 10000 == 0:
            print("Searching clusters. Progress: ", np.round(frame / ttree_mu3e.GetEntries() * 100), " % , Found: ",
                  cluster_counter, end='\r')

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
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[hit_tile].dt_truth
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[
                    neighbour_z_id].dt_truth
                dt = hit_time_2 - hit_time_1

                if abs(dt) <= threshold:
                    # Add to mean
                    if dt_z.get(hit_tile):
                        dt_z[hit_tile] = (dt_z[hit_tile] * dt_z_n[hit_tile] + dt) / (dt_z_n[hit_tile] + 1)
                        dt_z_n[hit_tile] += 1
                    else:
                        dt_z[hit_tile] = dt
                        dt_z_n[hit_tile] = 1

                    cluster_counter += 1

            # -----------------------------
            # Look for clusters in phi-dir
            neighbour_z_id = __detector__.TileDetector.getNeighbour(hit_tile, "up")
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
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[
                    hit_tile].dt_truth
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[
                    neighbour_z_id].dt_truth
                dt = hit_time_2 - hit_time_1

                if abs(dt) <= threshold:
                    # Add to mean
                    if dt_phi.get(hit_tile):
                        dt_phi[hit_tile] = (dt + dt_phi[hit_tile] * dt_phi_n[hit_tile]) / (dt_phi_n[hit_tile] + 1)
                        dt_phi_n[hit_tile] += 1
                    else:
                        dt_phi[hit_tile] = dt
                        dt_phi_n[hit_tile] = 1
                    cluster_counter += 1
    print("Searching clusters. Progress: ", 100, " % , Found: ", cluster_counter)
    return dt_z, dt_phi
