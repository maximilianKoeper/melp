import warnings

import ROOT
import numpy as np
from scipy.optimize import minimize

from melp import Detector
# fast index lookup
from melp.libs.timer import Timer
from melp.taft.corrections.misc_corrections import loop_correction_phi, loop_adv_correction_phi, loop_adv_correction_z
# different functions for calibration
from melp.taft.corrections.tof_corrections import tof_correction_z
from melp.taft.corrections.global_fit import correct_z_two_event
from melp.taft.utils.root_helper import save_histo, read_histo, fill_dt_histos

# ---------------------------------------------------------------------
#  Define global variables and functions to select Detector
# ---------------------------------------------------------------------

__detector__: Detector = None


def select(selection: Detector):
    global __detector__
    __detector__ = selection


# ---------------------------------------------------------------------
# calibrate functions calls each step for the calibration and writes result
# to the detector object
#
# after the calibration the calibrated Flag in the detector is set to True
# ---------------------------------------------------------------------

def calibrate(**kwargs):
    # TODO: tidy up and streamline options
    global __detector__

    # ------------------------------------------------------
    # checking if everything is set up correctly
    if __detector__ is None:
        raise RuntimeError("ERROR: Detector not selected")

    # check detector calibrated Flag
    if __detector__.TileDetector.calibrated is True:
        if kwargs.get("overwrite"):
            warnings.warn("Warning: Overwriting old calibration data")
        else:
            warnings.warn("Warning: detector has already calibration data")
            return

    # ------------------------------------------------------
    # Start timers
    t = Timer()
    t.start()

    # TODO: workaround (not all functions return residuals)
    resid_z, resid_phi = (None, None)

    # ------------------------------------------------------
    # reading data / preparing data
    histogram = read_histo(kwargs["hist_file"])

    print("*Using dt mode: ", kwargs["dt_mode"])
    if kwargs["dt_mode"] == "median":
        # calculate residuals to truth
        # and dt between tiles (median of the histogram)
        dt_z_rel, dt_phi_rel, resid_z, resid_phi = get_median_from_hist(histogram, resid=True)
    elif kwargs["dt_mode"] == "mean":
        # calculate dt between tiles (mean)
        dt_z_rel, dt_phi_rel = get_mean_from_hist(histogram)
    elif kwargs["dt_mode"] == "gaus":
        dt_z_rel, dt_phi_rel = get_mu_gaus(histogram)
    del histogram
    # ------------------------------------------------------

    # ------------------------------------------------------
    # correction relative dts
    # correction for phi loop (going around the loop should result in sum dt = 0)
    loop_correction_phi(__detector__, dt_phi_rel, 200000)
    loop_correction_phi(__detector__, dt_phi_rel, 300000)

    # ------------------------------------------------------
    # correction relative dts
    # correction for phi loops
    if kwargs.get("phi_penalties") is not None:
        loop_adv_correction_phi(__detector__, dt_phi_rel, dt_z_rel.copy(), 200000, penalties=kwargs["phi_penalties"])
        loop_adv_correction_phi(__detector__, dt_phi_rel, dt_z_rel.copy(), 300000, penalties=kwargs["phi_penalties"])

    # ------------------------------------------------------
    # correction relative dts
    # correction for z
    if kwargs.get("z_penalties") is not None:
        loop_adv_correction_z(__detector__, dt_phi_rel.copy(), dt_z_rel, 200000, penalties=kwargs["z_penalties"])
        loop_adv_correction_z(__detector__, dt_phi_rel.copy(), dt_z_rel, 300000, penalties=kwargs["z_penalties"])

    # ------------------------------------------------------
    # correcting for tof in z direction
    tof_correction_z(__detector__, dt_z_rel, station_offset=200000, tof_mode=kwargs["tof"])
    tof_correction_z(__detector__, dt_z_rel, station_offset=300000, tof_mode=kwargs["tof"])
    # ------------------------------------------------------

    # ------------------------------------------------------
    # calculating absolute values
    align_timings(dt_phi_rel, dt_z_rel, 200000)
    align_timings(dt_phi_rel, dt_z_rel, 300000)
    # ------------------------------------------------------
    #correct_z_error_prop(dt_z_rel)
    # ------------------------------------------------------
    # correcting tof with cosmic data (z - direction)
    # cosmic_correction_z(__detector__, **kwargs)
    popt_1, popt_2 = (None, None)
    if kwargs.get("cosmic_correction"):
        popt_1 = correct_z_two_event(__detector__, cosmic_station=1, **kwargs)
        popt_2 = correct_z_two_event(__detector__, cosmic_station=2, **kwargs)
    # ------------------------------------------------------

    # ------------------------------------------------------
    # finalizing calibration
    # set calibrated Flag in detector
    __detector__.TileDetector.calibrated = True

    # ------------------------------------------------------
    # deleting cached function calls
    __detector__.TileDetector.getNeighbour.cache_clear()
    __detector__.TileDetector.row_ids.cache_clear()
    __detector__.TileDetector.column_ids.cache_clear()
    __detector__.TileDetector.id_from_row_col.cache_clear()

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

    return resid_z, resid_phi, cal, cal_phi, popt_1, popt_2


# ------------------------------------------------------
# generates Histogram file
# if possible use ROOT Macro (melp/taft/cpp)
def generate_hist(filename: str, **kwargs):
    global __detector__
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


# ------------------------------------------------------
# currently used function to minimize error in z direction
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
# prepares data for minimizer and applies result to detector
def align_timings(dt_phi_rel: dict, dt_z_rel: dict, station_offset: int):
    global __detector__
    print(f"*Calculating absolute timing offsets to master tile: {station_offset}")

    result = {}
    # loop over all z - columns
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
            dt_phi_arr_tmp.append(dt_phi_rel[tile_id])
            dt_z_arr_tmp.append(dt_z_rel[tile_id])

        for tile_id in column_ids_for_current_z_plus_one:
            dt_phi_plus_one_arr_tmp.append(dt_phi_rel[tile_id])

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
        # ATTENTION: directly accessing timing data (can lead to bugs)
        __detector__.TileDetector.tile[phi_id].dt_cal = absolute_timing_offset
        absolute_timing_offset += dt_phi_rel[phi_id]

    # all remaining rows
    for z_position in range(1, len(__detector__.TileDetector.row_ids(0, station_offset))):
        last_master_id = int(__detector__.TileDetector.row_ids(0, station_offset)[z_position - 1])
        last_master_offset = __detector__.TileDetector.tile[last_master_id].get_calibrated_offset()
        # tmp_id = __detector__.TileDetector.column_ids (0, station_offset)[z_position]
        current_master_offset = last_master_offset + result[z_position - 1]

        absolute_timing_offset = current_master_offset
        for phi_id in __detector__.TileDetector.column_ids(z_position, station_offset):
            # print(phi_id, dt_phi_rel[phi_id])
            __detector__.TileDetector.tile[phi_id].dt_cal = absolute_timing_offset
            absolute_timing_offset += dt_phi_rel[phi_id]

    return


# ------------------------------------------------------
# calculates the median from the histogram dict returned by fill_dt_histogram()
# returns residuals to true dt if resid=True

def get_median_from_hist(histogram: dict, resid=False):
    global __detector__
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
                __detector__.TileDetector.tile[neighbour_z_id].get_offset() - __detector__.TileDetector.tile[
            tile_entry].get_offset())
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
                __detector__.TileDetector.tile[neighbour_phi_id].get_offset() - __detector__.TileDetector.tile[
            tile_entry].get_offset())
        if resid:
            resid_dt = q
            resid_phi.append(resid_dt)

    return dt_z, dt_phi, resid_z, resid_phi


# ------------------------------------------------------
# fits a gaussian to every histogram to get mu
def get_mu_gaus(histogram: dict) -> (dict, dict):
    global __detector__
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
        # g1.SetParameters(h1.GetMaximum(), h1.GetMean(), h1.GetRMS())
        # g1.SetParLimits(0, 1, 1000)
        histogram[tile_entry][0].Fit(g1, "WMEGQR", "0")
        # if g1.GetParameter(1) > 1:
        #    print("z", tile_entry, g1.GetParameter(0), g1.GetParameter(1), g1.GetParameter(2))
        dt_z[tile_entry] = g1.GetParameter(1) + (
                __detector__.TileDetector.tile[neighbour_z_id].get_offset() - __detector__.TileDetector.tile[
            tile_entry].get_offset())
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
        # g1.SetParameters(h1.GetMaximum(), h1.GetMean(), h1.GetRMS())
        # g1.SetParLimits(0, 1, 1000)
        histogram[tile_entry][1].Fit(g1, "WMEGQR", "0")
        # if g1.GetParameter(1) > 1:
        #    print("phi", tile_entry, g1.GetParameter(1))
        dt_phi[tile_entry] = g1.GetParameter(1) + (
                __detector__.TileDetector.tile[neighbour_phi_id].get_offset() - __detector__.TileDetector.tile[
            tile_entry].get_offset())

        del g1

    return dt_z, dt_phi


# ------------------------------------------------------
# gets the mean from each histogram
def get_mean_from_hist(histogram: dict) -> (dict, dict):
    global __detector__
    # z dir:
    dt_z = {}
    for tile_entry in histogram:
        neighbour_z_id = __detector__.TileDetector.getNeighbour(tile_entry, "right")
        if neighbour_z_id is False:
            continue

        if histogram[tile_entry][0].ComputeIntegral() == 0:
            # warnings.warn(f"WARNING: INTEGRAL = 0, tile: {tile_entry}")
            continue

        mean = histogram[tile_entry][0].GetMean()

        dt_z[tile_entry] = mean + (
                __detector__.TileDetector.tile[neighbour_z_id].get_offset() - __detector__.TileDetector.tile[
            tile_entry].get_offset())

    # phi dir:
    dt_phi = {}
    for tile_entry in histogram:

        neighbour_phi_id = __detector__.TileDetector.getNeighbour(tile_entry, "up")

        if histogram[tile_entry][1].ComputeIntegral() == 0:
            # print("WARNING: INTEGRAL = 0", tile_entry)
            continue

        mean = histogram[tile_entry][1].GetMean()

        dt_phi[tile_entry] = mean + (
                __detector__.TileDetector.tile[neighbour_phi_id].get_offset() - __detector__.TileDetector.tile[
            tile_entry].get_offset())

    return dt_z, dt_phi


# --------------------------------------------------------
# ---------------------  DEPRECATED  ---------------------
# --------------------------------------------------------

# ------------------------------------------------------

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


# ------------------------------------------------------

def check_cal_phi(cal_data, station):
    global __detector__
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


# ------------------------------------------------------
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


# ------------------------------------------------------

def check_cal_z(cal_data, station):
    global __detector__
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


# --------------------------------------------------------
# ---------------------    UNUSED    ---------------------
# --------------------------------------------------------


def correct_z_error_prop(dt_z):
    global __detector__

    Tile = __detector__.TileDetector.tile
    time_diff = Tile[202856].dt_cal - Tile[200000].dt_cal

    Tile_ids_row = __detector__.TileDetector.row_ids(0, 200000)
    tmp_time_diff = 0
    for ids in Tile_ids_row:
        try:
            tmp_time_diff += dt_z[ids]
        except:
            pass

    time_diff_to_correct = time_diff - tmp_time_diff
    time_diff_to_correct /= 52

    for z_pos in range(52):
        for phi_pos in __detector__.TileDetector.column_ids(z_pos, 200000):
            __detector__.TileDetector.tile[phi_pos].dt_cal -= time_diff_to_correct * z_pos

    print("Correcting accumulated z-error: ", np.round(time_diff, 4), " >> ", np.round(tmp_time_diff, 4))

    time_diff = Tile[302857].dt_cal - Tile[300001].dt_cal

    Tile_ids_row = __detector__.TileDetector.row_ids(1, 300000)
    tmp_time_diff = 0
    for ids in Tile_ids_row:
        try:
            tmp_time_diff += dt_z[ids]
        except:
            pass

    time_diff_to_correct = time_diff - tmp_time_diff
    time_diff_to_correct /= 52

    for z_pos in range(52):
        for phi_pos in __detector__.TileDetector.column_ids(z_pos, 300000):
            __detector__.TileDetector.tile[phi_pos].dt_cal -= time_diff_to_correct * z_pos

    print("Correcting accumulated z-error: ", np.round(time_diff,4), " >> ", np.round(tmp_time_diff,4))
