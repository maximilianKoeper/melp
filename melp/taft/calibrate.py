import ROOT

import numpy as np

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

def calibrate(filename: str, **kwargs):
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

    # calibrating along z and return difference from truth
    cal = check_cal_z(dt_z)

    return resid_z, resid_phi, cal


# ---------------------------------------


def check_cal_z(dt_z):
    cal = {}
    for phi_row in range(55):
        dt_cal = [0]
        dt_truth = [0]
        for tile in range(0, 51):
            try:
                dt_tmp = dt_z[(200000 + phi_row + tile * 56)]
            except:
                dt_tmp = 0
                print("WARNING: Not enough data for calibration", tile)

            dt_cal.append(dt_cal[-1] + dt_tmp)

            dt_tmp = (__detector__.TileDetector.tile[200000 + phi_row + tile * 56].dt -
                      __detector__.TileDetector.tile[200000 + phi_row + (tile + 1) * 56].dt)
            dt_truth.append(dt_truth[-1] + dt_tmp)

        cal[phi_row] = np.array(dt_cal) + np.array(dt_truth)

    return cal


# ---------------------------------------


def get_median_from_hist(histogram: dict, resid=False):
    resid_z = []
    dt_z = {}
    for tile_entry in histogram:
        if (tile_entry + 56) not in __detector__.TileDetector.tile.keys():
            continue
        # setting prob to 0.5 for median
        prob = np.array([0.5])
        # q for result (has to be np.array because of type casting)
        q = np.array([0.])
        if histogram[tile_entry][0].ComputeIntegral() == 0:
            print("WARNING: INTEGRAL = 0", tile_entry)
            continue
        histogram[tile_entry][0].GetQuantiles(1, q, prob)
        dt_z[tile_entry] = q
        if resid:
            resid_dt = q - (__detector__.TileDetector.tile[tile_entry + 56].dt - __detector__.TileDetector.tile[tile_entry].dt)
            resid_z.append(resid_dt)

    # phi dir:
    resid_phi = []
    dt_phi = {}
    for tile_entry in histogram:

        # -----------------
        # checking for 56th tile (56 +1)th tile is not in the same ring anymore !!!
        next_tile_phi = 0
        tmp_hittile = tile_entry - 200000
        if tmp_hittile >= 100000:
            tmp_hittile -= 100000
        if (tmp_hittile + 1) % 56 == 0:
            next_tile_phi = tile_entry - 55
        else:
            next_tile_phi = tile_entry + 1
        # -----------------

        # setting prob to 0.5 for median
        prob = np.array([0.5])
        # q for result (has to be np.array because of type casting)
        q = np.array([0.])
        if histogram[tile_entry][1].ComputeIntegral() == 0:
            print("WARNING: INTEGRAL = 0", tile_entry)
            continue
        histogram[tile_entry][1].GetQuantiles(1, q, prob)
        dt_phi[tile_entry] = q
        if resid:
            resid_dt = q - (__detector__.TileDetector.tile[next_tile_phi].dt - __detector__.TileDetector.tile[tile_entry].dt)
            resid_phi.append(resid_dt)

    return dt_z, dt_phi, resid_z, resid_phi


# ---------------------------------------
#
# Generates dictionary with ROOT TH1 Histogramms
#   -> dict[tileid] = [hist_z, hist_pih]
#

def fill_dt_histos(ttree_mu3e) -> dict:
    hist_dict = {}
    nbins, lo, hi = 1000, -100, 100
    # Generating empty histos:
    for tile in __detector__.TileDetector.tile:
        histo_name_z = str(tile) + "_z"
        histo_name_phi = str(tile) + "_phi"
        hist_dict[tile] = [ROOT.TH1D(histo_name_z, 'my hist', nbins, lo, hi),
                      ROOT.TH1D(histo_name_phi, 'my hist', nbins, lo, hi)]

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        # TODO: index_finder cant handle multiple events on one tile in one frame!!!
        #
        # Analyzing frame
        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

            # Look for clusters in z-dir
            if (hit_tile + 56) in ttree_mu3e.tilehit_tile:
                # find associated tile hit
                hit_tile_assoc = index_finder(list(ttree_mu3e.tilehit_tile), (hit_tile + 56))

                try:
                    hit_tile_assoc = int(*hit_tile_assoc)
                except:
                    continue

                # calculate dt
                # TODO: TOF has to be added
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[hit_tile].dt
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[hit_tile + 56].dt
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist_dict[hit_tile][0].Fill(dt)

            # TODO: index_finder cant handle multiple events on one tile in one frame!!!
            #
            # Look for clusters in phi-dir

            # -----------------
            # checking for 56th tile (56 +1)th tile is not in the same ring anymore !!!
            next_tile_phi = 0
            tmp_hittile = hit_tile - 200000
            if tmp_hittile >= 100000:
                tmp_hittile -= 100000
            if (tmp_hittile+1) % 56 == 0:
                next_tile_phi = hit_tile - 55
            else:
                next_tile_phi = hit_tile + 1
            # -----------------

            if next_tile_phi in ttree_mu3e.tilehit_tile:
                hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]
                # find associated tile hit
                hit_tile_assoc = index_finder(list(ttree_mu3e.tilehit_tile), next_tile_phi)

                try:
                    hit_tile_assoc = int(*hit_tile_assoc)
                except:
                    continue

                # calculate dt
                # TODO: TOF has to be added
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[hit_tile].dt
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[next_tile_phi].dt
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist_dict[hit_tile][1].Fill(dt)

    return hist_dict
