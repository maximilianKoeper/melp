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

    # Get pre align histogram
    histogram = fill_dt_histos(ttree_mu3e)

    # calculate residuals to truth
    # and dt between tiles

    # z dir:
    resid_z = []
    dt_z = {}
    for entry in histogram:
        prob = np.array([0.5])
        q = np.array([0.])
        if histogram[entry][0].ComputeIntegral() == 0:
            print("WARNING: INTEGRAL = 0")
            continue
        histogram[entry][0].GetQuantiles(1, q, prob)
        dt_z[entry] = q
        resid_dt = q - (__detector__.TileDetector.tile[entry + 56].dt - __detector__.TileDetector.tile[entry].dt)
        resid_z.append(resid_dt)

    # phi dir:
    resid_phi = []
    dt_phi = {}
    for entry in histogram:
        prob = np.array([0.5])
        q = np.array([0.])
        if histogram[entry][1].ComputeIntegral() == 0:
            print("WARNING: INTEGRAL = 0")
            continue
        histogram[entry][1].GetQuantiles(1, q, prob)
        dt_phi[entry] = q
        resid_dt = q - (__detector__.TileDetector.tile[entry + 1].dt - __detector__.TileDetector.tile[entry].dt)
        resid_phi.append(resid_dt)

    # calibrating along z and return difference from truth
    cal = {}
    for phi_row in range(56):
        dt_cal = [0]
        dt_truth = [0]
        for tile in range(0, 51):
            try:
                dt_tmp = dt_z[(200000 + phi_row + tile * 56)]
            except:
                dt_tmp = 0
                print("WARNING: Not enough data for calibration")

            dt_cal.append(dt_cal[-1] + dt_tmp)

            dt_tmp = (__detector__.TileDetector.tile[200000 + phi_row + tile * 56].dt -
                       __detector__.TileDetector.tile[200000 + phi_row + (tile+1) * 56].dt)
            dt_truth.append(dt_truth[-1] + dt_tmp)

        cal[phi_row] = np.array(dt_cal) + np.array(dt_truth)

    return resid_z, resid_phi, cal


# ---------------------------------------

def fill_dt_histos(ttree_mu3e) -> dict:
    hist = {}
    nbins, lo, hi = 1000, -100, 100
    # Generating empty histos:
    for tile in __detector__.TileDetector.tile:
        histo_name_z = str(tile) + "_z"
        histo_name_phi = str(tile) + "_phi"
        hist[tile] = [ROOT.TH1F(histo_name_z, 'my hist', nbins, lo, hi),
                      ROOT.TH1F(histo_name_phi, 'my hist', nbins, lo, hi)]

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
                hist[hit_tile][0].Fill(dt)

            # TODO: index_finder cant handle multiple events on one tile in one frame!!!
            #
            # Look for clusters in phi-dir
            if (hit_tile + 1) in ttree_mu3e.tilehit_tile:
                hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]
                # find associated tile hit
                hit_tile_assoc = index_finder(list(ttree_mu3e.tilehit_tile), (hit_tile + 56))

                try:
                    hit_tile_assoc = int(*hit_tile_assoc)
                except:
                    continue

                # calculate dt
                # TODO: TOF has to be added
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index] + __detector__.TileDetector.tile[hit_tile].dt
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc] + __detector__.TileDetector.tile[hit_tile + 1].dt
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist[hit_tile][1].Fill(dt)

    return hist
