import ROOT
import numpy as np

# fast index lookup
from melp.libs.misc import index_finder


def save_histo(filename: str, dt_dict: dict):
    histo_file = ROOT.TFile.Open(filename, "RECREATE")

    for keys in dt_dict.keys():
        name_z = str(keys) + "z"
        name_phi = str(keys) + "phi"
        histo_file.WriteObject(dt_dict[keys][0], name_z)
        histo_file.WriteObject(dt_dict[keys][1], name_phi)


def read_histo(filename: str) -> dict:
    global histo_file
    histo_file = ROOT.TFile.Open(filename, "READ")

    dt_dict = {}

    for key in histo_file.GetListOfKeys():
        h = key.ReadObj()
        name = h.GetName()
        dict_key = name.replace("_z", "")
        dict_key = int(dict_key.replace("_phi", ""))

        if dict_key not in dt_dict.keys():
            dt_dict[dict_key] = [None, None]

        if "z" in name:
            dt_dict[dict_key][0] = h
            # print(h)
        elif "phi" in name:
            dt_dict[dict_key][1] = h

    return dt_dict


# ---------------------------------------
#
# Generates dictionary with ROOT TH1D Histogramms
#   -> dict[tileid] = [hist_z, hist_pih]
#

def fill_dt_histos(detector, ttree_mu3e, histo_options: tuple) -> dict:
    cluster_counter = 0

    hist_dict = {}
    nbins, lo, hi = histo_options
    # Generating empty histos:
    for tile in detector.TileDetector.tile:
        histo_name_z = str(tile) + "_z"
        histo_name_phi = str(tile) + "_phi"
        hist_dict[tile] = [ROOT.TH1D(histo_name_z, histo_name_z, nbins, lo, hi),
                           ROOT.TH1D(histo_name_phi, histo_name_phi, nbins, lo, hi)]

    # tilehits = ROOT.vector('int')()
    # tilehitstime = ROOT.vector('double')()
    # ttree_mu3e.SetBranchStatus("tilehit_tile", 1)
    # ttree_mu3e.SetBranchStatus("tilehit_time", 1)

    # ttree_mu3e.SetBranchAddress("tilehit_tile", tilehits)
    # ttree_mu3e.SetBranchAddress("tilehit_time", tilehitstime)

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        # Printing status info
        if frame % 10000 == 0:
            print("Searching clusters. Progress: ", np.round(frame / ttree_mu3e.GetEntries() * 100), " % , Found: ",
                  cluster_counter, end='\r')

        # TODO: index_finder cant handle multiple events on one tile in one frame!!!
        #       --> skipping frame (looses some data)
        # Analyzing frame
        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

            # -----------------------------
            # Look for clusters in z-dir
            neighbour_z_id = detector.TileDetector.getNeighbour(hit_tile, "right")
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
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index]  # + detector.TileDetector.tile[hit_tile].dt_truth
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc]  # + detector.TileDetector.tile[
                # neighbour_z_id].dt_truth
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist_dict[hit_tile][0].Fill(dt)
                cluster_counter += 1

            # -----------------------------
            # Look for clusters in phi-dir
            neighbour_phi_id = detector.TileDetector.getNeighbour(hit_tile, "up")
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
                hit_time_1 = ttree_mu3e.tilehit_time[hit_tile_index]  # + detector.TileDetector.tile[hit_tile].dt_truth
                hit_time_2 = ttree_mu3e.tilehit_time[hit_tile_assoc]  # + detector.TileDetector.tile[
                # neighbour_phi_id].dt_truth
                dt = hit_time_2 - hit_time_1

                # Fill histogram
                hist_dict[hit_tile][1].Fill(dt)
                cluster_counter += 1

    print("Searching clusters. Progress: ", 100, " % , Found: ", cluster_counter)
    return hist_dict
