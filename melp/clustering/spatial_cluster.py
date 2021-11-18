import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

"""
def spatial_truth_clusters_frame(filename,frame,threshold):
    clusters = {}
    cluster_counter = 0
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")

    return clusters
"""

def build_mask(filename, frame, mask_type = "medium"): #build mask around hit
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    mask = {}
    ttree_mu3e.GetEntry(frame)
    if mask_type == "medium":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = tile_centre - 1
            tile_centre_bottom = tile_centre + 1
            tile_left_top = tile_centre - 53
            tile_left_centre = tile_centre - 52
            tile_left_bottom = tile_centre - 51
            tile_right_top = tile_centre + 51
            tile_right_centre = tile_centre + 52
            tile_right_bootom = tile_centre + 53
            mask_tmp.append(tile_centre)
            mask_tmp.append(tile_centre_top)
            mask_tmp.append(tile_centre_bottom)
            mask_tmp.append(tile_left_top)
            mask_tmp.append(tile_left_centre)
            mask_tmp.append(tile_left_bottom)
            mask_tmp.append(tile_right_top)
            mask_tmp.append(tile_right_centre)
            mask_tmp.append(tile_right_bootom)
            mask[tile_centre] = mask_tmp

    else:
        print("Fail")

    return mask




