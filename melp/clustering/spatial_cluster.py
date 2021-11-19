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

#------------------------------------------
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
            tile_left_top = tile_centre - 57
            tile_left_centre = tile_centre - 56
            tile_left_bottom = tile_centre - 55
            tile_right_top = tile_centre + 55
            tile_right_centre = tile_centre + 56
            tile_right_bootom = tile_centre + 57
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

    return mask

#-----------------------------------------------
def build_mask_detector_class(filename, frame, mu3e_detector: melp.Detector, mask_type):
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    mask = {}
    ttree_mu3e.GetEntry(frame)
    if mask_type == "small":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = neighbour_wiper(mu3e_detector, tile_centre, "up")
            tile_centre_bottom = neighbour_wiper(mu3e_detector, tile_centre, "down")
            tile_left_centre = neighbour_wiper(mu3e_detector, tile_centre, "left")
            tile_right_centre = neighbour_wiper(mu3e_detector, tile_centre, "right")
            mask_tmp.append(tile_centre)
            mask_tmp.append(tile_centre_top)
            mask_tmp.append(tile_centre_bottom)
            mask_tmp.append(tile_left_centre)
            mask_tmp.append(tile_right_centre)

            if -1 in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != -1]

            mask[tile_centre] = mask_tmp

    if mask_type == "medium":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = int(neighbour_wiper(mu3e_detector, tile_centre, "up"))
            tile_centre_bottom = int(neighbour_wiper(mu3e_detector, tile_centre, "down"))
            tile_left_centre = int(neighbour_wiper(mu3e_detector, tile_centre, "left"))
            tile_left_top =  int(neighbour_wiper(mu3e_detector, tile_left_centre, "up"))
            tile_left_bottom = int(neighbour_wiper(mu3e_detector, tile_left_centre, "down"))
            tile_right_centre = int(neighbour_wiper(mu3e_detector, tile_centre, "right"))
            tile_right_top = int(neighbour_wiper(mu3e_detector, tile_right_centre, "up"))
            tile_right_bootom = int(neighbour_wiper(mu3e_detector, tile_right_centre, "down"))
            mask_tmp.append(tile_centre)
            mask_tmp.append(tile_centre_top)
            mask_tmp.append(tile_centre_bottom)
            mask_tmp.append(tile_left_top)
            mask_tmp.append(tile_left_centre)
            mask_tmp.append(tile_left_bottom)
            mask_tmp.append(tile_right_top)
            mask_tmp.append(tile_right_centre)
            mask_tmp.append(tile_right_bootom)
            if -1 in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != -1]
            mask[tile_centre] = mask_tmp

    if mask_type == "big":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = neighbour_wiper(mu3e_detector, tile_centre, "up")
            tile_centre_bottom = neighbour_wiper(mu3e_detector, tile_centre, "down")
            tile_left_centre = neighbour_wiper(mu3e_detector, tile_centre, "left")
            tile_left_top =  neighbour_wiper(mu3e_detector, tile_left_centre, "up")
            tile_left_bottom = neighbour_wiper(mu3e_detector, tile_left_centre, "down")
            tile_right_centre = neighbour_wiper(mu3e_detector, tile_centre, "right")
            tile_right_top = neighbour_wiper(mu3e_detector, tile_right_centre, "up")
            tile_right_bottom = neighbour_wiper(mu3e_detector, tile_right_centre, "down")

            tile_centre_far_top = neighbour_wiper(mu3e_detector, tile_centre_top, "up")
            tile_centre_far_bottom = neighbour_wiper(mu3e_detector, tile_centre_bottom, "down")
            tile_left_far_centre = neighbour_wiper(mu3e_detector, tile_left_centre, "left")
            tile_left_far_top =  neighbour_wiper(mu3e_detector, tile_left_top, "up")
            tile_left_far_bottom = neighbour_wiper(mu3e_detector, tile_left_bottom, "down")
            tile_right_far_centre = neighbour_wiper(mu3e_detector, tile_right_centre, "right")
            tile_right_far_top = neighbour_wiper(mu3e_detector, tile_right_top, "up")
            tile_right_far_bottom = neighbour_wiper(mu3e_detector, tile_right_bottom, "down")

            tile_far_left_top = neighbour_wiper(mu3e_detector, tile_left_top, "left")
            tile_far_left_bottom = neighbour_wiper(mu3e_detector, tile_left_bottom, "left")
            tile_far_right_top = neighbour_wiper(mu3e_detector, tile_right_top, "right")
            tile_far_right_bottom = neighbour_wiper(mu3e_detector, tile_right_bottom, "right")

            mask_tmp.append(tile_centre)
            mask_tmp.append(tile_centre_top)
            mask_tmp.append(tile_centre_bottom)
            mask_tmp.append(tile_left_top)
            mask_tmp.append(tile_left_centre)
            mask_tmp.append(tile_left_bottom)
            mask_tmp.append(tile_right_top)
            mask_tmp.append(tile_right_centre)
            mask_tmp.append(tile_right_bottom)

            mask_tmp.append(tile_centre_far_top)
            mask_tmp.append(tile_centre_far_bottom)
            mask_tmp.append(tile_left_far_centre)
            mask_tmp.append(tile_left_far_top)
            mask_tmp.append(tile_left_far_bottom)
            mask_tmp.append(tile_right_far_centre)
            mask_tmp.append(tile_right_far_top)
            mask_tmp.append(tile_right_far_bottom)

            mask_tmp.append(tile_far_left_top)
            mask_tmp.append(tile_far_left_bottom)
            mask_tmp.append(tile_far_right_top)
            mask_tmp.append(tile_far_right_bottom)

            if -1 in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != -1]

            mask[tile_centre] = mask_tmp

    return mask

#-------------------------------------------
def build_cluster_mask(filename, frame):
    pass