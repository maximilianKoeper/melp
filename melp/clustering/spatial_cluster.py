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
#build mask around every hit in frame using tile ids (outdated)
def build_mask(filename, frame, mask_type = "medium"): 
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
#builds masks around every hit in frame using the get neighbour function, which works if hit is at the edge of the detector.
def build_mask_detector_class(filename, frame, mu3e_detector: melp.Detector, mask_type):
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    mask = {}
    ttree_mu3e.GetEntry(frame)
    if mask_type == "small":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = mu3e_detector.TileDetector.getNeighbour(tile_centre, "up")
            tile_centre_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre, "down")
            tile_left_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "left")
            tile_right_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "right")
            mask_tmp.append(tile_centre)
            mask_tmp.append(tile_centre_top)
            mask_tmp.append(tile_centre_bottom)
            mask_tmp.append(tile_left_centre)
            mask_tmp.append(tile_right_centre)

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    if mask_type == "medium":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = mu3e_detector.TileDetector.getNeighbour(tile_centre, "up")
            tile_centre_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre, "down")
            tile_left_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "left")
            tile_left_top = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "up")
            tile_left_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "down")
            tile_right_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "right")
            tile_right_top = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "up")
            tile_right_bootom = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "down")
            mask_tmp.append(tile_centre)
            mask_tmp.append(tile_centre_top)
            mask_tmp.append(tile_centre_bottom)
            mask_tmp.append(tile_left_top)
            mask_tmp.append(tile_left_centre)
            mask_tmp.append(tile_left_bottom)
            mask_tmp.append(tile_right_top)
            mask_tmp.append(tile_right_centre)
            mask_tmp.append(tile_right_bootom)

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    if mask_type == "big":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            mask_tmp = []
            tile_centre = ttree_mu3e.tilehit_tile[i]
            tile_centre_top = mu3e_detector.TileDetector.getNeighbour(tile_centre, "up")
            tile_centre_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre, "down")
            tile_left_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "left")
            tile_left_top =  mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "up")
            tile_left_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "down")
            tile_right_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "right")
            tile_right_top = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "up")
            tile_right_bottom = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "down")

            tile_centre_far_top = mu3e_detector.TileDetector.getNeighbour(tile_centre_top, "up")
            tile_centre_far_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre_bottom, "down")
            tile_left_far_centre = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "left")
            tile_left_far_top =  mu3e_detector.TileDetector.getNeighbour(tile_left_top, "up")
            tile_left_far_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_bottom, "down")
            tile_right_far_centre = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "right")
            tile_right_far_top = mu3e_detector.TileDetector.getNeighbour(tile_right_top, "up")
            tile_right_far_bottom = mu3e_detector.TileDetector.getNeighbour(tile_right_bottom, "down")

            tile_far_left_top = mu3e_detector.TileDetector.getNeighbour(tile_left_top, "left")
            tile_far_left_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_bottom, "left")
            tile_far_right_top = mu3e_detector.TileDetector.getNeighbour(tile_right_top, "right")
            tile_far_right_bottom = mu3e_detector.TileDetector.getNeighbour(tile_right_bottom, "right")

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

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    return mask

#------------------------------------------------
#builds masks around every hit in frame with hid=1,-1 using the get neighbour function, which works if hit is at the edge of the detector.
def build_mask_around_cluster_primary(filename, frame, mu3e_detector: melp.Detector, mask_type):
    primaries = get_cluster_primary_truth_frame(filename, frame)

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    mask = {}
    ttree_mu3e.GetEntry(frame)
    if mask_type == "small":
        for i in range(len(primaries)):
            tile_centre = primaries[i]
            tile_centre_top = mu3e_detector.TileDetector.getNeighbour(tile_centre, "up")
            tile_centre_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre, "down")
            tile_left_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "left")
            tile_right_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "right")

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, 
                                tile_left_centre, tile_right_centre])

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    if mask_type == "medium":
        for i in range(len(primaries)):
            tile_centre = primaries[i]
            tile_centre_top = mu3e_detector.TileDetector.getNeighbour(tile_centre, "up")
            tile_centre_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre, "down")
            tile_left_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "left")
            tile_left_top = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "up")
            tile_left_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "down")
            tile_right_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "right")
            tile_right_top = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "up")
            tile_right_bootom = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "down")

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, 
                                tile_left_centre, tile_left_bottom, tile_right_top, tile_right_centre, 
                                tile_right_bootom])

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    if mask_type == "big":
        for i in range(len(primaries)):
            tile_centre = primaries[i]
            tile_centre_top = mu3e_detector.TileDetector.getNeighbour(tile_centre, "up")
            tile_centre_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre, "down")
            tile_left_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "left")
            tile_left_top =  mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "up")
            tile_left_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "down")
            tile_right_centre = mu3e_detector.TileDetector.getNeighbour(tile_centre, "right")
            tile_right_top = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "up")
            tile_right_bottom = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "down")

            tile_centre_far_top = mu3e_detector.TileDetector.getNeighbour(tile_centre_top, "up")
            tile_centre_far_bottom = mu3e_detector.TileDetector.getNeighbour(tile_centre_bottom, "down")
            tile_left_far_centre = mu3e_detector.TileDetector.getNeighbour(tile_left_centre, "left")
            tile_left_far_top =  mu3e_detector.TileDetector.getNeighbour(tile_left_top, "up")
            tile_left_far_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_bottom, "down")
            tile_right_far_centre = mu3e_detector.TileDetector.getNeighbour(tile_right_centre, "right")
            tile_right_far_top = mu3e_detector.TileDetector.getNeighbour(tile_right_top, "up")
            tile_right_far_bottom = mu3e_detector.TileDetector.getNeighbour(tile_right_bottom, "down")

            tile_far_left_top = mu3e_detector.TileDetector.getNeighbour(tile_left_top, "left")
            tile_far_left_bottom = mu3e_detector.TileDetector.getNeighbour(tile_left_bottom, "left")
            tile_far_right_top = mu3e_detector.TileDetector.getNeighbour(tile_right_top, "right")
            tile_far_right_bottom = mu3e_detector.TileDetector.getNeighbour(tile_right_bottom, "right")

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, tile_left_centre, tile_left_bottom,
                                   tile_right_top, tile_right_centre, tile_right_bottom, tile_centre_far_top, tile_centre_far_bottom, 
                                   tile_left_far_centre, tile_left_far_top, tile_left_far_top, tile_left_far_bottom, tile_right_far_centre, 
                                   tile_right_far_top, tile_right_far_bottom, tile_far_left_top, tile_far_left_bottom, tile_far_right_top, 
                                   tile_far_right_bottom])

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    return mask 

#----------------------------------------------------
#builds clusters in the masks around hit with hid=1,-1 according to primaries.
def build_clusters_in_masks(filename, frame, mu3e_detector: melp.Detector, mask_type):
    primary_masks = build_mask_around_cluster_primary(filename, frame, mu3e_detector, mask_type)

    clusters = {} #keys: "master"-tile; values: rest of cluster

    keys = []
    values = []
    for key in primary_masks.keys():
        keys.append(key)
        values.append(primary_masks[key])

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e.GetEntry(frame)

    hit_tiles = []
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])


    for tile_id in hit_tiles:
        if tile_id < 300000:
            cluster_tmp = []
            cluster_primary_tmp = 0
            if tile_id not in primary_masks.keys(): #if not primary
                for i in range(len(values)):
                    if tile_id in values[i]:
                        cluster_tmp.append(tile_id)
                        cluster_primary_tmp = keys[i]



            if cluster_primary_tmp != 0:
                if cluster_primary_tmp not in clusters.keys() and len(cluster_tmp) > 0:              
                    clusters[cluster_primary_tmp] = cluster_tmp
                elif len(cluster_tmp) > 0:
                    clusters[cluster_primary_tmp].append(cluster_tmp[0])

    for i in keys:
        if i not in clusters.keys():
            clusters[i] = []

    return clusters

#-----------------------------------------------------
#builds clusters where dict-key is the primary of "master"-tile and not "master" tile and value is the whole cluster
def build_cluster_with_truth_primary(filename, frame, mu3e_detector: melp.Detector, mask_type):
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e.GetEntry(frame)

    #get clusters
    clusters_frame = build_clusters_in_masks(filename, frame, mu3e_detector, mask_type)
    whole_clusters = []
    clusters_with_primaries = {} #returned: keys:primary of "master"-tile; values:whole
    for key in clusters_frame.keys():
        whole_clusters_tmp = []
        whole_clusters_tmp.append(key)
        for i in clusters_frame[key]:
            whole_clusters_tmp.append(i)
        whole_clusters.append(np.array(whole_clusters_tmp))

        tilehit_primary = []
        for i in range(len(ttree_mu3e.tilehit_tile)):
            tile = ttree_mu3e.tilehit_tile[i]
            if tile == key:
                primary = ttree_mu3e.tilehit_primary[i]
                tilehit_primary.append(primary)

        for j in tilehit_primary:
            clusters_with_primaries[j] = whole_clusters_tmp
    
    return clusters_with_primaries
