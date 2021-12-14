import ROOT
import numpy as np
import melp
from melp import Detector
from melp.clustering.misc import*




#-----------------------------------------------
#builds masks around every hit in frame using the get neighbour function, which works if hit is at the edge of the detector.
def build_mask_detector_class(ttree_mu3e, mu3e_detector: melp.Detector, mask_type):
    mask = {}

    if mask_type == "small":
        for i in range(len(ttree_mu3e.tilehit_tile)):
            tile_centre = ttree_mu3e.tilehit_tile[i]
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
        for i in range(len(ttree_mu3e.tilehit_tile)):
            tile_centre = ttree_mu3e.tilehit_tile[i]
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
        for i in range(len(ttree_mu3e.tilehit_tile)):
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

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, tile_left_centre, tile_left_bottom,
                                            tile_right_top, tile_right_centre, tile_right_bottom, tile_centre_far_top, tile_centre_far_bottom, 
                                            tile_left_far_centre, tile_left_far_top, tile_left_far_bottom, tile_right_far_centre, 
                                            tile_right_far_top, tile_right_far_bottom, tile_far_left_top, tile_far_left_bottom, tile_far_right_top, 
                                            tile_far_right_bottom])

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    return mask

#------------------------------------------------
#builds masks around every hit in frame with hid=1,-1 using the get neighbour function, which also works if hit is at the edge of the detector.
def build_mask_around_cluster_master(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    #select reconstruction/tracking method
    if rec_type == "pixelpixel":
        mask_masters = melp.clustering.tracking.get_mask_masters_hitAnglePixelRec_without_hid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, matching="nearest")
    elif rec_type == "pixelpixelcheck":
        mask_masters = melp.clustering.tracking.get_mask_masters_hitAnglePixelRec_with_hid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, matching="nearest")
    else:
        mask_masters = melp.clustering.spatial_cluster.get_cluster_primary_truth_frame(ttree_mu3e, ttree_mu3e_mc, frame)

    mask = {}
    if mask_type == "small":
        for i in range(len(mask_masters)):
            tile_centre = mask_masters[i]
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
        for i in range(len(mask_masters)):
            tile_centre = mask_masters[i]
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
        for i in range(len(mask_masters)):
            tile_centre = mask_masters[i]
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
                                   tile_left_far_centre, tile_left_far_top, tile_left_far_bottom, tile_right_far_centre, 
                                   tile_right_far_top, tile_right_far_bottom, tile_far_left_top, tile_far_left_bottom, tile_far_right_top, 
                                   tile_far_right_bottom])

            if False in mask_tmp:
                mask_tmp = [x for x in mask_tmp if x != False]

            mask[tile_centre] = mask_tmp

    return mask 

#------------------------------------------------
#builds masks around every hit in 3 frames with hid=1,-1 using the get neighbour function, which also works if hit is at the edge of the detector.
def build_mask_around_cluster_master_3_frames(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    #select reconstruction/tracking method
    if rec_type == "pixelpixel":
        mask_masters = melp.clustering.tracking.get_mask_masters_hitAnglePixelRec_without_hid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, matching="nearest")
    else:
        mask_masters = melp.clustering.spatial_cluster.get_cluster_primary_truth_3frames(ttree_mu3e, ttree_mu3e_mc, frame)

    mask = {}
    if mask_type == "small":
        for i in range(len(mask_masters)):
            tile_centre = [mask_masters[i][0], mask_masters[i][1]]
            tile_centre_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "up"), mask_masters[i][1]]
            tile_centre_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "down"), mask_masters[i][1]]
            tile_left_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "left"), mask_masters[i][1]]
            tile_right_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "right"), mask_masters[i][1]]

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, 
                                tile_left_centre, tile_right_centre])

            if False in mask_tmp:
                mask_tmp = np.array([list(x) for x in mask_tmp if x[0] != False])

            mask[tile_centre[0]] = mask_tmp

    if mask_type == "medium":
        for i in range(len(mask_masters)):
            tile_centre = [mask_masters[i][0], mask_masters[i][1]]
            tile_centre_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "up"), mask_masters[i][1]]
            tile_centre_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "down"), mask_masters[i][1]]
            tile_left_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "left"), mask_masters[i][1]]
            tile_left_top = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "up"), mask_masters[i][1]]
            tile_left_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "down"), mask_masters[i][1]]
            tile_right_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "right"), mask_masters[i][1]]
            tile_right_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "up"), mask_masters[i][1]]
            tile_right_bootom = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "down"), mask_masters[i][1]]

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, 
                                tile_left_centre, tile_left_bottom, tile_right_top, tile_right_centre, 
                                tile_right_bootom])

            if False in mask_tmp:
                mask_tmp = np.array([list(x) for x in mask_tmp if x[0] != False])

            mask[tile_centre[0]] = mask_tmp

    if mask_type == "big":
        for i in range(len(mask_masters)):
            tile_centre = [mask_masters[i][0], mask_masters[i][1]]
            tile_centre_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "up"), mask_masters[i][1]]
            tile_centre_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "down"), mask_masters[i][1]]
            tile_left_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "left"), mask_masters[i][1]]
            tile_left_top =  [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "up"), mask_masters[i][1]]
            tile_left_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "down"), mask_masters[i][1]]
            tile_right_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "right"), mask_masters[i][1]]
            tile_right_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "up"), mask_masters[i][1]]
            tile_right_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "down"), mask_masters[i][1]]

            tile_centre_far_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre_top[0], "up"), mask_masters[i][1]]
            tile_centre_far_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre_bottom[0], "down"), mask_masters[i][1]]
            tile_left_far_centre = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "left"), mask_masters[i][1]]
            tile_left_far_top =  [mu3e_detector.TileDetector.getNeighbour(tile_left_top[0], "up"), mask_masters[i][1]]
            tile_left_far_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_bottom[0], "down"), mask_masters[i][1]]
            tile_right_far_centre = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "right"), mask_masters[i][1]]
            tile_right_far_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_top[0], "up"), mask_masters[i][1]]
            tile_right_far_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_right_bottom[0], "down"), mask_masters[i][1]]

            tile_far_left_top = [mu3e_detector.TileDetector.getNeighbour(tile_left_top[0], "left"), mask_masters[i][1]]
            tile_far_left_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_bottom[0], "left"), mask_masters[i][1]]
            tile_far_right_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_top[0], "right"), mask_masters[i][1]]
            tile_far_right_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_right_bottom[0], "right"), mask_masters[i][1]]

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, tile_left_centre, tile_left_bottom,
                                   tile_right_top, tile_right_centre, tile_right_bottom, tile_centre_far_top, tile_centre_far_bottom, 
                                   tile_left_far_centre, tile_left_far_top, tile_left_far_bottom, tile_right_far_centre, 
                                   tile_right_far_top, tile_right_far_bottom, tile_far_left_top, tile_far_left_bottom, tile_far_right_top, 
                                   tile_far_right_bottom])
            
            if False in mask_tmp:
                mask_tmp = np.array([list(x) for x in mask_tmp if x[0] != False])

            mask[tile_centre[0]] = mask_tmp

    return mask 