import ROOT
import numpy as np

import melp
import melp.clustering as clump
from melp.clustering.misc import*

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster


###########################
#builds clusters in the masks around hit with hid=1,-1 according to primaries.
def build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    clusters = []
    #----------------------------
    #get masks around master tile
    #----------------------------
    master_masks, master_primaries = clump.masks.build_mask_around_cluster_master(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type)

    #-----------------------------------------------------------------------
    #get all tiles that have been hit in frame and their primaries and times
    #-----------------------------------------------------------------------
    hit_tiles_frame = []
    primaries_frame = []
    times_frame     = []
    for hit_tile_index in range(ttree_mu3e.Ntilehit):
        hit_tiles_frame.append(ttree_mu3e.tilehit_tile[hit_tile_index])
        primaries_frame.append(ttree_mu3e.tilehit_primary[hit_tile_index])
        times_frame.append(ttree_mu3e.tilehit_time[hit_tile_index])
    
    #---------------------------------------
    #build clusters around mask master tiles
    #---------------------------------------
    for key in master_masks.keys():
        cluster_tmp     = []
        cluster_tmp_ids = []
        i = 0
        while i in range(len(hit_tiles_frame)):
            if hit_tiles_frame[i] != "associated" and hit_tiles_frame[i] in master_masks[key]:
                cluster_tmp.append(ClusterHit(tile_id = hit_tiles_frame[i], frame_id = frame, primary = primaries_frame[i], time = times_frame[i]))
                cluster_tmp_ids.append(hit_tiles_frame[i])
                hit_tiles_frame[i] = "associated" #prevent hits from being in multiple clusters
            i += 1    

        if len(cluster_tmp) != 0 and key in cluster_tmp_ids:
            master_id         = key
            primary_of_master = master_primaries[key]

            clusters.append(Cluster(id = master_id, master_id=key, master_primary = primary_of_master, frame_id = frame, hits = cluster_tmp))

    return clusters


#######################
def build_truth_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    #-----------------------------------------------------------------------
    #get all tiles that have been hit in frame and their primaries and times
    #-----------------------------------------------------------------------
    hit_tiles_frame = []
    primaries_frame = []
    times_frame     = []
    mcis_frame      = []
    tids_frame      = []
    for hit_tile_index in range(ttree_mu3e.Ntilehit):
        hit_tiles_frame.append(ttree_mu3e.tilehit_tile[hit_tile_index])
        primaries_frame.append(ttree_mu3e.tilehit_primary[hit_tile_index])
        times_frame.append(ttree_mu3e.tilehit_time[hit_tile_index])
        mcis_frame.append(ttree_mu3e.tilehit_mc_i[hit_tile_index])
    
    for mc_i in mcis_frame:
        ttree_mu3e_mc.GetEntry(mc_i)
        tids_frame.append(ttree_mu3e_mc.tid)

    #---------------------------------------------
    #sort hits with same primary into same cluster
    #---------------------------------------------
    #clusters = []
    #added_primaries = []
    #for i in range(len(primaries_frame)):
    #    cluster_tmp = []
    #    for j in range(len(primaries_frame)):
    #        if primaries_frame[j] == primaries_frame[i] and primaries_frame[i] not in added_primaries:
    #            cluster_tmp.append(ClusterHit(tile_id = hit_tiles_frame[j], frame_id = frame, primary = primaries_frame[j], tid = tids_frame[j], time = times_frame[j]))
    #    added_primaries.append(primaries_frame[i])

    #    if len(cluster_tmp) > 0:
    #        clusters.append(Cluster(id=i, master_id=i, master_primary = primaries_frame[i], master_tid = tids_frame[i], frame_id = frame, hits = cluster_tmp))

    #---------------------------------------------
    #sort hits with same tid into same cluster
    #---------------------------------------------
    clusters = []
    added_tids = []
    for i in range(len(tids_frame)):
        cluster_tmp = []
        for j in range(len(tids_frame)):
            if tids_frame[j] == tids_frame[i] and tids_frame[i] not in added_tids:
                cluster_tmp.append(ClusterHit(tile_id = hit_tiles_frame[j], frame_id = frame, primary = primaries_frame[j], tid = tids_frame[j], time = times_frame[j]))
        added_tids.append(tids_frame[i])

        if len(cluster_tmp) > 0:
            clusters.append(Cluster(id=i, master_id=i, master_primary = primaries_frame[i], master_tid = tids_frame[i], frame_id = frame, hits = cluster_tmp))

    return clusters