import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

from melp import clustering as clump

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster

#########################
#simple threshold 1d time clustering 
def time_clustering_frame(ttree_mu3e, frame, printing = None):
    time_clusters = []

    #-------------------------------------------------------------
    #set maximum time between hits to get assigned to same cluster
    #-------------------------------------------------------------
    # 0.175ns is ideal for combined time and spatial clustering
    # 0.4ns is ideal for pure time clustering
    time_threshold = 0.4 #ns

    #get hittimes (and hit tiles) in frame
    hittimes_frame, primaries_frame = hittimes_and_primaries_in_frame (ttree_mu3e)

    #--------------
    #build clusters
    #--------------
    added_tiles = []
    for key1 in hittimes_frame.keys():
        if key1 not in added_tiles:
            time_cluster_tmp = []
            time_cluster_tmp.append(ClusterHit(tile_id=key1, time=hittimes_frame[key1][0], frame_id=frame, primary=primaries_frame[key1][0]))
            added_tiles.append(key1)
            for key2 in hittimes_frame.keys():
                if key2 not in added_tiles:
                    if np.abs(hittimes_frame[key1][0] - hittimes_frame[key2][0]) < time_threshold and key2 != key1:
                        time_cluster_tmp.append(ClusterHit(tile_id=key2, time=hittimes_frame[key2][0], frame_id=frame, primary=primaries_frame[key2][0]))
                        added_tiles.append(key2)
            #select shortest time in cluster as master 
            cluster_tmp = Cluster(id=-1, master_id = -1, frame_id=frame, hits=time_cluster_tmp)
            min_hit, __ = cluster_tmp.get_min_time()
            time_clusters.append(Cluster(id=min_hit.tile_id, master_id=min_hit.tile_id, master_primary=min_hit.primary, frame_id=frame, hits=time_cluster_tmp))

    #---------------------------
    #get highest and lowest time
    #---------------------------
    hittimes_frame, __ = hittimes_and_primaries_in_frame (ttree_mu3e)
    times_arr = []
    for key in hittimes_frame.keys():
        times_arr.append(hittimes_frame[key][0])

    if printing == True:
        print("Highest time: ", np.max(times_arr), " ns")
        print("Lowest time: ", np.min(times_arr), " ns")

    #------------------------------------------------------------------------------------
    #get average time difference between hits and the minimum and maximum time difference
    #------------------------------------------------------------------------------------
    times_arr_sorted = np.flip(np.sort(times_arr))
    delta_t = []
    for i in range(len(times_arr_sorted)-1):
        delta_t.append(times_arr_sorted[i]- times_arr_sorted[i+1])

    if printing == True:
        print("Average time difference: ", np.sum(delta_t)/len(delta_t), " ns")
        print("Minimum time difference: ", np.min(delta_t), " ns")
        print("Maximum time difference: ", np.max(delta_t), " ns")
        print("Sorted time differences: ", np.sort(delta_t))

    return time_clusters


###########################
#returns clusters from spatial_cluster with timestamp
def add_hit_time_to_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    time_spatial_clusters = {} #gets returned. Key: cluster_master_tile, Values: whole cluster with time info [[tile_id,time],[tile_id, time],...]
    #get spatial clusters
    spatial_clusters = clump.spatial_cluster.build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)

    #get hittimes for all hits in frame
    hittimes_frame, __ = hittimes_and_primaries_in_frame(ttree_mu3e)

    #--------------------------------
    #get hit time for hits in cluster
    #--------------------------------
    for key in spatial_clusters.keys():
        time_spatial_clusters_tmp = []
        for hit_tile in spatial_clusters[key]:
            if hit_tile in hittimes_frame.keys():
                hit_time = hittimes_frame[hit_tile][0]
                time_spatial_clusters_tmp.append([hit_tile, hit_time])

        time_spatial_clusters[key] = time_spatial_clusters_tmp

    return time_spatial_clusters


