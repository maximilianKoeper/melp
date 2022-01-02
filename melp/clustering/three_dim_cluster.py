import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

from melp import clustering as clump

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster

#######################################
def get_cluster_master_of_time_clusters(ttree_mu3e, frame, printing = None):
    #get time clusters
    time_clusters = clump.time_cluster.time_clustering_frame(ttree_mu3e, frame, printing = None)

    #get master and their primaries from time clusters
    time_cluster_masters =          []
    time_cluster_master_primaries = []
    for time_cluster in time_clusters:
        time_cluster_masters.append(time_cluster.master_id)
        time_cluster_master_primaries.append(time_cluster.master_primary)

    return time_cluster_masters, time_cluster_master_primaries



########################################
def spatial_clustering_for_time_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    #use time cluster as first rough cut / reference
    time_clusters = clump.time_cluster.time_clustering_frame(ttree_mu3e, frame, printing = None)

    #apply spatial cluster to those time clusters
    clusters = []
    #--------------------------------------------
    #get masks around master tile of time cluster
    #--------------------------------------------
    master_masks, master_primaries = clump.masks.build_mask_around_cluster_master(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type = "timethenspatial")


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



    #----------------------------------------------
    #check if spatial clusters are in time clusters
    #----------------------------------------------
    combined_clusters = time_clusters
    new_clusters = []
    #for time_cluster in time_clusters:
    for i in range(len(time_clusters)):
        #print("i = ",i)
        hits_not_in_spatial_clusters = [] #should be seperate cluster
        for hit in time_clusters[i].hits:
        #for j in range(time_clusters[i].__len__()-1):
        #    print("j = ",j)
            found_counter = 0
            #if hit in hits_not_in_spatial_clusters:
            #    continue
            for cluster in clusters:
                if hit in cluster.hits:
                #if time_clusters[i].hits[j] in cluster.hits:
                    found_counter += 1
                    break
            
            if found_counter == 0:
                hits_not_in_spatial_clusters.append(hit)
                #time_cluster.hits.remove(hit) #TODO: Check how to remove properly
                #del combined_clusters[i].hits[j]
                combined_clusters[i].hits.remove(hit)

        if len(hits_not_in_spatial_clusters) != 0:
            new_clusters.append(hits_not_in_spatial_clusters)

    #----------------------------------------
    #append new seperate clusters to clusters
    #----------------------------------------
    for i in range(len(new_clusters)):
        #time_clusters.append(Cluster(id=-1, master_id=-1, master_primary=-1, frame_id=frame, hits=new_clusters[i]))
        combined_clusters.append(Cluster(id=-1, master_id=-1, master_primary=-1, frame_id=frame, hits=new_clusters[i]))

    return combined_clusters #time_clusters

    


