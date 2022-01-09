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
    for i in range(len(time_clusters)):
        hits_not_in_spatial_clusters = [] #should be separate cluster
        for hit in time_clusters[i].hits:
            found_counter = 0
            for cluster in clusters:
                if hit in cluster.hits:
                    found_counter += 1
                    break

            if found_counter == 0:
                hits_not_in_spatial_clusters.append(hit)
                combined_clusters[i].hits.remove(hit)

        if len(hits_not_in_spatial_clusters) != 0:
            new_clusters.append(hits_not_in_spatial_clusters)

    #----------------------------------------
    #append new seperate clusters to clusters
    #----------------------------------------
    for i in range(len(new_clusters)):
        combined_clusters.append(Cluster(id=-1, master_id=-1, master_primary=-1, frame_id=frame, hits=new_clusters[i]))

    return combined_clusters 

    
#########################################
def iterative_masks_after_time_clustering(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    #-----------------------------------------------
    #use time cluster as first rough cut / reference
    #-----------------------------------------------
    time_clusters = clump.time_cluster.time_clustering_frame(ttree_mu3e, frame, printing = None)

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

    #-------------------------
    #build iterative masks
    #-------------------------
    added_hits = [] #just tile_ids
    new_clusters = []
    for time_cluster in time_clusters:  #loop over all clusters
        if time_cluster.__len__() == 1:
            new_clusters.append(time_cluster)
            added_hits.append(time_cluster.hits[0].tile_id)
        else:
            for i in range(len(time_cluster.hits)):  #loop over all hits in cluster  
                cluster_tmp = []  
                if time_cluster.hits[i].tile_id not in added_hits:
                    mask_tmp, master_primary_mask = clump.masks.build_mask_single_hit(time_cluster.hits[i], ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type = None)
                    #cluster_tmp = [] 
                    for j in range(len(hit_tiles_frame)):
                        if hit_tiles_frame[j] in mask_tmp and hit_tiles_frame[j] not in added_hits: 
                            cluster_tmp.append(ClusterHit(tile_id = hit_tiles_frame[j], frame_id = frame, primary = primaries_frame[j], time = times_frame[j]))
                            added_hits.append(hit_tiles_frame[j])
                    #build mask around hits in first cluster
                    if len(cluster_tmp) != 0:
                        for hit_tmp in cluster_tmp:
                            next_mask_tmp, __ = clump.masks.build_mask_single_hit(hit_tmp, ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type = None)
                            for m in range(len(hit_tiles_frame)):
                                if hit_tiles_frame[m] in next_mask_tmp and hit_tiles_frame[m] not in added_hits: 
                                    cluster_tmp.append(ClusterHit(tile_id = hit_tiles_frame[m], frame_id = frame, primary = primaries_frame[m], time = times_frame[m]))
                                    added_hits.append(hit_tiles_frame[m])

                #add leftover hits as separate cluster
                separate_cluster_tmp = []
                for n in range(len(time_cluster.hits)):
                    tile_ids_tmp = []
                    for hit_tmp in cluster_tmp:
                        tile_ids_tmp.append(hit_tmp.tile_id)
                    if time_cluster.hits[n].tile_id not in tile_ids_tmp and time_cluster.hits[n].tile_id not in added_hits: 
                        separate_cluster_tmp.append(time_cluster.hits[n])
                        added_hits.append(time_cluster.hits[n].tile_id)
                if len(separate_cluster_tmp) != 0:
                    new_clusters.append(Cluster(id = separate_cluster_tmp[0].tile_id, master_id = separate_cluster_tmp[0].tile_id, master_primary = separate_cluster_tmp[0].primary, frame_id = frame, hits = separate_cluster_tmp))


                if len(cluster_tmp) != 0:
                    new_clusters.append(Cluster(id = time_cluster.master_id, master_id = time_cluster.master_id, master_primary = master_primary_mask, frame_id = frame, hits = cluster_tmp))

    #######################################
    if frame < 100:
        cluster_hits_counter = 0
        for cluster in new_clusters:
            cluster_hits_counter += cluster.__len__()
        if cluster_hits_counter > len(hit_tiles_frame):
            print("Frame:", frame, "#Total hits:", len(hit_tiles_frame), "#Cluster hits:", cluster_hits_counter)
    
        

        #double_hit_counter = 0
        #all_hits_clusters = []
        #double_hits = []
        #for cluster in new_clusters:
        #    for hit in cluster.hits:
        #        all_hits_clusters.append(hit)
        #for i in all_hits_clusters:
        #    if all_hits_clusters.count(i) > 1:
        #        double_hit_counter +=1
        #        double_hits.append(i)
        
        #test_clusters = [Cluster(id = 0, master_id = 0, master_primary = 0, frame_id = frame, hits = double_hits)]
        #if double_hit_counter != 0:
        #    print("Frame:", frame, "#Total hits:", len(hit_tiles_frame), "#Double hits:", double_hit_counter, "\n")
    ###########################################

    return new_clusters

