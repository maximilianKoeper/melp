import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

from melp import clustering as clump

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster

#######################################
def get_cluster_master_of_time_clusters(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold, printing = None):
    #get time clusters
    #time_clusters = clump.time_cluster.time_clustering_frame(ttree_mu3e, frame, printing = None)
    time_clusters = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)

    #get master and their primaries from time clusters
    time_cluster_masters          = []
    time_cluster_master_primaries = []
    time_cluster_master_tids      = []
    for time_cluster in time_clusters:
        time_cluster_masters.append(time_cluster.master_id)
        time_cluster_master_primaries.append(time_cluster.master_primary)

    return time_cluster_masters, time_cluster_master_primaries



########################################
def spatial_clustering_for_time_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, time_threshold, mask_type, rec_type = None):
    #use time cluster as first rough cut / reference
    time_clusters = clump.time_cluster.time_clustering_frame(ttree_mu3e, ttree_mu3e_mc, frame, printing = None)

    #apply spatial cluster to those time clusters
    clusters = []
    #--------------------------------------------
    #get masks around master tile of time cluster
    #--------------------------------------------
    master_masks, master_primaries = clump.masks.build_mask_around_cluster_master(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, time_threshold, mask_type, rec_type = "timethenspatial")


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
def iterative_masks_after_time_clustering(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, time_threshold, mask_type, rec_type = None):
    #-----------------------------------------------
    #use time cluster as first rough cut / reference
    #-----------------------------------------------
    #without energy cut
    time_clusters = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)
    #with energy cut
    #time_clusters = clump.time_cluster.time_clustering_frame_improv_energy_cut(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)

    #-----------------------------------------------------------------------
    #get all tiles that have been hit in frame and their primaries and times
    #-----------------------------------------------------------------------
    hit_tiles_frame = []
    primaries_frame = []
    times_frame     = []
    mcis_frame      = []
    tids_frame      = []
    edep_frame      = []
    pdgs_frame      = []
    
    traj_PID_frame      = []
    traj_type_frame     = []
    traj_tlhid_frame    = []
    traj_ID_frame       = []

    for hit_tile_index in range(ttree_mu3e.Ntilehit):
        hit_tiles_frame.append(ttree_mu3e.tilehit_tile[hit_tile_index])
        primaries_frame.append(ttree_mu3e.tilehit_primary[hit_tile_index])
        times_frame.append(ttree_mu3e.tilehit_time[hit_tile_index])
        mcis_frame.append(ttree_mu3e.tilehit_mc_i[hit_tile_index])
        edep_frame.append(ttree_mu3e.tilehit_edep[hit_tile_index])
        traj_PID_frame.append(ttree_mu3e.traj_PID[hit_tile_index])
        traj_type_frame.append(ttree_mu3e.traj_type[hit_tile_index])
        traj_tlhid_frame.append(ttree_mu3e.traj_tlhid[hit_tile_index])
        traj_ID_frame.append(ttree_mu3e.traj_ID[hit_tile_index])
    
    for mc_i in mcis_frame:
        ttree_mu3e_mc.GetEntry(mc_i)
        tids_frame.append(ttree_mu3e_mc.tid)
        pdgs_frame.append(ttree_mu3e_mc.pdg)

    #-------------------------
    #build iterative masks
    #-------------------------
    new_clusters = []
    for time_cluster in time_clusters:  #loop over all clusters
        added_hits = [] #just tile_ids
        for i in range(len(time_cluster.hits)):  #loop over all hits in cluster  
            cluster_tmp = []  
            remaining_hits = []
            if time_cluster.hits[i].tile_id not in added_hits:
                mask_tmp, master_primary_mask, master_tid_mask = clump.masks.build_mask_single_hit(time_cluster.hits[i], ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, 
                                                                                                   mu3e_detector, frame, mask_type, rec_type = None)
                for j in range(len(time_cluster.hits)):
                    if time_cluster.hits[j].tile_id in mask_tmp and time_cluster.hits[j].tile_id not in added_hits: 
                        cluster_tmp.append(ClusterHit(tile_id = time_cluster.hits[j].tile_id, frame_id = frame, primary = time_cluster.hits[j].primary, 
                                                      time = time_cluster.hits[j].time, mc_i = time_cluster.hits[j].mc_i, tid = time_cluster.hits[j].tid, 
                                                      pdg = time_cluster.hits[j].pdg ,edep = time_cluster.hits[j].edep, traj_PID = time_cluster.hits[j].traj_PID, 
                                                      traj_type = time_cluster.hits[j].traj_type, traj_tlhid = time_cluster.hits[j].traj_tlhid, 
                                                      traj_ID = time_cluster.hits[j].traj_ID))
                        added_hits.append(time_cluster.hits[j].tile_id)
            #build mask around hits in first cluster
            cluster_tmp_2 = []
            if len(cluster_tmp) != 0:
                for hit_tmp in cluster_tmp:
                    next_mask_tmp, __, __ = clump.masks.build_mask_single_hit(hit_tmp, ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, 
                                                                              mu3e_detector, frame, mask_type, rec_type = None)
                    for m in range(len(time_cluster.hits)):
                        if time_cluster.hits[m].tile_id in next_mask_tmp and time_cluster.hits[m].tile_id not in added_hits: 
                            cluster_tmp.append(ClusterHit(tile_id = time_cluster.hits[m].tile_id, frame_id = frame, primary = time_cluster.hits[m].primary, 
                                                          time = time_cluster.hits[m].time, mc_i = time_cluster.hits[m].mc_i, tid = time_cluster.hits[m].tid, 
                                                          pdg = time_cluster.hits[m].pdg, edep = time_cluster.hits[m].edep, traj_PID = time_cluster.hits[m].traj_PID, 
                                                          traj_type = time_cluster.hits[m].traj_type, traj_tlhid = time_cluster.hits[m].traj_tlhid, 
                                                          traj_ID = time_cluster.hits[m].traj_ID))
                            cluster_tmp_2.append(ClusterHit(tile_id = time_cluster.hits[m].tile_id, frame_id = frame, primary = time_cluster.hits[m].primary, 
                                                            time = time_cluster.hits[m].time, mc_i = time_cluster.hits[m].mc_i, tid = time_cluster.hits[m].tid, 
                                                            pdg = time_cluster.hits[m].pdg, edep = time_cluster.hits[m].edep, traj_PID = time_cluster.hits[m].traj_PID, 
                                                            traj_type = time_cluster.hits[m].traj_type, traj_tlhid = time_cluster.hits[m].traj_tlhid, 
                                                            traj_ID = time_cluster.hits[m].traj_ID))
                            added_hits.append(time_cluster.hits[m].tile_id)

            #build mask around hits in second iteration clusters
            #cluster_tmp_3 = []
            #if len(cluster_tmp_2) != 0:
            #    for hit_tmp_2 in cluster_tmp_2:
            #        next_mask_tmp_2, __, __ = clump.masks.build_mask_single_hit(hit_tmp_2, ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type = None)
            #        for m in range(len(time_cluster.hits)):
            #            if time_cluster.hits[m].tile_id in next_mask_tmp_2 and time_cluster.hits[m].tile_id not in added_hits: 
            #                cluster_tmp.append(ClusterHit(tile_id = time_cluster.hits[m].tile_id, frame_id = frame, primary = time_cluster.hits[m].primary, time = time_cluster.hits[m].time, mc_i = time_cluster.hits[m].mc_i, tid = time_cluster.hits[m].tid))
            #                cluster_tmp_3.append(ClusterHit(tile_id = time_cluster.hits[m].tile_id, frame_id = frame, primary = time_cluster.hits[m].primary, time = time_cluster.hits[m].time, mc_i = time_cluster.hits[m].mc_i, tid = time_cluster.hits[m].tid))
            #                added_hits.append(time_cluster.hits[m].tile_id)

            if len(cluster_tmp) != 0:
                cluster_tmp_times = []
                for hit in cluster_tmp:
                    cluster_tmp_times.append(hit.time)
                index_min_time = np.argmin(cluster_tmp_times)
                new_clusters.append(Cluster(id = cluster_tmp[index_min_time].tile_id, master_id = cluster_tmp[index_min_time].tile_id, 
                                            master_primary = cluster_tmp[index_min_time].primary, master_tid = cluster_tmp[index_min_time].tid, 
                                            frame_id = frame, hits = cluster_tmp))


    """
    #----------------------------------------------------------------------------------------------
    #Apply more sensitive time cut (0.2ns) for the new clusters from tile to the neighbouring tiles
    #----------------------------------------------------------------------------------------------
    new_clusters_2 = []
    for n in range(len(new_clusters)):
        if len(new_clusters[n]) == 2:
            new_cluster_2_counter_tmp = 0 
            for i in range(1, len(new_clusters[n])):
                mask_tmp_2, __, __ = clump.masks.build_mask_single_hit(new_clusters[n].hits[i], ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type = "medium", rec_type = None)
                for j in range(len(new_clusters[n])):
                    if j != i and new_clusters[n].hits[j].tile_id in mask_tmp_2:
                        if abs(new_clusters[n].hits[i].time - new_clusters[n].hits[j].time) > 0.3 and new_cluster_2_counter_tmp == 0:
                            new_clusters_2.append(Cluster(id = new_clusters[n].hits[i].tile_id, master_id = new_clusters[n].hits[i].tile_id, master_primary = new_clusters[n].hits[i].primary, master_tid = new_clusters[n].hits[i].tid, frame_id = frame, hits = [new_clusters[n].hits[i]]))
                            new_clusters_2.append(Cluster(id = new_clusters[n].hits[j].tile_id, master_id = new_clusters[n].hits[j].tile_id, master_primary = new_clusters[n].hits[j].primary, master_tid = new_clusters[n].hits[j].tid, frame_id = frame, hits = [new_clusters[n].hits[j]]))
                            new_cluster_2_counter_tmp += 1
                            pass
            if new_cluster_2_counter_tmp == 0:
                new_clusters_2.append(new_clusters[n])
        else:
            new_clusters_2.append(new_clusters[n])
    """

    """
    #-----------------------------------------------------
    #check for cluster constituents in neighbouring frames 
    #-----------------------------------------------------
    #get first and last cluster in time
    min_times = []
    max_times = []
    for i in range(len(new_clusters_2)):
        min_times.append(min(new_clusters_2[i].get_times()))
        max_times.append(max(new_clusters_2[i].get_times()))

    if len(min_times) > 0 and len(max_times) > 0:
        min_time = min(min_times)
        max_time = max(max_times)
        cluster_min_time_index = min_times.index(min(min_times))
        cluster_max_time_index = max_times.index(max(max_times))

        #get clusters in frame-1 and frame+1
        ttree_mu3e.GetEntry(frame-1)
        time_clusters_minus = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame-1, time_threshold)
        ttree_mu3e.GetEntry(frame+1)
        time_clusters_plus = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame+1, time_threshold)
        ttree_mu3e.GetEntry(frame)

        #check if the maximum time of frame-1 is close to the min time of frame 0
        max_times_minus = []
        for i in range(len(time_clusters_minus)):
            max_times_minus.append(max(time_clusters_minus[i].get_times()))
        max_time_minus = max(max_times_minus)
        cluster_max_time_minus_index = max_times_minus.index(max(max_times_minus))

        #check if the minimum time of frame+1 is close to the max time of frame 0
        min_times_plus = []
        for i in range(len(time_clusters_plus)):
            min_times_plus.append(min(time_clusters_plus[i].get_times()))
        min_time_plus = min(min_times_plus)
        cluster_min_time_plus_index = min_times_plus.index(min(min_times_plus))

        #testing
        test = abs(abs(min_time - max_time_minus) - 0)
        if test < 0.3:
            print("Frame -1:", test)

        test2 = abs(abs(max_time - min_time_plus) - 0)
        if test2 < 0.3:
            print("Frame +1:", test2)
    """
    
    return new_clusters

