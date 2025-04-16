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
def time_clustering_frame(ttree_mu3e, ttree_mu3e_mc, frame, printing = None):
    time_clusters = []

    #-------------------------------------------------------------
    #set maximum time between hits to get assigned to same cluster
    #-------------------------------------------------------------
    # 0.175ns is ideal for combined time and spatial clustering
    # 0.4ns is ideal for pure time clustering
    time_threshold = 1.2 #ns

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

################################
def time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame: int, time_threshold: float = 0.4) -> dict:
    indices           = np.argsort(list(ttree_mu3e.tilehit_time))
    tilehit_times     = np.asarray(list(ttree_mu3e.tilehit_time))[indices]
    tilehit_ids       = np.asarray(list(ttree_mu3e.tilehit_tile))[indices]
    tilehit_primaries = np.asarray(list(ttree_mu3e.tilehit_primary))[indices]
    tilehit_mcis      = np.asarray(list(ttree_mu3e.tilehit_mc_i))[indices]
    tilehit_edep      = np.asarray(list(ttree_mu3e.tilehit_edep))[indices]

    if len(tilehit_times) == 0:
        return [Cluster(id=-1, frame_id=frame, master_id=-1, master_primary=-1, hits=[])]
    
    else:
        clusters_dict = {}
        tmp_time_reference = tilehit_times[0]
        index_start_track = 0
        for index in range(len(tilehit_times)):
            if abs(tmp_time_reference - tilehit_times[index]) > time_threshold:
                clusters_dict[index] = [tilehit_ids[index_start_track:index], tilehit_primaries[index_start_track:index], 
                                        tilehit_times[index_start_track:index], tilehit_mcis[index_start_track:index], 
                                        tilehit_edep[index_start_track:index]]
                index_start_track = index
                tmp_time_reference = tilehit_times[index]

        #fill up remaining event
        if index_start_track != len(tilehit_times):
            clusters_dict[len(tilehit_times)] = [tilehit_ids[index_start_track:],
                                                tilehit_primaries[index_start_track:],
                                                tilehit_times[index_start_track:],
                                                tilehit_mcis[index_start_track:],
                                                tilehit_edep[index_start_track:]]

        #convert to cluster object
        clusters = []
        for key in clusters_dict.keys():
            cluster_tmp = []
            for i in range(len(clusters_dict[key][0])):
                mc_i = clusters_dict[key][3][i]
                ttree_mu3e_mc.GetEntry(mc_i)
                tid = ttree_mu3e_mc.tid
                pdg = ttree_mu3e_mc.pdg
                cluster_tmp.append(ClusterHit(tile_id=clusters_dict[key][0][i], frame_id=frame, primary=clusters_dict[key][1][i], 
                                              time=clusters_dict[key][2][i], mc_i=clusters_dict[key][3][i], tid = tid, pdg = pdg, 
                                              edep = clusters_dict[key][4][i]))                
            clusters.append(Cluster(id=key, frame_id=frame, master_id=cluster_tmp[0].tile_id, master_primary=cluster_tmp[0].primary, 
                                    master_tid = cluster_tmp[0].tid, hits=cluster_tmp))

    return clusters


################################
def time_clustering_frame_improv_energy_cut(ttree_mu3e, ttree_mu3e_mc, frame: int, time_threshold: float = 0.4) -> dict:
    indices           = np.argsort(list(ttree_mu3e.tilehit_time))
    tilehit_times     = np.asarray(list(ttree_mu3e.tilehit_time))[indices]
    tilehit_ids       = np.asarray(list(ttree_mu3e.tilehit_tile))[indices]
    tilehit_primaries = np.asarray(list(ttree_mu3e.tilehit_primary))[indices]
    tilehit_mcis      = np.asarray(list(ttree_mu3e.tilehit_mc_i))[indices]
    tilehit_edep      = np.asarray(list(ttree_mu3e.tilehit_edep))[indices]

    if len(tilehit_times) == 0:
        return [Cluster(id=-1, frame_id=frame, master_id=-1, master_primary=-1, hits=[])]
    
    else:
        clusters_dict = {}
        tmp_time_reference = tilehit_times[0]
        index_start_track = 0
        for index in range(len(tilehit_times)):
            if abs(tmp_time_reference - tilehit_times[index]) > time_threshold:
                clusters_dict[index] = [tilehit_ids[index_start_track:index], 
                                        tilehit_primaries[index_start_track:index], 
                                        tilehit_times[index_start_track:index], 
                                        tilehit_mcis[index_start_track:index], 
                                        tilehit_edep[index_start_track:index]]
                index_start_track = index
                tmp_time_reference = tilehit_times[index]

        #-----------------------
        #fill up remaining event
        #-----------------------
        if index_start_track != len(tilehit_times):
            clusters_dict[len(tilehit_times)] = [tilehit_ids[index_start_track:],
                                                tilehit_primaries[index_start_track:],
                                                tilehit_times[index_start_track:],
                                                tilehit_mcis[index_start_track:],
                                                tilehit_edep[index_start_track:]]

        #-------------------------
        #convert to cluster object
        #-------------------------
        clusters = []
        for key in clusters_dict.keys():
            cluster_tmp = []
            for i in range(len(clusters_dict[key][0])):
                mc_i = clusters_dict[key][3][i]
                ttree_mu3e_mc.GetEntry(mc_i)
                tid = ttree_mu3e_mc.tid
                pdg = ttree_mu3e_mc.pdg
                cluster_tmp.append(ClusterHit(tile_id=clusters_dict[key][0][i], frame_id=frame, primary=clusters_dict[key][1][i], 
                                              time=clusters_dict[key][2][i], mc_i=clusters_dict[key][3][i], tid = tid, pdg = pdg,
                                              edep = clusters_dict[key][4][i]))                
            clusters.append(Cluster(id=key, frame_id=frame, master_id=cluster_tmp[0].tile_id, master_primary=cluster_tmp[0].primary, 
                                    master_tid = cluster_tmp[0].tid, hits=cluster_tmp))

        #----------------
        #apply energy cut
        #----------------
        for cluster in clusters:
            for cluster_hit in cluster.hits:
                if cluster_hit.edep < 0.3: #MeV
                    cluster.hits.remove(cluster_hit)

    return clusters


#############################
#returns average number of time cluster hits
def average_number_of_cluster_hits(ttree_mu3e, ttree_mu3e_mc, number_of_frames: int, time_threshold: float = 0.4):
    hits_in_clusters            = []
    number_of_clusters_in_frame = []
    number_of_hits_in_frame     = []
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames

    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #get clusters
        clusters = time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)

        #get number of clusters in frame
        number_of_clusters_in_frame.append(len(clusters))

        #get number of tile hits in frame
        number_of_hits_in_frame.append(ttree_mu3e.Ntilehit)

        #calculate number of hits
        for cluster in clusters:
            hits_in_clusters.append(len(cluster))

    print("Progress: 100 %","of ", frames_to_analyze, " frames")

    #calculate average number of hits and clusters
    mean_cluster_hits   = np.mean(hits_in_clusters)
    mean_clusters_frame = np.mean(number_of_clusters_in_frame)
    mean_hits_frame     = np.mean(number_of_hits_in_frame)

    return mean_cluster_hits, mean_clusters_frame, mean_hits_frame



