import ROOT
import numpy as np

import melp
from melp import Detector
from melp.clustering.misc import*

#----------------------------------------------------
#builds clusters in the masks around hit with hid=1,-1 according to primaries.
def build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    clusters = {} #keys: "master"-tile; values: rest of cluster

    #get masks around master tile
    primary_masks = melp.clustering.masks.build_mask_around_cluster_master(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type)

    keys = []
    values = []
    for key in primary_masks.keys():
        keys.append(key)
        values.append(primary_masks[key])

    #get all tiles that have been hit in frame
    hit_tiles = []
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])

    #build clusters around mask master tiles
    for key in primary_masks.keys():
        cluster_tmp = []
        for tile_id in hit_tiles:
            #if tile_id < 300000: #just left recurl station
            if tile_id in primary_masks[key]:
                cluster_tmp.append(tile_id)
                #hit_tiles.remove(tile_id) #prevent hits from being in multiple clusters
                for _ in range(hit_tiles.count(tile_id)):
                    hit_tiles.remove(tile_id)
        
        if len(cluster_tmp) != 0 and key in cluster_tmp:
            clusters[key] = cluster_tmp

    return clusters

#----------------------------------------------------
#builds clusters in the masks around hit with hid=1,-1 according to primaries.
def build_clusters_in_masks_3_frames(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    clusters = {} #keys: "master"-tile; values: rest of cluster

    #get masks around master tile
    primary_masks = melp.clustering.masks.build_mask_around_cluster_master_3_frames(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type)

    #print(primary_masks)

    keys = []
    values = []
    for key in primary_masks.keys():
        keys.append(key)
        values.append(primary_masks[key])

    #get all tiles that have been hit in frame
    hit_tiles = []
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])

    ttree_mu3e.GetEntry(frame+1)
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])

    ttree_mu3e.GetEntry(frame-1)
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])
    
    ttree_mu3e.GetEntry(frame)

    #build clusters around mask master tiles
    for tile_id in hit_tiles:
        #if tile_id < 300000: #just left recurl station
        cluster_tmp = []
        cluster_primary_tmp = 0
        if tile_id not in keys: #if not primary
            for i in range(len(values)):
                for j in range(len(values[i])):
                    if tile_id == values[i][j][0]:   
                        cluster_tmp.append(list(values[i][j]))
                        if list(values[i][0]) not in cluster_tmp:
                            cluster_tmp.insert(0,list(values[i][0]))
                        cluster_primary_tmp = keys[i]

            if cluster_primary_tmp != 0:
                if cluster_primary_tmp not in clusters.keys() and len(cluster_tmp) > 0:          
                    clusters[cluster_primary_tmp] = cluster_tmp
                elif len(cluster_tmp) > 0: 
                    clusters[cluster_primary_tmp].append(cluster_tmp[1])

    for i in keys:
        if i not in clusters.keys():
            index_tmp = keys.index(i)
            clusters[i] = [list(values[index_tmp][0])]

    return clusters

#-----------------------------------------------------
#builds clusters where dict-key is the primary of "master"-tile and not "master" tile and value is the whole cluster
def build_cluster_with_truth_primary(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector,frame, mask_type, rec_type = None, cluster_type = None):
    #get primaries
    #primaries = get_primary_frame(ttree_mu3e, ttree_mu3e_mc)
    primaries = get_mc_primary_for_hit_frame(ttree_mu3e)

    #get clusters
    if cluster_type == None:
        clusters_frame = build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector,frame, mask_type, rec_type)
    elif cluster_type == "time":
        __ , clusters_frame = melp.clustering.time_cluster.time_clustering_frame(ttree_mu3e, printing = None)

    #convert cluster tiles to primaries
    clusters_with_primaries = {} #gets returned: keys:primary of "master"-tile; values: primary of whole cluster

    for key in clusters_frame.keys(): #loop over all clusters
        primary_whole_clusters_tmp = []
        #primary_whole_clusters_tmp.append(primaries[key])
        for cluster_hit in clusters_frame[key]: #loop over all hits in cluster i
            primary_whole_clusters_tmp.append(primaries[cluster_hit])

        #clusters_with_primaries[primaries[key]] = primary_whole_clusters_tmp
        clusters_with_primaries[key] = primary_whole_clusters_tmp #TODO: Fix issue that efficiency is 100% for small thresholds
                                                                  #since all hits have then the correct tid. Master missing?

    len_counter_clusters = 0
    for key in clusters_frame.keys():
        len_counter_clusters += len(clusters_frame[key])
    len_counter_primary_clusters = 0
    for key in clusters_with_primaries.keys():
        len_counter_primary_clusters += len(clusters_with_primaries[key])

    return clusters_with_primaries

#-------------------------------------------------
#builds dict with tid of master tile as key and tids of cluster as value
def build_cluster_with_truth_tid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector,frame, mask_type, rec_type = None):
    #get tids
    tilehit_tid = get_tid_frame(ttree_mu3e, ttree_mu3e_mc)

    #get clusters
    clusters_frame = build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector,frame, mask_type, rec_type)

    #convert cluster tiles to tids
    clusters_with_tid = {} #gets returned: keys:tid of "master"-tile; values:tids of whole cluster
    for key in clusters_frame.keys():
        tid_whole_clusters_tmp = []
        tid_whole_clusters_tmp.append(tilehit_tid[key])
        for i in clusters_frame[key]:
            tid_whole_clusters_tmp.append(tilehit_tid[i])

        clusters_with_tid[tilehit_tid[key]] = tid_whole_clusters_tmp
    
    return clusters_with_tid

#-----------------------------------------------------
#returns number of hits in cluster and total number of hits
def count_hits_in_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_hits_counter = 0
    cluster_hits_counter = 0

    #counting
    for frame in range(frames_to_analyze):
        ttree_mu3e.GetEntry(frame)

        #Printing status info
        if frame % 2000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')
        
        #count total hits
        tot_hits_frame = len(ttree_mu3e.tilehit_tile)
        total_hits_counter += tot_hits_frame

        #count hits in clusters
        clusters_frame = build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, frame, rec_type)
        
        for key in clusters_frame.keys():
            cluster_hits_counter +=1
            for i in clusters_frame[key]:
                cluster_hits_counter +=1

    print("Progress: 100 %","of ", frames_to_analyze, " frames")

    return cluster_hits_counter, total_hits_counter

#------------------------------------------------------------------
#checks (for 3 frame clustering) if a cluster consists of hits in multiple frames
def check_for_multiple_frame_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_cluster_hits_counter = 0
    mult_frame_cluster_hits_counter = 0

    frac_mult_frame_cluster_hits = []

    #counting
    #for frame in range(frames_to_analyze):
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)

        mult_frame_cluster_hits_counter_tmp = 0
        total_cluster_hits_counter_tmp = 0

        #Printing status info
        if int(frame) % 250 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        clusters_frame = build_clusters_in_masks_3_frames(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)

        for key in clusters_frame.keys():
            value = np.array(clusters_frame[key])
            total_cluster_hits_counter_tmp += len(value)
            for i in range(len(value)): #checks if hit in cluster has different frame id than the "master" hit
                if value[i][1] != value[0][1]: 
                    mult_frame_cluster_hits_counter_tmp += 1

        total_cluster_hits_counter += total_cluster_hits_counter_tmp
        mult_frame_cluster_hits_counter += mult_frame_cluster_hits_counter_tmp

        if total_cluster_hits_counter_tmp != 0:
            frac_mult_frame_cluster_hits.append(mult_frame_cluster_hits_counter_tmp/total_cluster_hits_counter_tmp)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    print("Hits in cluster in different frame than master out of all hits in clusters: ", mult_frame_cluster_hits_counter/(total_cluster_hits_counter/100), "%")

    return frac_mult_frame_cluster_hits, total_cluster_hits_counter, mult_frame_cluster_hits_counter


