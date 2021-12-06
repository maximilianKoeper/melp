import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

import melp.clustering.spatial_cluster as sclump
import melp.clustering.three_frame_cluster as clump_3

#-------------------------------------------------
def compare_to_primary(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    frac_corr_frame = []
    frac_corr_clusters_frame = []
    frac_uncorr_frame = []
    total_hits_counter = 0
    cluster_hits_counter = 0

    tot_corr_counter = 0
    tot_uncorr_counter = 0

    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames

    for frame in range(frames_to_analyze):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        total_hits_frame = len(ttree_mu3e.tilehit_tile)
        total_hits_counter += total_hits_frame

        #set counters
        corr_counter = 0
        uncorr_counter = 0
        
        #get primaries
        primaries_frame = get_mc_primary_for_hit_frame(ttree_mu3e)
        primaries_frame_arr = []
        for key in primaries_frame.keys():
            primaries_frame_arr.append([key,primaries_frame[key]]) #[hit tile, primary for tile hit]

        #get clusters
        clusters_with_primaries = sclump.build_cluster_with_truth_primary(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_primaries_arr = []
        cluster_master_primaries_arr = []
        for key in clusters_with_primaries.keys():
            cluster_master_primaries_arr.append(key)
            cluster_primaries_arr.append(clusters_with_primaries[key])

        #count hits in clusters
        cluster_hits_counter_tmp = 0
        for key in clusters_with_primaries.keys():
            cluster_hits_counter_tmp +=1
            for i in clusters_with_primaries[key]:
                cluster_hits_counter_tmp +=1
        cluster_hits_counter += cluster_hits_counter_tmp

        #comparison
        for j in range(len(cluster_primaries_arr)): #loop over all clusters in frame
            for k in range(len(cluster_primaries_arr[j])): #loop over all primaries in cluster 
                if cluster_primaries_arr[j][k] == cluster_master_primaries_arr[j]: #if primary in cluster = primary of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1 

        #add #master tiles to corr_counter
        corr_counter += len(cluster_master_primaries_arr)

        #add to total corr and uncorr counters
        tot_corr_counter += corr_counter
        tot_uncorr_counter += uncorr_counter

        if corr_counter != 0:
            frac_corr_clusters_frame.append(corr_counter/cluster_hits_counter_tmp)
            frac_corr_frame.append(corr_counter/total_hits_frame)

        if uncorr_counter != 0:
            frac_uncorr_frame.append(uncorr_counter/cluster_hits_counter_tmp)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    print("Total #hits in frames/#hits in clusters = ", total_hits_counter/cluster_hits_counter)
    print("Correctly associated out of all hits", tot_corr_counter/(total_hits_counter/100),"%")
    print("Incorrectly associated out of all hits", tot_uncorr_counter/(total_hits_counter/100),"%")
        
    return frac_corr_frame, frac_corr_clusters_frame, frac_uncorr_frame, tot_corr_counter

#----------------------------------------------------
#compares the tids of hits in cluster to the of cluster master. Returns the fractions (correctly associated hits)/(all hits in clusters)
# and (incorrectly associated hits)/(all hits in clusters)
def compare_to_tid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    frac_corr_frame = []
    frac_uncorr_frame = []
    frac_corr_clusters_frame = []
    total_hits_counter = 0
    cluster_hits_counter = 0

    tot_corr_counter = 0
    tot_uncorr_counter = 0

    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames

    for frame in range(frames_to_analyze):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')
        
        #count total hits
        total_hits_frame = len(ttree_mu3e.tilehit_tile)
        total_hits_counter += total_hits_frame

        #set counters
        corr_counter = 0
        uncorr_counter = 0
        
        #get primaries
        tids_frame = get_tid_frame(ttree_mu3e, ttree_mu3e_mc)
        tids_frame_arr = []
        for key in tids_frame.keys():
            tids_frame_arr.append([key,tids_frame[key]]) #[hit tile, tid for tile hit]

        #get clusters
        clusters_with_tids = sclump.build_cluster_with_truth_tid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_tids_arr = []
        cluster_master_tids_arr = []
        for key in clusters_with_tids.keys():
            cluster_master_tids_arr.append(key)
            cluster_tids_arr.append(clusters_with_tids[key])

        #count hits in clusters
        cluster_hits_counter_tmp = 0
        for key in clusters_with_tids.keys():
            cluster_hits_counter_tmp +=1
            for i in clusters_with_tids[key]:
                cluster_hits_counter_tmp +=1
        cluster_hits_counter += cluster_hits_counter_tmp

        #comparison
        for j in range(len(cluster_master_tids_arr)): #loop over all clusters in frame
            for k in range(len(cluster_tids_arr[j])): #loop over all tids in cluster 
                if cluster_tids_arr[j][k] == cluster_master_tids_arr[j]: #if tid in cluster = tid of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1  

        #add #master tiles to corr_counter
        corr_counter += len(cluster_master_tids_arr) 

        if (corr_counter + uncorr_counter) != cluster_hits_counter_tmp:
            print("error: counters don't match",(corr_counter + uncorr_counter), cluster_hits_counter_tmp)

        #add to total corr and uncorr counters
        tot_corr_counter += corr_counter
        tot_uncorr_counter += uncorr_counter

        #calculate fractions
        if corr_counter != 0:
            frac_corr_clusters_frame.append(corr_counter/cluster_hits_counter_tmp)
            frac_corr_frame.append(corr_counter/total_hits_frame)

        if uncorr_counter != 0:
            frac_uncorr_frame.append(uncorr_counter/cluster_hits_counter_tmp)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    print("Total #hits in frames/#hits in clusters = ", total_hits_counter/cluster_hits_counter)
    print("Correctly associated out of all hits", tot_corr_counter/(total_hits_counter/100),"%")
    print("Incorrectly associated out of all hits", tot_uncorr_counter/(total_hits_counter/100),"%")
        
    return frac_corr_frame, frac_corr_clusters_frame, frac_uncorr_frame, tot_corr_counter

#-----------------------------------------------------
#returns fraction of number of hits in cluster and total number of hits
def get_hits_not_in_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_hits_counter = []
    cluster_hits_counter = []
    frac_not_in_cluster = []

    #counting
    for frame in range(frames_to_analyze):
        ttree_mu3e.GetEntry(frame)

        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')
        
        #count total hits
        tot_hits_frame = len(ttree_mu3e.tilehit_tile)
        total_hits_counter.append(tot_hits_frame)

        #count hits in clusters
        clusters_frame = sclump.build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_hits_counter_tmp = 0
        for key in clusters_frame.keys():
            cluster_hits_counter_tmp +=1
            for i in clusters_frame[key]:
                cluster_hits_counter_tmp +=1
        cluster_hits_counter.append(cluster_hits_counter_tmp)

        #calculate fraction
        if cluster_hits_counter_tmp != 0:
            frac_not_in_cluster.append((tot_hits_frame-cluster_hits_counter_tmp)/tot_hits_frame)  


    print("Progress: 100 %","of ", frames_to_analyze, " frames")

    print("Not associated hits out of all hits: ",(np.sum(total_hits_counter)-np.sum(cluster_hits_counter))/(np.sum(total_hits_counter)/100) , "%")

    return frac_not_in_cluster

#-----------------------------------------------------
#returns fraction of number of hits in cluster and total number of hits
def get_hits_not_in_cluster_3_frame(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_hits_counter = []
    cluster_hits_counter = []
    frac_not_in_cluster = []
    #############################
    over_counter = 0
    ###############################

    #get total hits
    hits_all_frames ,hits_all_frames_counter_after = clump_3.del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, number_of_frames, rec_type)


    #counting
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)

        #Printing status info
        if frame % 2000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        tot_hits_frame = len(hits_all_frames[frame])
        total_hits_counter.append(tot_hits_frame)

        #count hits in clusters and
        #remove double hits in clusters like in check_for_mult_hit_tiles_diff_frame
        clusters_frame = clump_3.build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_hits_counter_tmp = 0
        double_hit_counter_tmp = 0
        for key in clusters_frame.keys():
            cluster = clusters_frame[key]
            cluster_hits_counter_tmp += len(cluster)
            for hit1 in cluster:
                for hit2 in cluster:
                    if hit1[0] == hit2[0] and hit1[1] != hit2[1]:
                        double_hit_counter_tmp += 1
            #correct for moved cluster hits
            for hit in cluster:
                if hit[1] != frame:
                    tot_hits_frame += 1
        cluster_hits_counter.append(cluster_hits_counter_tmp - double_hit_counter_tmp)
        ##################################
        if (cluster_hits_counter_tmp - double_hit_counter_tmp) > tot_hits_frame:
            over_counter += (cluster_hits_counter_tmp - double_hit_counter_tmp)- tot_hits_frame
        #    print(frame, (cluster_hits_counter_tmp - double_hit_counter_tmp)- tot_hits_frame)
        ################################

        #calculate fraction
        if tot_hits_frame != 0:
            frac_not_in_cluster.append((tot_hits_frame - (cluster_hits_counter_tmp - double_hit_counter_tmp))/tot_hits_frame)  

    if np.sum(total_hits_counter) != hits_all_frames_counter_after:
        print("ERROR: Total hit counters don't match", np.sum(total_hits_counter), hits_all_frames_counter_after)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")

    print("Not associated hits out of all hits: ", np.sum(cluster_hits_counter)/(np.sum(total_hits_counter)/100), "%")

    ##########################
    print(over_counter)
    ##########################

    return frac_not_in_cluster

#---------------------------------------------------------
def compare_to_primary_3_frames(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    frac_corr_frame = []
    frac_corr_clusters_frame = []
    frac_uncorr_frame = []
    total_hits_counter = []
    cluster_hits_counter = 0

    tot_corr_counter = 0
    tot_uncorr_counter = 0

    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames

    #get total hits
    hits_all_frames ,hits_all_frames_counter_after = clump_3.del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, number_of_frames, rec_type)


    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        tot_hits_frame = len(hits_all_frames[frame])
        total_hits_counter.append(tot_hits_frame)

        #set counters
        corr_counter = 0
        uncorr_counter = 0
        
        #get primaries
        primaries_frame_0 = get_mc_primary_for_hit_frame(ttree_mu3e)
        ttree_mu3e.GetEntry(frame+1)
        primaries_frame_plus = get_mc_primary_for_hit_frame(ttree_mu3e)
        ttree_mu3e.GetEntry(frame-1)
        primaries_frame_minus = get_mc_primary_for_hit_frame(ttree_mu3e)
        ttree_mu3e.GetEntry(frame)

        primaries_frame_arr_0 = []
        for key in primaries_frame_0.keys():
            primaries_frame_arr_0.append([key,primaries_frame_0[key]]) #[hit tile, primary for tile hit]

        primaries_frame_arr_plus = []
        for key in primaries_frame_plus.keys():
            primaries_frame_arr_plus.append([key,primaries_frame_plus[key]]) #[hit tile, primary for tile hit]

        primaries_frame_arr_minus = []
        for key in primaries_frame_minus.keys():
            primaries_frame_arr_minus.append([key,primaries_frame_minus[key]]) #[hit tile, primary for tile hit]

        #get clusters
        clusters_with_primaries = clump_3.build_cluster_with_truth_primary_3_frame(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_primaries_arr = []
        cluster_master_primaries_arr = []
        for key in clusters_with_primaries.keys():
            cluster_master_primaries_arr.append(key)
            cluster_primaries_arr.append(clusters_with_primaries[key])

        #count hits in clusters
        cluster_hits_counter_tmp = 0
        for key in clusters_with_primaries.keys():
            cluster_hits_counter_tmp += len(clusters_with_primaries[key])
        cluster_hits_counter += cluster_hits_counter_tmp

        #comparison
        for j in range(len(cluster_primaries_arr)): #loop over all clusters in frame
            for k in range(len(cluster_primaries_arr[j])): #loop over all primaries in cluster 
                if cluster_primaries_arr[j][k] == cluster_master_primaries_arr[j]: #if primary in cluster = primary of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1 

        #add #master tiles to corr_counter
        corr_counter += len(cluster_master_primaries_arr)

        #add to total corr and uncorr counters
        tot_corr_counter += corr_counter
        tot_uncorr_counter += uncorr_counter

        if cluster_hits_counter_tmp != 0:
            frac_corr_clusters_frame.append(corr_counter/cluster_hits_counter_tmp)

        if tot_hits_frame != 0:
            frac_corr_frame.append(corr_counter/tot_hits_frame)

        if cluster_hits_counter_tmp != 0:
            frac_uncorr_frame.append(uncorr_counter/cluster_hits_counter_tmp)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    print("Total #hits in frames/#hits in clusters = ", np.sum(total_hits_counter)/cluster_hits_counter)
    print("Correctly associated out of all hits", tot_corr_counter/(np.sum(total_hits_counter)/100),"%")
    print("Incorrectly associated out of all hits", tot_uncorr_counter/(np.sum(total_hits_counter)/100),"%")
        
    return frac_corr_frame, frac_corr_clusters_frame, frac_uncorr_frame, tot_corr_counter