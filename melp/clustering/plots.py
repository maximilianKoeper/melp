import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

import melp.clustering.spatial_cluster as sclump

#-------------------------------------------------
def compare_to_primary(ttree_mu3e, ttree_mu3e_mc, mu3e_detector: melp.Detector, mask_type):
    frac_corr_frame = []
    frac_uncorr_frame = []
    total_hits = []
    #primary_found_counter = 0
    #primary_not_found_counter = 0

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        #if frame % 1000 == 0:
        #    print("Progress: ", np.round(frame / ttree_mu3e.GetEntries() * 100), " %", end='\r')

        
        tot_hits_frame = len(ttree_mu3e.tilehit_tile)
        corr_counter = 0
        uncorr_counter = 0
        
        #get primaries
        primaries_frame = get_mc_primary_for_hit_frame(ttree_mu3e)
        #primaries_frame = get_mc_primary_for_hit_frame(filename_sorted, frame)
        primaries_frame_arr = []
        for key in primaries_frame.keys():
            primaries_frame_arr.append([key,primaries_frame[key]]) #[hit tile, primary for tile hit]

        #get clusters
        clusters_with_primaries = sclump.build_cluster_with_truth_primary(ttree_mu3e, ttree_mu3e_mc, mu3e_detector, mask_type)
        #clusters_with_primaries = sclump.build_cluster_with_truth_primary(filename_sorted, frame, mu3e_detector, mask_type = "big")
        cluster_primaries_arr = []
        cluster_arr = []
        for key in clusters_with_primaries.keys():
            cluster_primaries_arr.append(key)
            cluster_arr.append(clusters_with_primaries[key])

        #comparison
        for i in range(len(primaries_frame_arr)):
            for j in range(len(cluster_primaries_arr)):
                if cluster_primaries_arr[j] == primaries_frame_arr[i][1]:
                    #primary_found_counter += 1
                    for k in range(len(cluster_arr[j])):
                        if cluster_arr[j][k] == primaries_frame_arr[i][0]:
                            corr_counter += 1
                        else:
                            uncorr_counter += 1
                #else:
                    #primary_not_found_counter += 1

        if corr_counter != 0:
            frac_corr_frame.append(corr_counter/tot_hits_frame)  

        if uncorr_counter != 0:
            frac_uncorr_frame.append(uncorr_counter/tot_hits_frame)

        total_hits.append(tot_hits_frame)
        
    return frac_corr_frame, frac_uncorr_frame

#----------------------------------------------------
#compares the tids of hits in cluster to the of cluster master. Returns the fractions (correctly associated hits)/(all hits in clusters)
# and (incorrectly associated hits)/(all hits in clusters)
def compare_to_tid(ttree_mu3e, ttree_mu3e_mc, mu3e_detector: melp.Detector, mask_type, number_of_frames = None):
    frac_corr_frame = []
    frac_uncorr_frame = []
    total_hits_counter = 0
    cluster_hits_counter = 0

    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames

    for frame in range(frames_to_analyze):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 2000 == 0:
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
        clusters_with_tids = sclump.build_cluster_with_truth_tid(ttree_mu3e, ttree_mu3e_mc, mu3e_detector, mask_type)
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
        for j in range(len(cluster_tids_arr)): #loop over all clusters in frame
            #if cluster_master_tids_arr[j] == tids_frame_arr[i][1]: #if tid of cluster master = a tid in frame
            for k in range(len(cluster_tids_arr[j])): #loop over all tids in cluster 
                if cluster_tids_arr[j][k] == cluster_master_tids_arr[j]: #if tid in cluster = tid of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1            

        #calculate fractions
        if corr_counter != 0:
            frac_corr_frame.append(corr_counter/cluster_hits_counter_tmp)  

        if uncorr_counter != 0:
            frac_uncorr_frame.append(uncorr_counter/cluster_hits_counter_tmp)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
        
    return frac_corr_frame, frac_uncorr_frame

