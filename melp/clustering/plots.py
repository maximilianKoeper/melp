import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

import melp.clustering.spatial_cluster as sclump

#-------------------------------------------------
def compare_to_primary(ttree_mu3e, ttree_mu3e_mc, mu3e_detector: melp.Detector, mask_type):
    #comparing to primary
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