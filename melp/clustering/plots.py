import ROOT
import numpy as np

import melp
import melp.clustering as clump

from melp.clustering.misc import*

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster

#-------------------------------------------------
def compare_to_primary(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None, cluster_type = None):
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

    for frame in range(frames_to_analyze):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        total_hits_frame = ttree_mu3e.Ntilehit 
        total_hits_counter.append(total_hits_frame)

        #set counters
        corr_counter = 0
        uncorr_counter = 0
        
        #get primaries
        #primaries_frame = get_mc_primary_for_hit_frame(ttree_mu3e)
        #primaries_frame_arr = []
        #for key in primaries_frame.keys():
        #    primaries_frame_arr.append([key,primaries_frame[key]]) #[hit tile, primary for tile hit]

        #get clusters
        #clusters_with_primaries = clump.spatial_cluster.build_cluster_with_truth_primary(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type, cluster_type)
        #cluster_primaries_arr = []
        #cluster_master_primaries_arr = []
        #for key in clusters_with_primaries.keys():
        #    cluster_master_primaries_arr.append(clusters_with_primaries[key][0])
        #    cluster_primaries_arr.append(clusters_with_primaries[key])
        clusters = clump.spatial_cluster.build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        
        #count hits in clusters
        #cluster_hits_counter_tmp = 0
        #for key in clusters_with_primaries.keys():
        #    cluster_hits_counter_tmp += len(clusters_with_primaries[key])
        #cluster_hits_counter += cluster_hits_counter_tmp
        cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            cluster_hits_counter_tmp += clusters[i].__len__()
        cluster_hits_counter += cluster_hits_counter_tmp

        #comparison hits in cluster
        #for j in range(len(cluster_primaries_arr)): #loop over all clusters in frame
        #    for k in range(len(cluster_primaries_arr[j])): #loop over all primaries in cluster 
        #        if cluster_primaries_arr[j][k] == cluster_master_primaries_arr[j]: #if primary in cluster = primary of cluster master
        #            corr_counter += 1
        #        else:
        #            uncorr_counter += 1 
        for j in range(len(clusters)): #loop over all clusters in frame
            cluster_primaries = clusters[j].get_primaries()
            for k in range(len(cluster_primaries)): #loop over all primaries in cluster 
                if cluster_primaries[k] == clusters[j].master_primary:#if primary in cluster = primary of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1

            ################################
            #if frame < 10:
            #    print("master primary: ", clusters[j].master_primary, "\n")
            #    print("Cluster primaries: ", cluster_primaries, "\n")
            ################################
        

        #comparison of different clusters for clustering methods without master
        if cluster_type == "time":
            new_corr_cluster_flags = []
            old_corr_cluster_flags = []
            checked_primaries = []
            #for i in range(len(cluster_primaries_arr)):
            for i in range(len(clusters)): #loop over all clusters in frame
                #master_primary = cluster_primaries_arr[i][0]
                #if master_primary not in checked_primaries:
                if clusters[i].master_primary not in checked_primaries:
                    cluster_primaries = clusters[i].get_primaries()
                    #number_of_primaries = cluster_primaries_arr[i].count(master_primary)
                    number_of_primaries = cluster_primaries.count(clusters[i].master_primary)
                    #checked_primaries.append(master_primary)
                    checked_primaries.append(clusters[i].master_primary)
                else:
                    continue
                #for j in range(len(cluster_primaries_arr)):
                for j in range(len(clusters)): #loop over all clusters in frame
                    number_of_primaries_comp = 0
                    cluster_primaries = clusters[j].get_primaries()
                    if j != i and j not in new_corr_cluster_flags:
                        #for k in range(len(cluster_primaries_arr[j])):
                        for k in range(clusters[j].__len__()):
                            #if cluster_primaries_arr[j][k] == master_primary:
                            if cluster_primaries[k] == clusters[j].master_primary:
                                number_of_primaries_comp += 1
                        if number_of_primaries_comp == 0: #if master primary of cluster i isn't found in cluster j do nothing
                            continue
                        elif number_of_primaries_comp <= number_of_primaries: #if correctly identified constituents are more in cluster i simply add cluster j as wrongly identified
                            #TODO: maybe split into < and = and decide for the correct cluster either via the smallest timestamp or by amount of wrong hits in cluster
                            corr_counter -= number_of_primaries_comp
                            uncorr_counter += number_of_primaries_comp
                        elif number_of_primaries_comp > number_of_primaries: #if cluster j has more correct primaries flag it as correct cluster and add cluster i to the incorrect counter
                            corr_counter -= number_of_primaries
                            uncorr_counter += number_of_primaries
                            new_corr_cluster_flags.append(j)
                            old_corr_cluster_flags.append(i)
            
                        
            #loop over old correct cluster flags
            checked_primaries_2 = []
            old_corr_cluster_flags_check = []
            for i in old_corr_cluster_flags:
                #master_primary = cluster_primaries_arr[i][0]
                #if master_primary not in checked_primaries_2:
                if clusters[i].master_primary not in checked_primaries_2:
                    cluster_primaries = clusters[i].get_primaries()
                    #number_of_primaries = cluster_primaries_arr[i].count(master_primary)
                    number_of_primaries = cluster_primaries.count(clusters[i].master_primary)
                    checked_primaries_2.append(clusters[i].master_primary)
                else:
                    continue
                #for j in range(len(cluster_primaries_arr)):
                for j in range(len(clusters)):
                    number_of_primaries_comp = 0
                    cluster_primaries = clusters[j].get_primaries()
                    if j != i and j not in new_corr_cluster_flags:
                        #for k in range(len(cluster_primaries_arr[j])):
                        for k in range(clusters[j].__len__()):
                            #if cluster_primaries_arr[j][k] == master_primary:
                            if cluster_primaries[k] == clusters[j].master_primary:
                                number_of_primaries_comp += 1
                        if number_of_primaries_comp == 0: #if master primary of cluster i isn't found in cluster j do nothing
                            continue
                        elif number_of_primaries_comp <= number_of_primaries: #if correctly identified constituents are more in cluster i simply add cluster j as wrongly identified
                            #TODO: maybe split into < and = and decide for the correct cluster either via the smallest timestamp or by amount of wrong hits in cluster
                            corr_counter -= number_of_primaries_comp
                            uncorr_counter += number_of_primaries_comp
                        elif number_of_primaries_comp > number_of_primaries: #if cluster j has more correct primaries flag it as correct cluster and add cluster i to the incorrect counter
                            corr_counter -= number_of_primaries
                            uncorr_counter += number_of_primaries
                            new_corr_cluster_flags.append(j)
                            old_corr_cluster_flags.append(i)
                            old_corr_cluster_flags_check.append(i)
        

        #add to total corr and uncorr counters
        tot_corr_counter += corr_counter
        tot_uncorr_counter += uncorr_counter

        if cluster_hits_counter_tmp != 0:
            frac_corr_clusters_frame.append(corr_counter/cluster_hits_counter_tmp)
            frac_uncorr_frame.append(uncorr_counter/cluster_hits_counter_tmp)

        if total_hits_frame != 0:
            frac_corr_frame.append(corr_counter/total_hits_frame)
            

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    print("Number of analyzed frames: ", len(total_hits_counter), "Number of correct counter fractions: ", len(frac_corr_frame))
    print("Total #hits in frames/#hits in clusters = ", np.sum(total_hits_counter)/cluster_hits_counter)
    print("Correctly associated out of all hits: ", tot_corr_counter/(np.sum(total_hits_counter)/100),"%")
    print("Correctly associated out of all hits in clusters: ", tot_corr_counter/(cluster_hits_counter/100),"%")
    print("Incorrectly associated out of all hits: ", tot_uncorr_counter/(np.sum(total_hits_counter)/100),"%")
    print("Incorrectly associated out of all hits in clusters: ", tot_uncorr_counter/(cluster_hits_counter/100),"%")
   
    return frac_corr_frame, frac_corr_clusters_frame, frac_uncorr_frame, tot_corr_counter


#-----------------------------------------------------
#returns fraction of number of hits in cluster and total number of hits
def get_hits_not_in_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None, cluster_type = None):
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
        tot_hits_frame = ttree_mu3e.Ntilehit 
        total_hits_counter.append(tot_hits_frame)

        #get clusters
        if cluster_type == None:
            clusters_frame = clump.spatial_cluster.build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector,frame, mask_type, rec_type)
        elif cluster_type == "time":
            __ , clusters_frame = clump.time_cluster.time_clustering_frame(ttree_mu3e, printing = None)

        #count hits in clusters
        #cluster_hits_counter_tmp = 0
        #for key in clusters_frame.keys():
        #    cluster_hits_counter_tmp += len(clusters_frame[key])
        #cluster_hits_counter.append(cluster_hits_counter_tmp)

        cluster_hits_counter_tmp = 0
        for i in range(len(clusters_frame)):
            cluster_hits_counter_tmp += clusters_frame[i].__len__()
        cluster_hits_counter.append(cluster_hits_counter_tmp)

        #calculate fraction
        if tot_hits_frame != 0:
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
    hits_all_frames ,hits_all_frames_counter_after = clump.three_frame_cluster.del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, number_of_frames, rec_type)


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
        clusters_frame = clump.three_frame_cluster.build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_hits_counter_tmp = 0
        double_hit_counter_tmp = 0

        for i in range(len(clusters_frame)):
            cluster_hits_counter_tmp += clusters_frame[i].__len__()
        
        #for key in clusters_frame.keys():
        for cluster in clusters_frame:
            #cluster_hits_counter_tmp += len(cluster.hits)
            for hit1 in cluster.hits:
                for hit2 in cluster.hits:
                    #if hit1[0] == hit2[0] and hit1[1] != hit2[1]:
                    if hit1.tile_id == hit2.tile_id and hit1.frame_id != hit2.frame_id:
                        double_hit_counter_tmp += 1
            #correct for moved cluster hits
            for hit in cluster.hits:
                #if hit[1] != frame:
                if hit.frame_id != frame:
                    tot_hits_frame += 1
        cluster_hits_counter.append(cluster_hits_counter_tmp - double_hit_counter_tmp)
        ##################################
        if tot_hits_frame != 0 and (cluster_hits_counter_tmp - double_hit_counter_tmp) > tot_hits_frame:
            over_counter += (cluster_hits_counter_tmp - double_hit_counter_tmp)- tot_hits_frame
            #print((cluster_hits_counter_tmp - double_hit_counter_tmp)- tot_hits_frame)
        ################################

        #calculate fraction
        if tot_hits_frame != 0:
            frac_not_in_cluster.append((tot_hits_frame - (cluster_hits_counter_tmp - double_hit_counter_tmp))/tot_hits_frame)  

    if np.sum(total_hits_counter) != hits_all_frames_counter_after:
        print("ERROR: Total hit counters don't match", np.sum(total_hits_counter), hits_all_frames_counter_after)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")

    print("Not associated hits out of all hits: ", np.sum(cluster_hits_counter)/(np.sum(total_hits_counter)/100), "%")

    ##########################
    print("Hits missed by the code: ", over_counter)
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
    hits_all_frames ,hits_all_frames_counter_after = clump.three_frame_cluster.del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, number_of_frames, rec_type)


    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        tot_hits_frame = ttree_mu3e.Ntilehit #len(hits_all_frames[frame])
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
        #clusters_with_primaries = clump.three_frame_cluster.build_cluster_with_truth_primary_3_frame(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        clusters = clump.three_frame_cluster.build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        #cluster_primaries_arr = []
        #cluster_master_primaries_arr = []
        #for key in clusters_with_primaries.keys():
        #for cluster in clusters:
        #    cluster_master_primaries_arr.append(key)
        #    cluster_primaries_arr.append(clusters_with_primaries[key])

        #count hits in clusters
        #cluster_hits_counter_tmp = 0
        #for key in clusters_with_primaries.keys():
        #for cluster in clusters:
            #cluster_hits_counter_tmp += len(clusters_with_primaries[key])
        #    cluster_hits_counter_tmp += cluster.__len__()
        #cluster_hits_counter += cluster_hits_counter_tmp

        cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            cluster_hits_counter_tmp += clusters[i].__len__()
        cluster_hits_counter += cluster_hits_counter_tmp

        #comparison
        #for j in range(len(cluster_primaries_arr)): #loop over all clusters in frame
        #    for k in range(len(cluster_primaries_arr[j])): #loop over all primaries in cluster 
        #        if cluster_primaries_arr[j][k] == cluster_master_primaries_arr[j]: #if primary in cluster = primary of cluster master
        #            corr_counter += 1
        #        else:
        #            uncorr_counter += 1 

        #comparison
        for j in range(len(clusters)): #loop over all clusters in frame
            cluster_primaries = clusters[j].get_primaries()
            for k in range(len(cluster_primaries)): #loop over all primaries in cluster 
                if cluster_primaries[k] == clusters[j].master_primary:#if primary in cluster = primary of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1

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
    print("Correctly associated out of all hits in clusters", tot_corr_counter/(cluster_hits_counter/100),"%")
    print("Incorrectly associated out of all hits", tot_uncorr_counter/(np.sum(total_hits_counter)/100),"%")
    print("Incorrectly associated out of all hits in clusters", tot_uncorr_counter/(cluster_hits_counter/100),"%")
        
    return frac_corr_frame, frac_corr_clusters_frame, frac_uncorr_frame, tot_corr_counter