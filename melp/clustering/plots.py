import ROOT
import numpy as np

import melp
import melp.clustering as clump

from melp.clustering.misc import*

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster

######################
#compares primary of the cluster hits to that of the cluster master and therefore provides efficiency data
def compare_to_primary(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, time_threshold, mask_type, number_of_frames = None, rec_type = None, cluster_type = None):
    frac_corr_frame          = []
    frac_corr_clusters_frame = []
    frac_uncorr_frame        = []
    total_hits_counter       = []
    cluster_hits_counter     = 0
    tot_corr_counter         = 0
    tot_uncorr_counter       = 0

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
        corr_counter   = 0
        uncorr_counter = 0

        #get clusters
        if cluster_type == "time":
            clusters = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)
        elif cluster_type == "timethenspatial":
            clusters = clump.three_dim_cluster.spatial_clustering_for_time_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, time_threshold, mask_type, rec_type)
        elif cluster_type == "timetheniterativespatial":
            clusters = clump.three_dim_cluster.iterative_masks_after_time_clustering(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, time_threshold, mask_type, rec_type)
        else:
            clusters = clump.spatial_cluster.build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        
        #----------------------
        #count hits in clusters
        #----------------------
        cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            cluster_hits_counter_tmp += clusters[i].__len__()
        cluster_hits_counter += cluster_hits_counter_tmp

        #--------------------------
        #comparison hits in cluster
        #--------------------------
        for j in range(len(clusters)): #loop over all clusters in frame
            cluster_primaries = clusters[j].get_primaries()
            for k in range(len(cluster_primaries)): #loop over all primaries in cluster 
                if cluster_primaries[k] == clusters[j].master_primary:#if primary in cluster = primary of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1

        #--------------------------------
        #comparison of different clusters
        #--------------------------------
        #define which cluster_types should be analyzed by this part of the algorithm
        sel_cluster_types = ["time", "timethenspatial", "timetheniterativespatial"]
        if cluster_type in sel_cluster_types:
            new_corr_cluster_flags = []
            old_corr_cluster_flags = []
            checked_primaries = []
            for i in range(len(clusters)):
                if clusters[i].master_primary not in checked_primaries:
                    number_of_primaries = 0
                    for hit in clusters[i].hits:
                        if hit.primary == clusters[i].master_primary:
                            number_of_primaries += 1
                    checked_primaries.append(clusters[i].master_primary)
                else:
                    continue
                for j in range(len(clusters)):
                    number_of_primaries_comp = 0
                    if j != i and j not in new_corr_cluster_flags:
                        for k in range(clusters[j].__len__()):
                            if clusters[j].hits[k].primary == clusters[i].master_primary:
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
                master_primary = clusters[i].master_primary
                if master_primary not in checked_primaries_2:
                    number_of_primaries = 0
                    for hit in clusters[i].hits:
                        if hit.primary == master_primary:
                            number_of_primaries += 1
                    checked_primaries_2.append(master_primary)
                else:
                    continue
                for j in range(len(clusters)):
                    number_of_primaries_comp = 0
                    if j != i and j not in new_corr_cluster_flags:
                        for k in range(len(clusters[j])):
                            if clusters[j].hits[k].primary == master_primary:
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

        #-------------------------------------
        #add to total corr and uncorr counters
        #-------------------------------------
        tot_corr_counter   += corr_counter
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

######################
#compares tid of the cluster hits to that of the cluster master and therefore provides efficiency data
#@blockPrinting
def compare_to_tid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, time_threshold, threshold_cluster_width, mask_type, number_of_frames = None, rec_type = None, cluster_type = None):
    frac_corr_frame            = []
    frac_corr_clusters_frame   = []
    frac_uncorr_frame          = []
    total_hits_counter         = []
    number_of_tids_0             = []
    cluster_hits_counter       = 0
    tot_corr_counter           = 0
    tot_uncorr_counter         = 0

    tot_cluster_counter        = 0
    double_tid_cluster_counter = 0
    double_tid_cluster_hits_counter = 0
    corr_double_tid_cluster_counter = 0
    long_time_between_cluster_hits_counter = 0

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
        corr_counter   = 0
        uncorr_counter = 0

        #get clusters
        if cluster_type == "time":
            clusters = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)
        elif cluster_type == "timethenspatial":
            clusters = clump.three_dim_cluster.spatial_clustering_for_time_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, time_threshold, mask_type, rec_type)
        elif cluster_type == "timetheniterativespatial":
            clusters = clump.three_dim_cluster.iterative_masks_after_time_clustering(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, time_threshold, mask_type, rec_type)
        elif cluster_type == "iterativespatial":
            clusters = clump.spatial_cluster.iterative_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        elif cluster_type == "truth":
            clusters = clump.spatial_cluster.build_truth_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        else:
            clusters = clump.spatial_cluster.build_clusters_in_masks(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        
        #----------------------
        #count hits in clusters
        #----------------------
        cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            tot_cluster_counter += 1
            cluster_hits_counter_tmp += len(clusters[i])
        cluster_hits_counter += cluster_hits_counter_tmp

        #---------------------------------------------------
        #count clusters with hits that are far apart in time
        #---------------------------------------------------
        long_time_between_cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            times = clusters[i].get_times()
            if len(times) > 0:
                min_time = min(times)
                max_time = max(times)
                if max_time - min_time > 0.5:
                    long_time_between_cluster_hits_counter_tmp +=1
        long_time_between_cluster_hits_counter += long_time_between_cluster_hits_counter_tmp

        #-----------------------------------------
        #Count number of different tids in cluster
        #-----------------------------------------
        number_of_tids_tmp_0 = []
        for i in range(len(clusters)):
            tids = clusters[i].get_tids()
            tids_checked = []
            diff_tid_counter = 0
            for tid in tids:
                if tid not in tids_checked:
                    diff_tid_counter += 1
                    tids_checked.append(tid)
            number_of_tids_tmp_0.append(diff_tid_counter)
        number_of_tids_0.extend(number_of_tids_tmp_0)

        #--------------------------
        #comparison hits in cluster
        #--------------------------
        for j in range(len(clusters)): #loop over all clusters in frame
            cluster_tids = clusters[j].get_tids()
            for k in range(len(cluster_tids)): #loop over all tids in cluster 
                if cluster_tids[k] == clusters[j].master_tid:#if tid in cluster = tid of cluster master
                    corr_counter += 1
                else:
                    uncorr_counter += 1
        
        #--------------------------------
        #comparison of different clusters
        #--------------------------------
        #define which cluster_types should be analyzed by this part of the algorithm
        sel_cluster_types = ["time", "timethenspatial", "timetheniterativespatial", "iterativespatial"]
        if cluster_type in sel_cluster_types:
            #count number of clusters that have same master_tid as other cluster
            for i in range(len(clusters)):
                master_tid_1 = clusters[i].master_tid
                for j in range(len(clusters)):
                    master_tid_2 = clusters[j].master_tid
                    if j != i and master_tid_1 == master_tid_2:
                        double_tid_cluster_counter += 1
                        double_tid_cluster_hits_counter += len(clusters[j])

            #correct for clusters with multiple tids
            for i in range(len(clusters)):
                master_tid_1 = clusters[i].master_tid
                double_tid_clusters_tmp = []
                min_times_clusters = [clusters[i]]
                for j in range(len(clusters)):
                    master_tid_2 = clusters[j].master_tid
                    if j != i and master_tid_1 == master_tid_2:
                        corr_double_tid_cluster_counter += 1
                        double_tid_clusters_tmp.append(clusters[j])
                        min_times_clusters.append(clusters[j])
                if len(double_tid_clusters_tmp) != 0:
                    min_times_tmp = []
                    for cluster_tmp in min_times_clusters:
                        min_times_tmp.append(min(cluster_tmp.get_times()))
                    index_first_cluster = min_times_tmp.index(min(min_times_tmp))
                    for k in range(len(min_times_clusters)):
                        if k != index_first_cluster:
                            pos_first_cluster = mu3e_detector.TileDetector.tile[min_times_clusters[index_first_cluster].hits[0].tile_id].pos
                            pos_double_cluster = mu3e_detector.TileDetector.tile[min_times_clusters[k].hits[0].tile_id].pos
                            distance = np.sqrt((pos_first_cluster[0] - pos_double_cluster[0]) ** 2 + (pos_first_cluster[1] - pos_double_cluster[1]) ** 2 + (pos_first_cluster[2] - pos_double_cluster[2]) ** 2) #mm
                            #applying size threshold
                            if cluster_type == "time":
                                #count hits in clusters that don't contain first hit with same tid
                                same_tid_counter_tmp = 0
                                for hit in min_times_clusters[k].hits:
                                    if hit.tid == master_tid_1:
                                        corr_counter   -= 1
                                        uncorr_counter += 1
                            else:
                                if distance < threshold_cluster_width: #mm
                                    #count hits in clusters that don't contain first hit with same tid
                                    same_tid_counter_tmp = 0
                                    for hit in min_times_clusters[k].hits:
                                        if hit.tid == master_tid_1:
                                            corr_counter   -= 1
                                            uncorr_counter += 1
                                    #corr_counter   -= len(min_times_clusters[k])
                                    #uncorr_counter += len(min_times_clusters[k])


        """
        #--------------------------------------------------------
        #comparison of different clusters from primary comparison
        #--------------------------------------------------------
        #define which cluster_types should be analyzed by this part of the algorithm
        sel_cluster_types = ["time", "timethenspatial", "timetheniterativespatial"]
        if cluster_type in sel_cluster_types:
            new_corr_cluster_flags = []
            old_corr_cluster_flags = []
            checked_tids = []
            for i in range(len(clusters)):
                if clusters[i].master_tid not in checked_tids:
                    number_of_tids = 0
                    for hit in clusters[i].hits:
                        if hit.tid == clusters[i].master_tid:
                            number_of_tids += 1
                    checked_tids.append(clusters[i].master_tid)
                else:
                    continue
                for j in range(len(clusters)):
                    number_of_tids_comp = 0
                    if j != i and j not in new_corr_cluster_flags:
                        for k in range(len(clusters[j])):
                            if clusters[j].hits[k].tid == clusters[i].master_tid:
                                number_of_tids_comp += 1
                        if number_of_tids_comp == 0: #if master tid of cluster i isn't found in cluster j do nothing
                            continue
                        elif number_of_tids_comp <= number_of_tids: #if correctly identified constituents are more in cluster i simply add cluster j as wrongly identified
                            #TODO: maybe split into < and = and decide for the correct cluster either via the smallest timestamp or by amount of wrong hits in cluster
                            corr_counter -= number_of_tids_comp
                            uncorr_counter += number_of_tids_comp
                        elif number_of_tids_comp > number_of_tids: #if cluster j has more correct primaries flag it as correct cluster and add cluster i to the incorrect counter
                            corr_counter -= number_of_tids
                            uncorr_counter += number_of_tids
                            new_corr_cluster_flags.append(j)
                            old_corr_cluster_flags.append(i)
            
                       
            #loop over old correct cluster flags
            checked_tids_2 = []
            old_corr_cluster_flags_check = []
            for i in old_corr_cluster_flags:
                master_tid = clusters[i].master_tid
                if master_tid not in checked_tids_2:
                    number_of_tids = 0
                    for hit in clusters[i].hits:
                        if hit.tid == master_tid:
                            number_of_tids += 1
                    checked_tids_2.append(master_tid)
                else:
                    continue
                for j in range(len(clusters)):
                    number_of_tids_comp = 0
                    if j != i and j not in new_corr_cluster_flags:
                        for k in range(len(clusters[j])):
                            if clusters[j].hits[k].tid == master_tid:
                                number_of_tids_comp += 1
                        if number_of_tids_comp == 0: #if master tid of cluster i isn't found in cluster j do nothing
                            continue
                        elif number_of_tids_comp <= number_of_tids: #if correctly identified constituents are more in cluster i simply add cluster j as wrongly identified
                            #TODO: maybe split into < and = and decide for the correct cluster either via the smallest timestamp or by amount of wrong hits in cluster
                            corr_counter -= number_of_tids_comp
                            uncorr_counter += number_of_tids_comp
                        elif number_of_tids_comp > number_of_tids: #if cluster j has more correct primaries flag it as correct cluster and add cluster i to the incorrect counter
                            corr_counter -= number_of_tids
                            uncorr_counter += number_of_tids
                            new_corr_cluster_flags.append(j)
                            old_corr_cluster_flags.append(i)
                            old_corr_cluster_flags_check.append(i)
        """
        #-------------------------------------
        #add to total corr and uncorr counters
        #-------------------------------------
        tot_corr_counter   += corr_counter
        tot_uncorr_counter += uncorr_counter

        if cluster_hits_counter_tmp != 0:
            frac_corr_clusters_frame.append(corr_counter/cluster_hits_counter_tmp)
            frac_uncorr_frame.append(uncorr_counter/cluster_hits_counter_tmp)

        if total_hits_frame != 0:
            frac_corr_frame.append(corr_counter/total_hits_frame)
            

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    if cluster_hits_counter > 0:
        print("Number of analyzed frames: ", len(total_hits_counter), "Number of correct counter fractions: ", len(frac_corr_frame))
        print("Total number of hits =",np.sum(total_hits_counter), ", Identified correctly + identified incorrectly =", tot_corr_counter + tot_uncorr_counter)
        print("Identified correctly:", tot_corr_counter)
        print("Identified incorrectly:", tot_uncorr_counter)
        print("Total #hits in frames/#hits in clusters = ", np.sum(total_hits_counter)/cluster_hits_counter)
        print("Total number of clusters:",tot_cluster_counter, ", Hits:",np.sum(total_hits_counter))
        print("Number of clusters with hits that are far apart in time:", long_time_between_cluster_hits_counter)
        print("Number of clusters where tid already exists:",double_tid_cluster_counter, ", Hits:", double_tid_cluster_hits_counter)
        print("Number of clusters where tid already exists, that are accounted for:", corr_double_tid_cluster_counter)
        print("Correctly associated out of all hits: ", tot_corr_counter/(np.sum(total_hits_counter)/100),"%")
        print("Correctly associated out of all hits in clusters: ", tot_corr_counter/(cluster_hits_counter/100),"%")
        print("Incorrectly associated out of all hits: ", tot_uncorr_counter/(np.sum(total_hits_counter)/100),"%")
        print("Incorrectly associated out of all hits in clusters: ", tot_uncorr_counter/(cluster_hits_counter/100),"%")
   
    return frac_corr_frame, frac_corr_clusters_frame, frac_uncorr_frame, tot_corr_counter, total_hits_counter, number_of_tids_0


###########################
#returns fraction of number of hits in cluster and total number of hits
def get_hits_not_in_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, time_threshold, mask_type, number_of_frames = None, rec_type = None, cluster_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_hits_counter   = []
    cluster_hits_counter = []
    frac_not_in_cluster  = []

    #--------
    #counting
    #--------
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
            clusters_frame = clump.time_cluster.time_clustering_frame_improv(ttree_mu3e, ttree_mu3e_mc, frame, time_threshold)
        elif cluster_type == "timethenspatial":
            clusters_frame = clump.three_dim_cluster.spatial_clustering_for_time_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, time_threshold, mask_type, rec_type)
        elif cluster_type == "timetheniterativespatial":
            clusters_frame = clump.three_dim_cluster.iterative_masks_after_time_clustering(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, time_threshold, mask_type, rec_type)
        
    

        #count hits in clusters
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

###################################
#returns fraction of number of hits in cluster and total number of hits
def get_hits_not_in_cluster_3_frame(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_hits_counter   = []
    cluster_hits_counter = []
    frac_not_in_cluster  = []
    over_counter         = 0

    #--------------
    #get total hits
    #--------------
    hits_all_frames ,hits_all_frames_counter_after = clump.three_frame_cluster.del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, number_of_frames, rec_type)

    #--------
    #counting
    #--------
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)

        #Printing status info
        if frame % 2000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        tot_hits_frame = len(hits_all_frames[frame])
        total_hits_counter.append(tot_hits_frame)

        #--------------------------------------------------------------------------
        #count hits in clusters and
        #remove double hits in clusters like in check_for_mult_hit_tiles_diff_frame
        #--------------------------------------------------------------------------
        clusters_frame = clump.three_frame_cluster.build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)
        cluster_hits_counter_tmp = 0
        double_hit_counter_tmp   = 0

        for i in range(len(clusters_frame)):
            cluster_hits_counter_tmp += clusters_frame[i].__len__()
        
        for cluster in clusters_frame:
            for hit1 in cluster.hits:
                for hit2 in cluster.hits:
                    if hit1.tile_id == hit2.tile_id and hit1.frame_id != hit2.frame_id:
                        double_hit_counter_tmp += 1
            #correct for moved cluster hits
            for hit in cluster.hits:
                if hit.frame_id != frame:
                    tot_hits_frame += 1
        cluster_hits_counter.append(cluster_hits_counter_tmp - double_hit_counter_tmp)

        if tot_hits_frame != 0 and (cluster_hits_counter_tmp - double_hit_counter_tmp) > tot_hits_frame:
            over_counter += (cluster_hits_counter_tmp - double_hit_counter_tmp)- tot_hits_frame

        #calculate fraction
        if tot_hits_frame != 0:
            frac_not_in_cluster.append((tot_hits_frame - (cluster_hits_counter_tmp - double_hit_counter_tmp))/tot_hits_frame)  

    if np.sum(total_hits_counter) != hits_all_frames_counter_after:
        print("ERROR: Total hit counters don't match", np.sum(total_hits_counter), hits_all_frames_counter_after)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")

    print("Not associated hits out of all hits: ", np.sum(cluster_hits_counter)/(np.sum(total_hits_counter)/100), "%")

    print("Hits missed by the code: ", over_counter)

    return frac_not_in_cluster

###############################
#compares primary of the cluster hits to that of the cluster master and therefore provides efficiency data
def compare_to_primary_3_frames(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    frac_corr_frame          = []
    frac_corr_clusters_frame = []
    frac_uncorr_frame        = []
    total_hits_counter       = []
    cluster_hits_counter     = 0

    tot_corr_counter         = 0
    tot_uncorr_counter       = 0

    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames

    #--------------
    #get total hits
    #--------------
    hits_all_frames, __ = clump.three_frame_cluster.del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, mask_type, number_of_frames, rec_type)


    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #count total hits
        tot_hits_frame = len(hits_all_frames[frame])
        total_hits_counter.append(tot_hits_frame)

        #set counters
        corr_counter   = 0
        uncorr_counter = 0
        
        #-------------
        #get primaries
        #-------------
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

        #------------
        #get clusters
        #------------
        clusters = clump.three_frame_cluster.build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)


        cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            cluster_hits_counter_tmp += clusters[i].__len__()
        cluster_hits_counter += cluster_hits_counter_tmp

        #----------
        #comparison
        #----------
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


###########################################
#Gives the efficiency as function of maximum cluster width (clusters with same tid that have been identified as seperate clusters aren't counted as wrongly identified if they are further apart then the given threshold)
def efficiency_as_function_of_cluster_width(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, time_threshold, mask_type, number_of_frames = None, rec_type = None, cluster_type = None):
    threshold_cluster_width_array = np.arange(10,100,2)
    efficiency = []
    finished = 0
    for threshold_cluster_width in threshold_cluster_width_array:
        #print status info
        print("Progress: ", finished,"of ", len(threshold_cluster_width_array), " thresholds at", number_of_frames, "frames each", end='\r')

        #get efficiency for the threshold
        with HiddenPrints():
            __, __, __, tot_corr_counter, total_hits_counter, __ = compare_to_tid(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, time_threshold, threshold_cluster_width, mask_type, number_of_frames, rec_type, cluster_type)
        efficiency.append(tot_corr_counter/(np.sum(total_hits_counter)/100))
        finished += 1

    print("Progress: ", len(threshold_cluster_width_array),"of ", len(threshold_cluster_width_array), " thresholds at", number_of_frames, "frames each")
    
    return efficiency, threshold_cluster_width_array


