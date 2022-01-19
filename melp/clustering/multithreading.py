# import modules
import melp
import multiprocessing as mp
import subprocess
from glob import glob
from functools import partial
import numpy as np
import time
import ROOT
from melp import Detector

from melp.clustering.misc import*

from melp import clustering as clump

from melp.src.cluster import ClusterHit
from melp.src.cluster import Cluster

import os, sys

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


#defining the functions with multithreading support
###############################
def compare_to_primary_filename(filename, time_threshold, mask_type, rec_type = None, cluster_type = None):
    #get file
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e_mc = file.Get("mu3e_mchits")
    ttree_sensor = file.Get("alignment/sensors")
    ttree_tiles = file.Get("alignment/tiles")

    with HiddenPrints():
        mu3e_detector = Detector.initFromROOT(filename)

    #define counters and arrays
    frac_corr_frame          = []
    frac_corr_clusters_frame = []
    frac_uncorr_frame        = []
    total_hits_counter       = []
    cluster_hits_counter     = 0
    tot_corr_counter         = 0
    tot_uncorr_counter       = 0

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        #count total hits
        total_hits_frame = ttree_mu3e.Ntilehit 
        total_hits_counter.append(total_hits_frame)

        #set counters
        corr_counter   = 0
        uncorr_counter = 0

        #get clusters
        if cluster_type == "time":
            clusters = clump.time_cluster.time_clustering_frame(ttree_mu3e, frame, printing = None)
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

    #write overall stats to txt file 
    f = open("./melp/clustering/results/efficiency_stats_primary_"+str(filename)[-8:-5]+".txt", "w")
    print("Number of analyzed frames: ", len(total_hits_counter), "Number of correct counter fractions: ", len(frac_corr_frame), file = f)
    print("Total #hits in frames/#hits in clusters = ", np.sum(total_hits_counter)/cluster_hits_counter, file = f)
    print("Correctly associated out of all hits: ", tot_corr_counter/(np.sum(total_hits_counter)/100),"%", file = f)
    print("Correctly associated out of all hits in clusters: ", tot_corr_counter/(cluster_hits_counter/100),"%", file = f)
    print("Incorrectly associated out of all hits: ", tot_uncorr_counter/(np.sum(total_hits_counter)/100),"%", file = f)
    print("Incorrectly associated out of all hits in clusters: ", tot_uncorr_counter/(cluster_hits_counter/100),"%", file = f)
    f.close()
    
    #save info as txt file
    np.savetxt("./melp/clustering/results/frac_corr_frame_primary_"+str(filename)[-8:-5]+".txt", np.array(frac_corr_frame))
    np.savetxt("./melp/clustering/results/frac_corr_clusters_frame_primary_"+str(filename)[-8:-5]+".txt", np.array(frac_corr_clusters_frame))
    np.savetxt("./melp/clustering/results/frac_uncorr_frame_primary_"+str(filename)[-8:-5]+".txt", np.array(frac_uncorr_frame))


###########################
def compare_to_tid_filename(filename, time_threshold, threshold_cluster_width, mask_type, rec_type = None, cluster_type = None):
    #get file
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e_mc = file.Get("mu3e_mchits")
    ttree_sensor = file.Get("alignment/sensors")
    ttree_tiles = file.Get("alignment/tiles")

    with HiddenPrints():
        mu3e_detector = Detector.initFromROOT(filename)

    #define counters and arrays
    frac_corr_frame            = []
    frac_corr_clusters_frame   = []
    frac_uncorr_frame          = []
    total_hits_counter         = []
    cluster_hits_counter       = 0
    tot_corr_counter           = 0
    tot_uncorr_counter         = 0

    tot_cluster_counter        = 0
    double_tid_cluster_counter = 0
    double_tid_cluster_hits_counter = 0
    corr_double_tid_cluster_counter = 0
    long_time_between_cluster_hits_counter = 0


    for frame in range(1000):#(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

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


        #---------------------------------------------------
        #count clusters with hits that are far apart in time
        #---------------------------------------------------
        long_time_between_cluster_hits_counter_tmp = 0
        for i in range(len(clusters)):
            times = clusters[i].get_times()
            min_time = min(times)
            max_time = max(times)
            if max_time - min_time > 0.5:
                long_time_between_cluster_hits_counter_tmp +=1
        long_time_between_cluster_hits_counter += long_time_between_cluster_hits_counter_tmp

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
        sel_cluster_types = ["time", "timethenspatial", "timetheniterativespatial"]
        if cluster_type in sel_cluster_types:
            #count number of clusters that have same master_tid as other cluster
            for i in range(len(clusters)):
                tot_cluster_counter += 1
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
                            if distance < threshold_cluster_width: #mm
                                corr_counter   -= len(min_times_clusters[k])
                                uncorr_counter += len(min_times_clusters[k])

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
            
    #write overall stats to txt file 
    f = open("./melp/clustering/results/efficiency_stats_tid_"+str(filename)[-8:-5]+".txt", "w")
    print("Number of analyzed frames: ", len(total_hits_counter), "Number of correct counter fractions: ", len(frac_corr_frame), file = f)
    print("Total number of hits =",np.sum(total_hits_counter), ", Identified correctly + identified incorrectly =", tot_corr_counter + tot_uncorr_counter, file = f)
    print("Identified correctly:", tot_corr_counter, file = f)
    print("Identified incorrectly:", tot_uncorr_counter, file = f)
    print("Total #hits in frames/#hits in clusters = ", np.sum(total_hits_counter)/cluster_hits_counter, file = f)
    print("Total number of clusters:",tot_cluster_counter, ", Hits:",np.sum(total_hits_counter), file = f)
    print("Number of clusters with hits that are far apart in time:", long_time_between_cluster_hits_counter, file = f)
    print("Number of clusters where tid already exists:",double_tid_cluster_counter, ", Hits:", double_tid_cluster_hits_counter, file = f)
    print("Number of clusters where tid already exists, that are accounted for:", corr_double_tid_cluster_counter, file = f)
    print("Correctly associated out of all hits: ", tot_corr_counter/(np.sum(total_hits_counter)/100),"%", file = f)
    print("Correctly associated out of all hits in clusters: ", tot_corr_counter/(cluster_hits_counter/100),"%", file = f)
    print("Incorrectly associated out of all hits: ", tot_uncorr_counter/(np.sum(total_hits_counter)/100),"%", file = f)
    print("Incorrectly associated out of all hits in clusters: ", tot_uncorr_counter/(cluster_hits_counter/100),"%", file = f)
    f.close()

    #save info as txt file
    np.savetxt("./melp/clustering/results/frac_corr_frame_tid_"+str(filename)[-8:-5]+".txt", np.array(frac_corr_frame))
    np.savetxt("./melp/clustering/results/frac_corr_clusters_frame_tid_"+str(filename)[-8:-5]+".txt", np.array(frac_corr_clusters_frame))
    np.savetxt("./melp/clustering/results/frac_uncorr_frame_tid_"+str(filename)[-8:-5]+".txt", np.array(frac_uncorr_frame))


#########################
def mt_compare_to_primary(input_files, time_threshold, mask_type, rec_type, cluster_type, i):
    print("read file ", i+1)
    print("started thread ", i+1)
    compare_to_primary_filename(input_files[0][i], time_threshold, mask_type, rec_type, cluster_type)   


#########################
def mt_compare_to_tid(input_files, time_threshold, threshold_cluster_width, mask_type, rec_type, cluster_type, i):
    print("read file ", i+1)
    print("started thread ", i+1)
    compare_to_tid_filename(input_files[0][i], time_threshold, threshold_cluster_width, mask_type, rec_type, cluster_type) 


##########
def run_mt(function_str, src, args):
    #function_str: pass function as str, src: directory of root files to analyse, args: arguments for parallelized function
    # set used threads
    unused_threads = 2 #set the number of threads you don't want to use
    print("-----------------------")
    print("Available threads = ",mp.cpu_count())
    print("Used threads = ",mp.cpu_count() - unused_threads)
    print("-----------------------")

    input_files = []

    # get list of input files
    input_files.append(glob(src))

    # map multiprocessing pool
    pool = mp.Pool(mp.cpu_count() - unused_threads)
    if function_str == "mt_compare_to_primary":
        func = partial(mt_compare_to_primary, input_files, *args)
        pool.map(func, [i for i in range(len(input_files[0]))])
    elif function_str == "mt_compare_to_tid":
        func = partial(mt_compare_to_tid, input_files, *args)
        pool.map(func, [i for i in range(len(input_files[0]))])
    else:
        raise ValueError("Function not found")

    pool.close()