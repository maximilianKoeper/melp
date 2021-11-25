import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

'''
def time_clustering_frame(filename,frame,threshold):
    clusters = {}
    cluster_counter = 0
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")

    hittimes = hittimes_in_frame (filename, frame)

    for hit_time in hittimes.values():
        for hit_time_2 in hittimes.values():
            if abs(hit_time[0] - hit_time_2[0]) < threshold:
                if cluster_counter in clusters.keys():
                    for val_arr in clusters.values():
                        if any(hit_time_2[0] in sl for sl in val_arr):
                            if [get_key_for_value(hittimes, hit_time_2[0]), hit_time_2[0]] not in clusters.values():
                                clusters[cluster_counter].append([get_key_for_value(hittimes, hit_time_2[0]), hit_time_2[0]])
                            else:
                                continue
                                
                        else: 
                            continue
                    cluster_counter += 1

                else:
                    clusters[cluster_counter] = [[get_key_for_value(hittimes, hit_time[0]), hit_time[0]], [get_key_for_value(hittimes, hit_time_2[0]), hit_time_2[0]]]

            else:
                continue

    return clusters
'''

def time_clustering_frame(filename,frame,threshold):
    clusters = {}
    cluster_counter = 0
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")

    hittimes = hittimes_in_frame (filename, frame)
    i=0 
    j=i+1
    for i in range(len(hittimes.values())):
        for j in range(len(hittimes.values())-i):
            for key in hittimes.keys():
                if abs() < threshold:
                    clusters[cluster_counter] = [[get_key_for_value(hittimes, np.array(hittimes.values())[i]), np.array(hittimes.values())[i]], [get_key_for_value(hittimes, np.array(hittimes.values())[j]), np.array(hittimes.values())[j]]]

    return clusters