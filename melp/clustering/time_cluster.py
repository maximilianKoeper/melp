import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*


def time_clustering_frame(filename,frame,threshold):
    clusters = {}
    cluster_counter = 0
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")

    hittimes = hittimes_in_frame (filename, frame)
    for hit_time in hittimes.values():
        for hit_time_2 in hittimes.values():
            if (np.array(hit_time) - np.array(hit_time_2)) < threshold:
                cluster_counter += 1
                clusters[cluster_counter] = [hit_time, hit_time_2] 

            else:
                continue

    return clusters

