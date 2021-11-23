from melp import Detector
import ROOT
import numpy as np

#----------------------------------------
def hittimes_in_file (ttree_mu3e):
#def hittimes_in_file (filename):
    hittimes = {}
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hittime = ttree_mu3e.tilehit_timestamp[hit_tile_index]
            hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

            if hit_tile in hittimes.keys():
                hittimes[hit_tile].append(hittime)
            else:
                hittimes[hit_tile] = [hittime]

    return hittimes

#-------------------------------------------
def hittimes_in_frame (ttree_mu3e):
#def hittimes_in_frame (filename, frame):
    hittimes = {}
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    
    #ttree_mu3e.GetEntry(frame)

    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hittime = ttree_mu3e.tilehit_timestamp[hit_tile_index]
        hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

        if hit_tile in hittimes.keys():
            hittimes[hit_tile].append(hittime)
        else:
            hittimes[hit_tile] = [hittime]

    return hittimes

#------------------------------------------
def get_key_for_value(dict, val):
    for key, value in dict.items():
        if isinstance(value, list):
            if val == value[0]:
                return key
        else:
            if val == value:
                return key

    return "key doesn't exist"

#--------------------------------------------
def get_tid_file(ttree_mu3e, ttree_mu3e_mc):
#def get_tid_file(filename):
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    #ttree_mu3e_mc = file.Get("mu3e_mchits")
    tid = {}

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        for i in range(len(ttree_mu3e.tilehit_tile)):
            tile = ttree_mu3e.tilehit_tile[i]
            mc_i = ttree_mu3e.tilehit_mc_i[i]
            ttree_mu3e_mc.GetEntry(mc_i)


            tid[tile] = ttree_mu3e_mc.tid

    return tid

#--------------------------------------------------
def get_tid_frame(ttree_mu3e, ttree_mu3e_mc):
#def get_tid_frame(filename, frame):
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    #ttree_mu3e_mc = file.Get("mu3e_mchits")
    tid = {}
    #ttree_mu3e.GetEntry(frame)
    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile = ttree_mu3e.tilehit_tile[i]
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)


        tid[tile] = ttree_mu3e_mc.tid

    return tid

#-----------------------------------------------
def get_mc_primary_for_hit_frame(ttree_mu3e):
#def get_mc_primary_for_hit_frame(filename, frame):
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    tilehit_primary_dict = {}
    #ttree_mu3e.GetEntry(frame)
    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile = ttree_mu3e.tilehit_tile[i]
        primary = ttree_mu3e.tilehit_primary[i]
        
        tilehit_primary_dict[tile] = primary

    return tilehit_primary_dict

#-------------------------------------------------
def get_mc_primary_for_hit_array(ttree_mu3e, cluster_tiles):
#def get_mc_primary_for_hit_array(filename, frame, cluster_tiles):
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    tilehit_primary_dict = {}
    #ttree_mu3e.GetEntry(frame)
    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile = ttree_mu3e.tilehit_tile[i]
        primary = ttree_mu3e.tilehit_primary[i]
        
        if tile in cluster_tiles:
            tilehit_primary_dict[tile] = primary

    return tilehit_primary_dict

#----------------------------------------------
def frame_as_cluster(ttree_mu3e):
#def frame_as_cluster(filename, frame):
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    hit_tiles = {}
    #ttree_mu3e.GetEntry(frame)
    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile = ttree_mu3e.tilehit_tile[i]
        primary = ttree_mu3e.tilehit_primary[0] #take first primary for all hits in frame
        
        hit_tiles[tile] = primary

    return hit_tiles

#---------------------------------------------
"""
def get_cluster_primary_truth_frame(filename, frame):
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e_mc = file.Get("mu3e_mchits")

    mask_primary = []
    indices_hit_time = {}

    tile_arr = []
    traj_id_arr = []
    time_hit_arr = []

    indices = []

    ttree_mu3e.GetEntry(frame)
    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile = ttree_mu3e.tilehit_tile[i]
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        traj_id = ttree_mu3e_mc.tid
        time_hit = ttree_mu3e_mc.time
        tile_arr.append(tile)
        traj_id_arr.append(traj_id)
        time_hit_arr.append(time_hit)

    for i in range(len(traj_id_arr)):
        indices_tmp = []
        for j in range(len(traj_id_arr)):
            if traj_id_arr[j != i] == traj_id_arr[i]:
                indices_tmp.append(j)

        if len(indices_tmp) > 1:
            time_hit_arr_tmp = []
            for index in indices_tmp:
                time_hit_arr_tmp.append(time_hit_arr[index])

            min_index = indices_tmp[time_hit_arr_tmp.index(min(time_hit_arr_tmp))]
            indices.append(min_index)

        else:
            indices.append(indices_tmp[0])

    for ind in indices:
        mask_primary.append(tile_arr[ind])
        

    return mask_primary
"""

def get_cluster_primary_truth_frame(ttree_mu3e, ttree_mu3e_mc):
#def get_cluster_primary_truth_frame(filename, frame):
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    #ttree_mu3e_mc = file.Get("mu3e_mchits")

    cluster_primary = []

    #ttree_mu3e.GetEntry(frame)
    for i in range(len(ttree_mu3e.tilehit_tile)):
        #tile = ttree_mu3e.tilehit_tile[i]
        #mc_i = ttree_mu3e.tilehit_mc_i[i]
        #ttree_mu3e_mc.GetEntry(mc_i)
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        hid = ttree_mu3e_mc.hid

        if np.abs(hid) == 1:
            tile = ttree_mu3e.tilehit_tile[i]
            cluster_primary.append(tile)

    return cluster_primary

#------------------------------------------
def hit_tiles_in_frame(ttree_mu3e):
#def hit_tiles_in_frame(filename, frame):
    hit_tiles = []
    #file = ROOT.TFile(filename)
    #ttree_mu3e = file.Get("mu3e")
    
    #ttree_mu3e.GetEntry(frame)

    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])

    return hit_tiles

