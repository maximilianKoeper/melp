from melp import Detector
import ROOT
import numpy as np
import os
import sys


####################
def hittimes_in_file (ttree_mu3e):
    hittimes = {}

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hittime  = ttree_mu3e.tilehit_timestamp[hit_tile_index]
            hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

            if hit_tile in hittimes.keys():
                hittimes[hit_tile].append(hittime)
            else:
                hittimes[hit_tile] = [hittime]

    return hittimes


###################################
def hittimes_and_primaries_in_frame (ttree_mu3e):
    hittimes  = {}
    primaries = {}

    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hittime  = ttree_mu3e.tilehit_time[hit_tile_index]
        hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]
        primary  = ttree_mu3e.tilehit_primary[hit_tile_index]

        if hit_tile in hittimes.keys():
            hittimes[hit_tile].append(hittime)
            primaries[hit_tile].append(primary)
        else:
            hittimes[hit_tile]  = [hittime]
            primaries[hit_tile] = [primary]

    return hittimes, primaries


#####################
def get_key_for_value(dict, val):
    for key, value in dict.items():
        if isinstance(value, list):
            if val == value[0]:
                return key
        else:
            if val == value:
                return key

    return "key doesn't exist"


################
def get_tid_file(ttree_mu3e, ttree_mu3e_mc):
    tid = {}

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        for i in range(len(ttree_mu3e.tilehit_tile)):
            tile = ttree_mu3e.tilehit_tile[i]
            mc_i = ttree_mu3e.tilehit_mc_i[i]
            ttree_mu3e_mc.GetEntry(mc_i)
            tid[tile] = ttree_mu3e_mc.tid

    return tid


#################
def get_tid_frame(ttree_mu3e, ttree_mu3e_mc):
    tid = {}
    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile = ttree_mu3e.tilehit_tile[i]
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        tid[tile] = ttree_mu3e_mc.tid

    return tid


################################
def get_mc_primary_for_hit_frame(ttree_mu3e):
    tilehit_primary_dict = {}

    for i in range(ttree_mu3e.Ntilehit):
        tile    = ttree_mu3e.tilehit_tile[i]
        primary = ttree_mu3e.tilehit_primary[i]
        tilehit_primary_dict[tile] = primary

    return tilehit_primary_dict


################################
def get_mc_primary_for_hit_array(ttree_mu3e, cluster_tiles):
    tilehit_primary_dict = {}

    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile    = ttree_mu3e.tilehit_tile[i]
        primary = ttree_mu3e.tilehit_primary[i]
        if tile in cluster_tiles:
            tilehit_primary_dict[tile] = primary

    return tilehit_primary_dict

####################
def frame_as_cluster(ttree_mu3e):
    hit_tiles = {}

    for i in range(len(ttree_mu3e.tilehit_tile)):
        tile    = ttree_mu3e.tilehit_tile[i]
        primary = ttree_mu3e.tilehit_primary[0] #take first primary for all hits in frame
        hit_tiles[tile] = primary

    return hit_tiles


##################################
#returns all hits in frame with hid = -1,+1
def get_cluster_master_truth_frame(ttree_mu3e, ttree_mu3e_mc, frame):
    cluster_master = []
    cluster_master_primary = []

    for i in range(len(ttree_mu3e.tilehit_tile)):
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        hid  = ttree_mu3e_mc.hid
        if np.abs(hid) == 1:
            tile    = ttree_mu3e.tilehit_tile[i]
            primary = ttree_mu3e.tilehit_primary[i]
            cluster_master.append(tile)
            cluster_master_primary.append(primary)

    return cluster_master, cluster_master_primary


#####################################
#returns all hits with hid=-1,+1 in 3 frames (overlapping)
def get_cluster_primary_truth_3frames(ttree_mu3e, ttree_mu3e_mc, frame):
    cluster_primary   = []
    tilehits_3_frames = []

    #add hits in "middle" frame
    tilehits_3_frames.append(np.array(ttree_mu3e.tilehit_tile))

    #add hits in the frame before
    if frame != 1:
        ttree_mu3e.GetEntry(frame-1)
        tilehits_3_frames.append(np.array(ttree_mu3e.tilehit_tile))

    #add hits in the frame after
    if frame != ttree_mu3e.GetEntries():
        ttree_mu3e.GetEntry(frame+1)
        tilehits_3_frames.append(np.array(ttree_mu3e.tilehit_tile))

    #check frame before
    ttree_mu3e.GetEntry(frame-1)
    for i in range(len(tilehits_3_frames[1])):
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        hid  = ttree_mu3e_mc.hid
        if np.abs(hid) == 1:
            tile = ttree_mu3e.tilehit_tile[i]
            cluster_primary.append([tile, frame-1])

    #check frame after
    ttree_mu3e.GetEntry(frame+1)
    for i in range(len(tilehits_3_frames[2])):
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        hid  = ttree_mu3e_mc.hid
        if np.abs(hid) == 1:
            tile = ttree_mu3e.tilehit_tile[i]
            cluster_primary.append([tile, frame+1])

    #check middle frame
    ttree_mu3e.GetEntry(frame)
    for i in range(len(tilehits_3_frames[0])):
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        hid  = ttree_mu3e_mc.hid
        if np.abs(hid) == 1:
            tile = ttree_mu3e.tilehit_tile[i]
            cluster_primary.append([tile, frame])


    return cluster_primary


#########################################
#returns all hits with hid=-1,+1 in a single frame and the frame id 
def get_cluster_master_truth_and_frame_id(ttree_mu3e, ttree_mu3e_mc, frame):
    cluster_master         = []
    cluster_master_primary = []

    for i in range(len(ttree_mu3e.tilehit_tile)):
        mc_i = ttree_mu3e.tilehit_mc_i[i]
        ttree_mu3e_mc.GetEntry(mc_i)
        hid  = ttree_mu3e_mc.hid

        if np.abs(hid) == 1:
            tile    = ttree_mu3e.tilehit_tile[i]
            primary = ttree_mu3e.tilehit_primary[i]
            cluster_master.append([tile, frame])
            cluster_master_primary.append(primary)

    return cluster_master, cluster_master_primary


######################
def hit_tiles_in_frame(ttree_mu3e):
    hit_tiles = []

    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append(ttree_mu3e.tilehit_tile[hit_tile_index])

    return hit_tiles


######################
def get_hit_data_frame(ttree_mu3e, ttree_mu3e_mc, frames):
    hit_data = {}

    for frame in frames:
        ttree_mu3e.GetEntry(frame)
        for i in range(len(ttree_mu3e.tilehit_tile)):
            data_tmp = []
            tile     = ttree_mu3e.tilehit_tile[i]
            mc_i     = ttree_mu3e.tilehit_mc_i[i]
            ttree_mu3e_mc.GetEntry(mc_i)
            hid      = ttree_mu3e_mc.hid
            tid      = ttree_mu3e_mc.tid

            data_tmp.append(hid)
            data_tmp.append(tid)
            data_tmp.append(frame)

            hit_data[tile] = data_tmp

    hit_data_sorted_tid = sorted(hit_data.items(), key=lambda x: x[1][1], reverse=True)

    for i in range(len(hit_data_sorted_tid)):
        print("Tile: ", hit_data_sorted_tid[i][0], " hid: ",hit_data_sorted_tid[i][1][0]," tid: ", hit_data_sorted_tid[i][1][1], "frame_id: ", hit_data_sorted_tid[i][1][2])

    return hit_data_sorted_tid


#################
# decorater used to block function printing to the console
def blockPrinting(func):
    def func_wrapper(*args, **kwargs):
        # block all printing to the console
        sys.stdout = open(os.devnull, 'w')
        # call the method in question
        value = func(*args, **kwargs)
        # enable all printing to the console
        sys.stdout = sys.__stdout__
        # pass the return value of the method back
        return value

    return func_wrapper