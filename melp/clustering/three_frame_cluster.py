import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*
import melp.clustering.tracking as clump_tr

#clustering that checks for every frame if in the frame before or after are constituents of clusters from "middle" frame



#---------------------------------------------------
#builds masks around every hit in frame with hid=1,-1 using the get neighbour function, which also works if hit is at the edge of the detector. Also returns the frame id.
def build_mask_around_cluster_master_frame_id(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    #select reconstruction/tracking method
    if rec_type == "pixelpixel":
        mask_masters = clump_tr.get_mask_masters_hitAnglePixelRec(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, matching="nearest")
    else:
        mask_masters = get_cluster_primary_truth_and_frame_id(ttree_mu3e, ttree_mu3e_mc, frame)

    mask = {}
    if mask_type == "small":
        for i in range(len(mask_masters)):
            tile_centre = [mask_masters[i][0], mask_masters[i][1]]
            tile_centre_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "up"), mask_masters[i][1]]
            tile_centre_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "down"), mask_masters[i][1]]
            tile_left_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "left"), mask_masters[i][1]]
            tile_right_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "right"), mask_masters[i][1]]

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, 
                                tile_left_centre, tile_right_centre])

            if False in mask_tmp:
                mask_tmp = np.array([list(x) for x in mask_tmp if x[0] != False])

            mask[tile_centre[0]] = mask_tmp

    if mask_type == "medium":
        for i in range(len(mask_masters)):
            tile_centre = [mask_masters[i][0], mask_masters[i][1]]
            tile_centre_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "up"), mask_masters[i][1]]
            tile_centre_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "down"), mask_masters[i][1]]
            tile_left_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "left"), mask_masters[i][1]]
            tile_left_top = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "up"), mask_masters[i][1]]
            tile_left_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "down"), mask_masters[i][1]]
            tile_right_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "right"), mask_masters[i][1]]
            tile_right_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "up"), mask_masters[i][1]]
            tile_right_bootom = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "down"), mask_masters[i][1]]

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, 
                                tile_left_centre, tile_left_bottom, tile_right_top, tile_right_centre, 
                                tile_right_bootom])

            if False in mask_tmp:
                mask_tmp = np.array([list(x) for x in mask_tmp if x[0] != False])

            mask[tile_centre[0]] = mask_tmp

    if mask_type == "big":
        for i in range(len(mask_masters)):
            tile_centre = [mask_masters[i][0], mask_masters[i][1]]
            tile_centre_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "up"), mask_masters[i][1]]
            tile_centre_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "down"), mask_masters[i][1]]
            tile_left_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "left"), mask_masters[i][1]]
            tile_left_top =  [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "up"), mask_masters[i][1]]
            tile_left_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "down"), mask_masters[i][1]]
            tile_right_centre = [mu3e_detector.TileDetector.getNeighbour(tile_centre[0], "right"), mask_masters[i][1]]
            tile_right_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "up"), mask_masters[i][1]]
            tile_right_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "down"), mask_masters[i][1]]

            tile_centre_far_top = [mu3e_detector.TileDetector.getNeighbour(tile_centre_top[0], "up"), mask_masters[i][1]]
            tile_centre_far_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_centre_bottom[0], "down"), mask_masters[i][1]]
            tile_left_far_centre = [mu3e_detector.TileDetector.getNeighbour(tile_left_centre[0], "left"), mask_masters[i][1]]
            tile_left_far_top =  [mu3e_detector.TileDetector.getNeighbour(tile_left_top[0], "up"), mask_masters[i][1]]
            tile_left_far_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_bottom[0], "down"), mask_masters[i][1]]
            tile_right_far_centre = [mu3e_detector.TileDetector.getNeighbour(tile_right_centre[0], "right"), mask_masters[i][1]]
            tile_right_far_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_top[0], "up"), mask_masters[i][1]]
            tile_right_far_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_right_bottom[0], "down"), mask_masters[i][1]]

            tile_far_left_top = [mu3e_detector.TileDetector.getNeighbour(tile_left_top[0], "left"), mask_masters[i][1]]
            tile_far_left_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_left_bottom[0], "left"), mask_masters[i][1]]
            tile_far_right_top = [mu3e_detector.TileDetector.getNeighbour(tile_right_top[0], "right"), mask_masters[i][1]]
            tile_far_right_bottom = [mu3e_detector.TileDetector.getNeighbour(tile_right_bottom[0], "right"), mask_masters[i][1]]

            mask_tmp = np.array([tile_centre, tile_centre_top, tile_centre_bottom, tile_left_top, tile_left_centre, tile_left_bottom,
                                   tile_right_top, tile_right_centre, tile_right_bottom, tile_centre_far_top, tile_centre_far_bottom, 
                                   tile_left_far_centre, tile_left_far_top, tile_left_far_bottom, tile_right_far_centre, 
                                   tile_right_far_top, tile_right_far_bottom, tile_far_left_top, tile_far_left_bottom, tile_far_right_top, 
                                   tile_far_right_bottom])
            
            if False in mask_tmp:
                mask_tmp = np.array([list(x) for x in mask_tmp if x[0] != False])

            mask[tile_centre[0]] = mask_tmp

    return mask 

#----------------------------------------------------
#builds clusters in the masks around hit with hid=1,-1 according to primaries. Also checks for cluster constituents in neighbouring frames
def build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector: melp.Detector, frame, mask_type, rec_type = None):
    clusters = {} #keys: "master"-tile; values: rest of cluster

    #get masks around master tile
    primary_masks = build_mask_around_cluster_master_frame_id(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector, frame, mask_type, rec_type)

    keys = []
    values = []
    for key in primary_masks.keys():
        keys.append(key)
        values.append(primary_masks[key])

    #get all tiles that have been hit in frame
    hit_tiles = []
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append([ttree_mu3e.tilehit_tile[hit_tile_index], frame])

    ttree_mu3e.GetEntry(frame+1)
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append([ttree_mu3e.tilehit_tile[hit_tile_index], frame+1])

    ttree_mu3e.GetEntry(frame-1)
    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hit_tiles.append([ttree_mu3e.tilehit_tile[hit_tile_index], frame-1])
    
    ttree_mu3e.GetEntry(frame)

    #build clusters around mask master tiles
    for tile in hit_tiles:
        if tile[0] < 300000: #just left recurl station
            cluster_tmp = []
            cluster_primary_tmp = 0
            if tile[0] not in keys: #if not primary
                for i in range(len(values)):
                    for j in range(len(values[i])):
                        if tile[0] == values[i][j][0]: 
                            #cluster_tmp.append(tile)
                            if len(cluster_tmp) != 0:
                                for k in cluster_tmp: #TODO: check if hit is already there from another frame (so the same tile isn't added twice when hit in multiple frames)
                                    if tile[0] == k[0]:
                                        continue
                                    else:
                                        cluster_tmp.append(tile)
                            elif len(cluster_tmp) == 0:
                                cluster_tmp.append(tile)
                            if list(values[i][0]) not in cluster_tmp:
                                cluster_tmp.insert(0,list(values[i][0]))
                            cluster_primary_tmp = keys[i]

            if cluster_primary_tmp != 0:
                if cluster_primary_tmp not in clusters.keys() and len(cluster_tmp) > 0:          
                    clusters[cluster_primary_tmp] = cluster_tmp
                elif len(cluster_tmp) > 0: 
                    clusters[cluster_primary_tmp].append(cluster_tmp[1])

    for i in keys:
        if i not in clusters.keys():
            index_tmp = keys.index(i)
            clusters[i] = [list(values[index_tmp][0])]

    return clusters

#------------------------------------------------------------------
#checks if a cluster consists of hits in multiple frames
def check_for_multiple_frame_clusters(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_cluster_hits_counter = 0
    mult_frame_cluster_hits_counter = 0

    frac_mult_frame_cluster_hits = []

    #counting
    #for frame in range(frames_to_analyze):
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)

        mult_frame_cluster_hits_counter_tmp = 0
        total_cluster_hits_counter_tmp = 0

        #Printing status info
        if int(frame) % 2000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        clusters_frame = build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)

        for key in clusters_frame.keys():
            value = np.array(clusters_frame[key])
            total_cluster_hits_counter_tmp += len(value)
            for i in range(len(value)): #checks if hit in cluster has different frame id than the "master" hit
                if value[i][1] != value[0][1]:
                    mult_frame_cluster_hits_counter_tmp += 1

        total_cluster_hits_counter += total_cluster_hits_counter_tmp
        mult_frame_cluster_hits_counter += mult_frame_cluster_hits_counter_tmp

        if total_cluster_hits_counter_tmp != 0:
            frac_mult_frame_cluster_hits.append(mult_frame_cluster_hits_counter_tmp/total_cluster_hits_counter_tmp)

    print("Progress: 100 %","of ", frames_to_analyze, " frames")
    print("Hits in cluster in different frame than master out of all hits in clusters: ", mult_frame_cluster_hits_counter/(total_cluster_hits_counter/100), "%")

    return frac_mult_frame_cluster_hits, total_cluster_hits_counter, mult_frame_cluster_hits_counter

#-------------------------------------------
#checks for clusters with multiple hits in one tile if the hits come from different frames
def check_for_mult_hit_tiles_diff_frame(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    
    #set counters
    total_cluster_hits_counter = 0
    double_hit_counter = 0

    #counting
    #for frame in range(frames_to_analyze):
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)

        double_hit_counter_tmp = 0
        total_cluster_hits_counter_tmp = 0

        #Printing status info
        if int(frame) % 2000 == 0:
            print("Progress: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        clusters_frame = build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)

        for key in clusters_frame.keys():
            cluster = clusters_frame[key]
            total_cluster_hits_counter_tmp += len(cluster)
            for hit1 in cluster:
                for hit2 in cluster:
                    if hit1[0] == hit2[0] and hit1[1] != hit2[1]:
                        double_hit_counter_tmp += 1

        double_hit_counter += double_hit_counter_tmp
        total_cluster_hits_counter += total_cluster_hits_counter_tmp

    if double_hit_counter != 0:
        print("Tiles in cluster hit in multiple frames out of all hits in clusters: ", double_hit_counter/(total_cluster_hits_counter/100), "%")

    return double_hit_counter, total_cluster_hits_counter

#----------------------------------------------
#checks clusters retuned by build_clusters_in_masks_with_neighbours and deletes hits in overlapping frames if hit belongs to cluster in 
#"main" frame. Then returns clusters without double hits to reduce the number of hits that can't get associated
def del_double_hits_in_3_frame_cluster(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, mask_type, number_of_frames = None, rec_type = None):
    #set frame number
    if number_of_frames == None:
        frames_to_analyze = ttree_mu3e.GetEntries()
    else:
        frames_to_analyze = number_of_frames
    

    hits_all_frames = {} #contains all hits from all frames. Key: frame_id, Value: List of all hits in frame with frame_id
    hits_all_frames_counter_before = 0
    removed_hits_counter = 0

    #build dictionary with hits from all frames
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 5000 == 0:
            print("Progress of building dictionary with all hits: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        hit_tiles_with_frame_id = []
        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hit_tiles_with_frame_id.append([ttree_mu3e.tilehit_tile[hit_tile_index], frame])
        
        hits_all_frames[frame] = hit_tiles_with_frame_id
        hits_all_frames_counter_before += len(hit_tiles_with_frame_id)

    print("Progress of building dictionary with all hits: 100 %","of ", frames_to_analyze, " frames")

    #loop over frames again
    for frame in np.arange(2, frames_to_analyze-2, 1):
        ttree_mu3e.GetEntry(frame)
        #Printing status info
        if frame % 1000 == 0:
            print("Progress of removing double hits: ", np.round(frame / frames_to_analyze * 100), " %","of ", frames_to_analyze, " frames", end='\r')

        #get clusters for single frame 
        clusters_frame = build_clusters_in_masks_with_neighbours(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles,  mu3e_detector, frame, mask_type, rec_type)

        #check if there is hit from other frame and add them to hits_from_other_frame_in_cluster
        hits_from_other_frame_in_cluster = []
        for key in clusters_frame.keys():
            value = np.array(clusters_frame[key])
            for i in range(len(value)): #checks if hit in cluster has different frame id than the "master" hit
                if value[i][1] != value[0][1]:
                    hits_from_other_frame_in_cluster.append([value[i][0], value[i][1]]) #add hits in different frame to list

        #get hits from other frames if hits_from_other_frame_in_cluster is not empty
        if len(hits_from_other_frame_in_cluster) != 0:
            #check which frame to get (faster if its just one)
            for hit in hits_from_other_frame_in_cluster:
                if hit[1] == frame+1:
                    #get hits in frame +1
                    hits_frame_plus = hits_all_frames[frame+1]
                    if hit in hits_frame_plus:
                        hits_frame_plus.remove(hit) #remove double hit
                        removed_hits_counter += 1
                    hits_all_frames[frame+1] = hits_frame_plus #replace the value in hits_all_frames with new value without double hits
                elif hit[1] == frame-1:
                    #get hits in frame -1
                    hits_frame_minus = hits_all_frames[frame-1]
                    if hit in hits_frame_minus:
                        hits_frame_minus.remove(hit) #remove double hit
                        removed_hits_counter += 1
                    hits_all_frames[frame-1] = hits_frame_minus #replace the value in hits_all_frames with new value without double hits
                else:
                    print("ERROR: Difference of frame_ids > 1")

        else:
            continue

    print("Progress of removing double hits: 100 %","of ", frames_to_analyze, " frames")
    
    #print number of hits before removing doubles
    print("Number of hits before removing doubles: ", hits_all_frames_counter_before)
    #print number of hits after removing doubles
    hits_all_frames_counter_after = hits_all_frames_counter_before - removed_hits_counter
    print("Number of hits after removing doubles: ", hits_all_frames_counter_after)
    #print number of removed hits
    print("Number of removed hits: ", removed_hits_counter)

    return hits_all_frames, hits_all_frames_counter_after