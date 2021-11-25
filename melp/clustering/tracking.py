import ROOT
import numpy as np
import melp
from melp import Detector

from melp.clustering.misc import*

import melp.clustering.spatial_cluster as sclump

# ------------------------------------
def __Get_HID_from_MC_I (ttree_mu3e_mc, mc_i):
    ttree_mu3e_mc.GetEntry(mc_i)
    hid = ttree_mu3e_mc.hid

    return hid

# ------------------------------------
def __Get_TID_from_MC_I (ttree_mu3e_mc, mc_i):
    ttree_mu3e_mc.GetEntry(mc_i)
    tid = ttree_mu3e_mc.tid

    return tid

# ------------------------------------
def __Get_Sensor_IDs_from_Frame_ID (ttree_mu3e, frame):
    ttree_mu3e.GetEntry(frame)

    # Get Sensor IDs
    hit_id_arr = []

    sensor_id_arr = []
    for j in ttree_mu3e.hit_pixelid:
        sensor_id_arr.append(j)

    # Get MC index
    mc_i_arr = []
    for j in ttree_mu3e.hit_mc_i:
        mc_i_arr.append(j)

    return sensor_id_arr, mc_i_arr


#--------------------------------------------
def hitAnglePixelRec(ttree_mu3e, ttree_mu3e_mc, ttree_sensor, ttree_tiles, mu3e_detector: melp.Detector, number_of_frames = None, matching="nearest"):
        hit_type = "primary"
        ana_tpye = "SensorMatching" + matching


        # initialize dictionaries
        tilehit_tile_dic = {}
        tile_mc_i = {}
        tile_id_pos      = {}
        tile_id_dir      = {}
        sensor_id_index  = {}

        for i in range(ttree_mu3e.GetEntries()):
            ttree_mu3e.GetEntry(i)
            temp_arr_1 = []
            for j in ttree_mu3e.tilehit_mc_i:
                temp_arr_1.append(j)
            tile_mc_i[i] = temp_arr_1

            temp_arr_2 = []
            for j in ttree_mu3e.tilehit_tile:
                temp_arr_2.append(j)
            tilehit_tile_dic[i] = temp_arr_2

        for i in range(ttree_tiles.GetEntries()):
            ttree_tiles.GetEntry(i)
            # direction
            xyz = []
            xyz.append(ttree_tiles.dirx)
            xyz.append(ttree_tiles.diry)
            xyz.append(ttree_tiles.dirz)

            tile_id_dir[ttree_tiles.sensor] = xyz

            # position
            tile_xyz = []
            tile_xyz.append(ttree_tiles.posx)
            tile_xyz.append(ttree_tiles.posy)
            tile_xyz.append(ttree_tiles.posz)
            tile_id_pos[ttree_tiles.sensor] = tile_xyz

        for i in range(ttree_sensor.GetEntries()):
            ttree_sensor.GetEntry(i)
            sensor_id_index[ttree_sensor.sensor] = i
  
        # counters
        hid_discard = 0
        hid_ok      = 0
        tid_discard = 0
        tid_ok      = 0

        #set frame number
        if number_of_frames == None:
            frames_to_analyze = ttree_mu3e.GetEntries()
        else:
            frames_to_analyze = number_of_frames

        # Define Arrays for result
        z_arr  = []
        id_arr = []

        # loop over all Root frames
        for i in range(frames_to_analyze):
            # loop over all tile hits in one Root frame
            for u in range(len(tilehit_tile_dic[i])):
                ##################################
                # HID CHECK
                # only primary hit gets analyzed
                ##################################
                tile_id  = tilehit_tile_dic[i][u]
                hid_test = __Get_HID_from_MC_I(ttree_mu3e_mc, tile_mc_i[i][u])
                if hid_test != 1:
                    hid_discard += 1
                    continue
                hid_ok += 1

                sensor_ids_tmp, sensor_frame_mc_i_tmp = __Get_Sensor_IDs_from_Frame_ID(ttree_mu3e, i)

                tmp_distance_tile_to_pixel_layer2    = []
                tmp_distance_pixel_to_pixel          = []

                tile_pos = tile_id_pos[tile_id]

                # TID for tile
                tid_tile_test   = __Get_TID_from_MC_I(ttree_mu3e_mc, tile_mc_i[i][u])

                ##################################
                # TID CHECK
                # check for matching sensor and tile hits
                ##################################
                sensor_id_tid   = []
                sensor_mc_i_tid = []
                for g in range(len(sensor_ids_tmp)):
                    tid_sensor_test = __Get_TID_from_MC_I(ttree_mu3e_mc, sensor_frame_mc_i_tmp[g])
                    if tid_sensor_test != tid_tile_test:
                        tid_discard += 1
                        continue
                    tid_ok += 1
                    sensor_id_tid.append(sensor_ids_tmp[g])
                    sensor_mc_i_tid.append(sensor_frame_mc_i_tmp[g])

                #split pixel ids into different pixel layers
                pixel_ids = []
                for l in sensor_id_tid:
                    pixel_ids.append(l >> 16)

                sensor_ids_layer2        = [] #pixel ids 2000 <= ID < 3000 || 14000 <= ID < 15200
                sensor_frame_mc_i_layer2 = []
                sensor_ids_layer3        = [] #pixel ids 3000 <= ID < 4000 || 15200 <= ID < 16500
                sensor_frame_mc_i_layer3 = []

                for k in pixel_ids:
                    #if (k >= 2000 and k < 3000) or (k >= 14000 and k < 15200):
                    if (k >= 10000 and k < 11500) or (k >= 14000 and k < 15200):
                        index_id_2 = np.where(np.array(pixel_ids) == k)
                        sensor_ids_layer2.append(sensor_id_tid[index_id_2[0][0]])
                        sensor_frame_mc_i_layer2.append(sensor_mc_i_tid[index_id_2[0][0]])

                    #elif (k >= 3000 and k < 4000) or (k >= 15200 and k < 16500):
                    elif (k >= 11500 and k < 12500) or (k >= 15200 and k < 16500):
                        index_id_3 = np.where(np.array(pixel_ids) == k)
                        sensor_ids_layer3.append(sensor_id_tid[index_id_3[0][0]])
                        sensor_frame_mc_i_layer3.append(sensor_mc_i_tid[index_id_3[0][0]])


                # loop over all pixel hits in one Root frame
                pixel_pos_layer2 = []
                pixel_pos_layer3 = []
                # find distance tile to pixel (in layer 2)
                for v in range(len(sensor_ids_layer2)):
                    pixel_id_layer2  = sensor_ids_layer2[v]
                    pixel_pos_layer2 = __Get_Sensor_Pos_from_Pixel_ID(pixel_id_layer2)
                    distance_layer2  = np.sqrt((tile_pos[0]-pixel_pos_layer2[0])**2 + (tile_pos[1]-pixel_pos_layer2[1])**2 + (tile_pos[2]-pixel_pos_layer2[2])**2)

                    tmp_distance_tile_to_pixel_layer2.append(distance_layer2)
                    ##################################
                    # the nearest matching sensor hits are used to approximate the trajectory
                    ##################################
                    # tmp_distance_tile_to_pixel can be zero!
                if len(tmp_distance_tile_to_pixel_layer2) != 0:
                    index_2     = np.where(tmp_distance_tile_to_pixel_layer2 == np.min(tmp_distance_tile_to_pixel_layer2))[0][0]
                    sensor_id_layer2 = sensor_ids_layer2[index_2]
                    pixel_pos_layer2 = __Get_Sensor_Pos_from_Pixel_ID(sensor_id_layer2)

                else:
                    continue

                # find distance pixel to pixel
                for w in range(len(sensor_ids_layer3)):
                    pixel_id_layer3  = sensor_ids_layer3[w]
                    pixel_pos_layer3 = __Get_Sensor_Pos_from_Pixel_ID(pixel_id_layer3)
                    distance_pixel  = np.sqrt((pixel_pos_layer2[0]-pixel_pos_layer3[0])**2 + (pixel_pos_layer2[1]-pixel_pos_layer3[1])**2 + (pixel_pos_layer2[2]-pixel_pos_layer3[2])**2)

                    tmp_distance_pixel_to_pixel.append(distance_pixel)
                    ##################################
                    # the nearest matching sensor hits are used to approximate the trajectory
                    ##################################
                    # tmp_distance_tile_to_pixel can be zero!
                if len(tmp_distance_pixel_to_pixel) != 0:
                    index_3     = np.where(tmp_distance_pixel_to_pixel == np.min(tmp_distance_pixel_to_pixel))[0][0]
                    sensor_id_layer3 = sensor_ids_layer3[index_3]
                    pixel_pos_layer3 = __Get_Sensor_Pos_from_Pixel_ID(sensor_id_layer3)

                if np.array(pixel_pos_layer3).size == 0:
                    continue
                
                vector_sensor_layers = np.array(pixel_pos_layer3) - np.array(pixel_pos_layer2)
                
                tile_norm = np.array(tile_id_dir[tile_id])

                z_arr.append(tile_pos[2])
                id_arr.append(tile_id)

            # Print progress
            if i % 1000 == 0 and i != 0:
                print(round((i/n)*100,2), "%")
        print("100%")

        print("HID CHECK: ", hid_ok, " of " , hid_ok+ hid_discard, "ok")
        print("TID CHECK: ", tid_ok, " of " , tid_ok+ hid_discard, "ok")
        print("Total Events with matching Tile and Sensor Hit: ", len(z_arr), " of: ", hid_ok, " primary Tile hits")

        result_z     = np.array(z_arr)
        result_id    = np.array(id_arr)

        return result_z, result_id

