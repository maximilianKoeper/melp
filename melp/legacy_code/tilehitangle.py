import ROOT
import numpy as np
from melp.libs import mathfunctions as mf
from melp.libs import helices as hl

import gc


class TileHitAngle:
    def __init__(self, filename, output):
        self.filename = filename
        self.output = output

        self.result_z = np.zeros(0)
        self.result_angle = np.zeros(0)
        self.result_id = np.zeros(0)

        self.file = ROOT.TFile(filename)
        self.mu3e_mchits = self.file.Get("mu3e_mchits")
        self.mu3e = self.file.Get("mu3e")
        self.sensor = self.file.Get("alignment/sensors")
        self.tiles = self.file.Get("alignment/tiles")

        self.tilehit_tile_dic = {}
        self.tile_mc_i = {}
        self.tile_id_pos = {}
        self.tile_id_dir = {}

        self.sensor_id_index = {}

        self.hit_type = ""
        self.ana_tpye = ""
        self.angle = ""

        # initialize dictionaries
        for i in range(self.mu3e.GetEntries()):
            self.mu3e.GetEntry(i)
            temp_arr_1 = []
            for j in self.mu3e.tilehit_mc_i:
                temp_arr_1.append(j)
            self.tile_mc_i[i] = temp_arr_1

            temp_arr_2 = []
            for j in self.mu3e.tilehit_tile:
                temp_arr_2.append(j)
            self.tilehit_tile_dic[i] = temp_arr_2

        for i in range(self.tiles.GetEntries()):
            self.tiles.GetEntry(i)
            # direction
            xyz = [self.tiles.dirx, self.tiles.diry, self.tiles.dirz]

            self.tile_id_dir[self.tiles.sensor] = xyz

            # position
            tile_xyz = [self.tiles.posx, self.tiles.posy, self.tiles.posz]
            self.tile_id_pos[self.tiles.sensor] = tile_xyz

        for i in range(self.sensor.GetEntries()):
            self.sensor.GetEntry(i)
            self.sensor_id_index[self.sensor.sensor] = i

    #####################
    # private functions #
    #####################

    def __Get_TID_from_Frame_ID(self, mc_i):
        self.mu3e_mchits.GetEntry(mc_i)
        tid = self.mu3e_mchits.tid

        return tid

    # ------------------------------------
    def __Get_HID_from_MC_I(self, mc_i):
        self.mu3e_mchits.GetEntry(mc_i)
        hid = self.mu3e_mchits.hid

        return hid

    # ------------------------------------
    def __Get_TID_from_MC_I(self, mc_i):
        self.mu3e_mchits.GetEntry(mc_i)
        tid = self.mu3e_mchits.tid

        return tid

    # ------------------------------------
    def __Get_Sensor_IDs_from_Frame_ID(self, entry):
        self.mu3e.GetEntry(entry)

        # Get Sensor IDs
        hit_id_arr = []

        sensor_id_arr = []
        for j in self.mu3e.hit_pixelid:
            sensor_id_arr.append(j)

        # Get MC index
        mc_i_arr = []
        for j in self.mu3e.hit_mc_i:
            mc_i_arr.append(j)

        return sensor_id_arr, mc_i_arr

    # ------------------------------------
    def __Get_Sensor_Pos_from_Pixel_ID(self, pixelid):
        pixel = pixelid >> 16
        pixel_index = self.sensor_id_index[pixel]
        self.sensor.GetEntry(pixel_index)

        row_param = pixelid & 0xFF
        col_param = (pixelid >> 8) & 0xFF

        sensor_pos_vxyz = [self.sensor.vx, self.sensor.vy, self.sensor.vz]

        sensor_pos_col = [self.sensor.colx, self.sensor.coly, self.sensor.colz]

        sensor_pos_row = [self.sensor.rowx, self.sensor.rowy, self.sensor.rowz]

        pos = np.array(sensor_pos_vxyz) + (col_param + 0.5) * np.array(sensor_pos_col) + (row_param + 0.5) * np.array(
            sensor_pos_row)
        return pos

    # ------------------------------------
    def __Get_Traj_from_TID(self, entry, tid):
        self.mu3e.GetEntry(entry)
        traj_id = np.array(self.mu3e.traj_ID)
        index = np.where(traj_id == tid)[0]
        try:
            index = int(index[0])
        except:
            return None, None, None, None, None, None, None

        vx = self.mu3e.traj_vx
        vy = self.mu3e.traj_vy
        vz = self.mu3e.traj_vz
        px = self.mu3e.traj_px
        py = self.mu3e.traj_py
        pz = self.mu3e.traj_pz
        type = self.mu3e.traj_type

        return vx[index], vy[index], vz[index], px[index], py[index], pz[index], type[index]

    # ------------------------------------
    def __Get_Momentum_from_MC_I(self, mc_i):
        self.mu3e_mchits.GetEntry(mc_i)

        px = self.mu3e_mchits.p_in_x
        py = self.mu3e_mchits.p_in_y
        pz = self.mu3e_mchits.p_in_z

        p_xyz = np.zeros(3)
        p_xyz[0] = px
        p_xyz[1] = py
        p_xyz[2] = pz

        return p_xyz

    # ------------------------------------
    def __Get_PDG_from_MC_I(self, mc_i):
        self.mu3e_mchits.GetEntry(mc_i)

        return self.mu3e_mchits.pdg

    #####################
    # public  functions #
    #####################

    def hitAngleTruth(self, n=0, angle="norm", hit_type="primary", particle_type="all"):

        self.hit_type = hit_type
        self.ana_tpye = "Truth" + particle_type
        self.angle = angle

        if n > len(self.tilehit_tile_dic) or n == 0:
            n = len(self.tilehit_tile_dic)
        print("Frames to analyze: ", n, " of ", len(self.tilehit_tile_dic))

        # Define Arrays for result
        angle_arr = []
        z_arr = []
        id_arr = []

        # loop over all Root frames
        for i in range(n):
            # for i in range(100):
            # loop over all tile hits in one Root frame
            for u in range(len(self.tilehit_tile_dic[i])):

                tile_id = self.tilehit_tile_dic[i][u]
                #############
                # HID CHECK #
                #############
                hid_test = self.__Get_HID_from_MC_I(self.tile_mc_i[i][u])
                if hit_type == "primary":
                    # only primary hit gets analyzed
                    if hid_test != 1:
                        continue
                elif hit_type == "secondary":
                    if hid_test != 1:
                        continue
                elif hit_type == "all":
                    pass
                else:
                    raise ValueError("hit_type: not supported")

                #############
                # PDG Check #
                #############
                pdg = self.__Get_PDG_from_MC_I(self.tile_mc_i[i][u])
                if particle_type == "electron":
                    if pdg != 11:
                        continue
                elif particle_type == "positron":
                    if pdg != -11:
                        continue
                elif particle_type == "all":
                    pass
                else:
                    raise ValueError("particle_type: not supported")

                tile_pos = self.tile_id_pos[tile_id]

                p_xyz = self.__Get_Momentum_from_MC_I(self.tile_mc_i[i][u])

                if angle == "norm":
                    angle_arr.append(mf.angle_between(p_xyz, self.tile_id_dir[tile_id]))
                elif angle == "theta":
                    angle_arr.append(mf.angle_between(p_xyz, np.array([0, 0, -1])))
                elif angle == "phi":
                    vector = -np.array(self.tile_id_dir[tile_id])
                    angle_arr.append(-mf.angle_between_phi(p_xyz[0:2], vector[0:2]))
                else:
                    raise ValueError('ERROR: angle != [norm, theta, phi]')

                z_arr.append(tile_pos[2])
                id_arr.append(tile_id)

            if i % 1000 == 0 and i != 0:
                print(round((i / n) * 100, 2), "%  |  Hits:", len(z_arr))
        print("100%")

        self.result_z = np.array(z_arr)
        self.result_angle = np.array(angle_arr)
        self.result_id = np.array(id_arr)

        return self.result_z, self.result_angle, self.result_id

    # ------------------------------------

    def hitAngleRec(self, n=0, angle="norm", matching="nearest"):
        """
            TODO:
                - add new options for sensor tile matching (sensor cluster)
        """
        self.hit_type = "primary"
        self.ana_tpye = "SensorMatching" + matching
        self.angle = angle

        # counters
        hid_discard = 0
        hid_ok = 0
        tid_discard = 0
        tid_ok = 0

        # Check Argument
        if n > len(self.tilehit_tile_dic) or n == 0:
            n = len(self.tilehit_tile_dic)
        print("Frames to analyze: ", n, " of ", len(self.tilehit_tile_dic))

        # Define Arrays for result
        angle_sensor_tile = []
        z_arr = []
        id_arr = []

        # loop over all Root frames
        for i in range(n):
            # for i in range(100):
            # loop over all tile hits in one Root frame
            for u in range(len(self.tilehit_tile_dic[i])):

                ##################################
                # HID CHECK
                # only primary hit gets analyzed
                ##################################
                tile_id = self.tilehit_tile_dic[i][u]
                hid_test = self.__Get_HID_from_MC_I(self.tile_mc_i[i][u])
                if hid_test != 1:
                    hid_discard += 1
                    continue
                hid_ok += 1

                tile_pos = self.tile_id_pos[tile_id]

                tmp_distance_tile_to_pixel = []
                sensor_ids, sensor_frame_mc_i = self.__Get_Sensor_IDs_from_Frame_ID(i)

                # TID for tile
                tid_tile_test = self.__Get_TID_from_MC_I(self.tile_mc_i[i][u])

                # loop over all pixel hits in one Root frame
                for v in range(len(sensor_ids)):
                    ##################################
                    # TID CHECK
                    # check for matching sensor and tile hits
                    ##################################
                    tid_sensor_test = self.__Get_TID_from_MC_I(sensor_frame_mc_i[v])
                    if tid_sensor_test != tid_tile_test:
                        tid_discard = tid_discard + 1
                        tmp_distance_tile_to_pixel.append(10000)
                        continue
                    tid_ok = tid_ok + 1

                    pixel_id = sensor_ids[v]
                    pixel_pos = self.__Get_Sensor_Pos_from_Pixel_ID(pixel_id)
                    distance = np.sqrt((tile_pos[0] - pixel_pos[0]) ** 2 + (tile_pos[1] - pixel_pos[1]) ** 2 + (
                            tile_pos[2] - pixel_pos[2]) ** 2)

                    tmp_distance_tile_to_pixel.append(distance)

                # print(tmp_distance_tile_to_pixel)
                ##################################
                # the nearest matching sensor hit is used to approximate the trajectory
                ##################################
                # tmp_distance_tile_to_pixel can be zero!
                if len(tmp_distance_tile_to_pixel) != 0 and min(tmp_distance_tile_to_pixel) < 150:
                    index = np.where(tmp_distance_tile_to_pixel == min(tmp_distance_tile_to_pixel))[0][0]
                    sensor_id = sensor_ids[index]

                    pixel_pos = self.__Get_Sensor_Pos_from_Pixel_ID(sensor_id)

                    vector_sensor_tile = np.array(pixel_pos) - np.array(tile_pos)

                    if angle == "norm":
                        angle_sensor_tile.append(mf.angle_between(vector_sensor_tile, self.tile_id_dir[tile_id]))
                    elif angle == "theta":
                        angle_sensor_tile.append(mf.angle_between(vector_sensor_tile, np.array([0, 0, 1])))
                    elif angle == "phi":
                        vector = np.array(self.tile_id_dir[tile_id])
                        angle_sensor_tile.append(-mf.angle_between_phi(vector_sensor_tile[0:2], vector[0:2]))
                    else:
                        raise ValueError('ERROR: angle != [norm, theta, phi]')
                    z_arr.append(tile_pos[2])
                    id_arr.append(tile_id)

            # Print progress
            if i % 1000 == 0 and i != 0:
                print(round((i / n) * 100, 2), "%")
        print("100%")

        print("HID CHECK: ", hid_ok, " of ", hid_ok + hid_discard, "ok")
        print("TID CHECK: ", tid_ok, " of ", tid_ok + hid_discard, "ok")
        print("Total Events with matching Tile and Sensor Hit: ", len(z_arr), " of: ", hid_ok, " primary Tile hits")

        self.result_z = np.array(z_arr)
        self.result_angle = np.array(angle_sensor_tile)
        self.result_id = np.array(id_arr)

        return self.result_z, self.result_angle, self.result_id

    # ------------------------------------

    def hitAngleHelix(self, n=0, angle="norm"):
        """
            TODO:
                [done] get trajectory ID for tilehit
                [done] get trajectory information from ID
                [done] get angle from helices
                [done] add angle "theta" and "norm"
                [...] testing
                    [...] phi
                [use mt] improve speed
        """
        self.hit_type = "primary"
        self.ana_tpye = "Helix"
        self.angle = angle

        # counters
        hid_discard = 0
        hid_ok = 0
        no_traj = 0

        if n > len(self.tilehit_tile_dic) or n == 0:
            n = len(self.tilehit_tile_dic)
        print("Frames to analyze: ", n, " of ", len(self.tilehit_tile_dic))

        # Define Arrays for result
        angle_arr = []
        z_arr = []
        id_arr = []

        # loop over all Root frames
        for i in range(n):
            # for i in range(100):
            # loop over all tile hits in one Root frame
            for u in range(len(self.tilehit_tile_dic[i])):

                ##################################
                # HID CHECK
                # only primary hit gets analyzed
                ##################################
                tile_id = self.tilehit_tile_dic[i][u]
                hid_test = self.__Get_HID_from_MC_I(self.tile_mc_i[i][u])
                if hid_test != 1:
                    hid_discard += 1
                    continue
                hid_ok += 1

                tile_pos = self.tile_id_pos[tile_id]

                tid_tile_test = self.__Get_TID_from_MC_I(self.tile_mc_i[i][u])

                vx, vy, vz, px, py, pz, traj_type = self.__Get_Traj_from_TID(i, tid_tile_test)
                if traj_type is None:
                    # No matching trajectory was found
                    no_traj += 1
                    continue

                # electron or position
                # if type == 2 or type == 3:
                # TODO: dont mix electrons with positions

                type_1 = abs(int(repr(traj_type)[-1]))
                if type_1 == 1 or type_1 == 2:
                    helix = hl.Helices(vx, vy, vz, px, py, pz, type_1, tile_pos)
                    angle_arr.append(helix.hitAngle(self.tile_id_dir[tile_id], angle))
                    z_arr.append(tile_pos[2])
                    id_arr.append(tile_id)
                    del helix
                    gc.collect()

            if i % 1000 == 0 and i != 0:
                print(round((i / n) * 100, 2), "%  |  Hits:", len(z_arr), "  |  Hits without matching trajectory: ",
                      no_traj)
        print("100%")

        self.result_z = np.array(z_arr)
        self.result_angle = np.array(angle_arr)
        self.result_id = np.array(id_arr)

        return self.result_z, self.result_angle, self.result_id

    # ------------------------------------

    def getBinned(self):
        return np.histogram2d(self.result_z, self.result_angle, bins=[220, 180])

    # ------------------------------------

    def saveBinned(self):
        binned_data, xedges, yedges = np.histogram2d(self.result_z, self.result_angle, bins=[220, 180])
        np.savez(self.output + self.angle + self.ana_tpye + self.hit_type + "_binned", data=binned_data, xedges=xedges,
                 yedges=yedges)

    # ------------------------------------

    def getResult(self):
        return self.result_z, self.result_angle, self.result_id

    # ------------------------------------

    def saveTxt(self):
        np.savetxt(self.output + self.angle + self.ana_tpye + self.hit_type + "_z.txt", self.result_z)
        np.savetxt(self.output + self.angle + self.ana_tpye + self.hit_type + "_angle.txt", self.result_angle)
        np.savetxt(self.output + self.angle + self.ana_tpye + self.hit_type + "_id.txt", self.result_id)

    # ------------------------------------

    def saveCompressed(self):
        np.savez_compressed(self.output + self.angle + self.ana_tpye + self.hit_type + "_compressed", z=self.result_z,
                            angle=self.result_angle, id=self.result_id)

    # ------------------------------------

    def saveNpz(self):
        np.savez(self.output + self.angle + self.ana_tpye + self.hit_type, z=self.result_z, angle=self.result_angle,
                 id=self.result_id)
