import ROOT
import numpy as np


class TileHitRate:
    def __init__(self, filename, output_z, output_rate):
        self.filename = filename
        self.output_z = output_z
        self.output_rate = output_rate

        self.result_z = np.zeros(0)
        self.result_rate = np.zeros(0)

        self.file = ROOT.TFile(filename)
        self.mu3e_mchits = self.file.Get("mu3e_mchits")
        self.mu3e = self.file.Get("mu3e")
        self.tiles = self.file.Get("alignment/tiles")

        self.tilehit_tile_dic = {}
        self.tile_mc_i = {}
        self.tile_id_pos = {}
        self.tilehit_edep_dic = {}

        self.tilehit_z = []
        self.tilehit_edep = []

        self.z_total_arr = []
        self.z_primary_arr = []
        self.z_secondary_arr = []
        self.z_tertiary_arr = []

        self.edep_total_arr = []
        self.edep_primary_arr = []
        self.edep_secondary_arr = []
        self.edep_tertiary_arr = []

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

            temp_arr_3 = []
            for j in self.mu3e.tilehit_edep:
                temp_arr_3.append(j)
            self.tilehit_edep_dic[i] = temp_arr_3

        for i in range(self.tiles.GetEntries()):
            self.tiles.GetEntry(i)

            # position
            tile_xyz = [self.tiles.posx, self.tiles.posy, self.tiles.posz]
            self.tile_id_pos[self.tiles.sensor] = tile_xyz

    #####################
    # private functions #
    #####################
    def __Get_HID_from_MC_I(self, mc_i):
        self.mu3e_mchits.GetEntry(mc_i)
        hid = self.mu3e_mchits.hid

        return hid

    #####################
    # public  functions #
    #####################

    def tileHitRate(self, n):
        self.tilehit_z = []
        for i in range(self.mu3e.GetEntries()):
            for j in self.tilehit_tile_dic[i]:
                self.tilehit_z.append(self.tile_id_pos[j][2])

        self.tilehit_edep = []
        for i in range(self.mu3e.GetEntries()):
            for j in self.tilehit_edep_dic[i]:
                self.tilehit_edep.append(j)

        z_arr = np.array(self.tilehit_z)
        edep_arr = np.array(self.tilehit_edep)

        return z_arr, edep_arr

    # ----
    def tileHitRateHID(self, n=0):

        tilehit_z_total = []
        tilehit_z_primary = []
        tilehit_z_secondary = []
        tilehit_z_tertiary = []

        tilehit_edep_total = []
        tilehit_edep_primary = []
        tilehit_edep_secondary = []
        tilehit_edep_tertiary = []

        # Check Argument
        if n > len(self.tilehit_tile_dic) or n == 0:
            n = len(self.tilehit_tile_dic)
        print(n, " of ", len(self.tilehit_tile_dic))

        # loop over all Root frames
        for i in range(n):

            # loop over all tile hits in one Root frame
            for u in range(len(self.tilehit_tile_dic[i])):

                tile_id = self.tilehit_tile_dic[i][u]
                hid_test = self.__Get_HID_from_MC_I(self.tile_mc_i[i][u])
                j = self.tilehit_tile_dic[i][u]
                k = self.tilehit_edep_dic[i][u]
                tilehit_z_total.append(self.tile_id_pos[j][2])
                tilehit_edep_total.append(k)
                if np.abs(hid_test) == 1:
                    tilehit_z_primary.append(self.tile_id_pos[j][2])
                    tilehit_edep_primary.append(k)
                elif np.abs(hid_test) == 2:
                    tilehit_z_secondary.append(self.tile_id_pos[j][2])
                    tilehit_edep_secondary.append(k)
                elif np.abs(hid_test) == 3:
                    tilehit_z_tertiary.append(self.tile_id_pos[j][2])
                    tilehit_edep_tertiary.append(k)

            # Print progress
            if i % 1000 == 0 and i != 0:
                print(round((i / n) * 100, 2), "%")
        print("100%")

        self.z_total_arr = np.array(tilehit_z_total)
        self.z_primary_arr = np.array(tilehit_z_primary)
        self.z_secondary_arr = np.array(tilehit_z_secondary)
        self.z_tertiary_arr = np.array(tilehit_z_tertiary)

        self.edep_total_arr = np.array(tilehit_edep_total)
        self.edep_primary_arr = np.array(tilehit_edep_primary)
        self.edep_secondary_arr = np.array(tilehit_edep_secondary)
        self.edep_tertiary_arr = np.array(tilehit_edep_tertiary)

        return self.z_total_arr, self.z_primary_arr, self.z_secondary_arr, self.z_tertiary_arr, self.edep_total_arr, self.edep_primary_arr, self.edep_secondary_arr, self.edep_tertiary_arr

        # ----

    def getResult(self):
        return self.z_total_arr, self.z_primary_arr, self.z_secondary_arr, self.z_tertiary_arr, self.edep_total_arr, self.edep_primary_arr, self.edep_secondary_arr, self.edep_tertiary_arr

        # ----

    def saveNpz(self):
        np.savez(self.output_z, total=self.z_total_arr, primary=self.z_primary_arr, secondary=self.z_secondary_arr,
                 tertiary=self.z_tertiary_arr)
        np.savez(self.output_rate, total=self.edep_total_arr, primary=self.edep_primary_arr,
                 secondary=self.edep_secondary_arr, tertiary=self.edep_tertiary_arr)
