import ROOT
import numpy as np


class TileHitRate:
    def __init__ (self, filename, output_z, output_rate):
        self.filename     = filename
        self.output_z     = output_z
        self.output_rate  = output_rate

        self.result_z     = np.zeros(0)
        self.result_rate  = np.zeros(0)

        self.file         = ROOT.TFile(filename)
        self.mu3e_mchits  = self.file.Get("mu3e_mchits")
        self.mu3e         = self.file.Get("mu3e")
        self.tiles        = self.file.Get("alignment/tiles")

        self.tilehit_tile_dic = {}
        self.tile_mc_i        = {}
        self.tile_id_pos      = {}
        self.tilehit_edep_dic = {}

        self.tilehit_z        = []
        self.tilehit_edep     = []

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
            tile_xyz = []
            tile_xyz.append(self.tiles.posx)
            tile_xyz.append(self.tiles.posy)
            tile_xyz.append(self.tiles.posz)
            self.tile_id_pos[self.tiles.sensor] = tile_xyz

    #####################
    # private functions #
    #####################



    #####################
    # public  functions #
    #####################

    def Tile_Hit_Rate (self, n):
        self.tilehit_z = []
        for i in range(self.mu3e.GetEntries()):
            for j in self.tilehit_tile_dic[i]:
                self.tilehit_z.append(self.tile_id_pos[j][2])

        self.tilehit_edep = []
        for i in range(self.mu3e.GetEntries()):
            for j in self.tilehit_edep_dic[i]:
                self.tilehit_edep.append(j)

        z_arr    = np.array(self.tilehit_z)
        edep_arr = np.array(self.tilehit_edep)

        return z_arr, edep_arr

    # ----
    def getResult(self):
        return np.array(self.tilehit_z), np.array(self.tilehit_edep)
