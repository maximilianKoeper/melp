#---------------------------------------------------------------------
#  TILE CLASS
#       - pos (position)
#       - dir (normal vector to tile surface)
#       - id
#---------------------------------------------------------------------
import ROOT

import numpy as np

from melp.libs import mathfunctions as mf

class Tile():
    def __init__ (self, tile_pos, tile_dir, tile_id):
        self.pos  = tile_pos
        self.dir  = tile_dir
        self.id   = tile_id

        self.hits       = []
        self.impact_vec = []


class TileDetector():
    def __init__ (self, tiles):
        self.tile   = tiles

        self.result_z     = np.zeros(0)
        self.result_angle = np.zeros(0)

    #-----------------------------------------
    #  public functions
    #-----------------------------------------

    #-----------------------------------------
    def rateId(self):
        hit_rate = {}
        for tileID in self.tile:
            hit_rate[tileID] = 0
            for hit in self.tile[tileID].hits:
                hit_rate[tileID] += 1
        return hit_rate

    #-----------------------------------------

    def rateZ(self):
        dict_rate = self.rateId()

        z_arr   = []
        hit_arr = []

        for tileID in dict_rate:
            z_arr.append(self.tile[tileID].pos[2])
            hit_arr.append(dict_rate[tileID])

        return z_arr, hit_arr

    #-----------------------------------------

    def addHit(self, tile, hit):
        self.tile[tile].hits.append(hit)

    #-----------------------------------------

    def getPos(self, tileID):
        return self.tile[tileID].pos

    #-----------------------------------------
    def calcTruthImpactVec(self, filename):

        file          = ROOT.TFile(filename)
        ttree_mu3e    = file.Get("mu3e")
        ttree_mu3e_mc = file.Get("mu3e_mchits")

        for frame in range(ttree_mu3e.GetEntries()):
            ttree_mu3e.GetEntry(frame)
            for tilehit in range(len(ttree_mu3e.tilehit_tile)):
                tileId   = ttree_mu3e.tilehit_tile[tilehit]
                tile_pos = self.tile[tileId]

                tile_mc_i = ttree_mu3e.tilehit_mc_i[tilehit]

                ttree_mu3e_mc.GetEntry(tile_mc_i)

                px   = ttree_mu3e_mc.p_in_x
                py   = ttree_mu3e_mc.p_in_y
                pz   = ttree_mu3e_mc.p_in_z

                p_xyz    = np.zeros(3)
                p_xyz[0] = px
                p_xyz[1] = py
                p_xyz[2] = pz

                self.tile[tileId].impact_vec.append(p_xyz)

    #-----------------------------------------
    def calcAngleTruth_byId(self, angle="theta"):
        hit_angle = {}

        for tileID in self.tile:
            hit_angle_tmp= []
            for vec in self.tile[tileID].impact_vec:
                hit_angle_tmp.append(mf.angle_between(vec, np.array([0,0,-1])))
            hit_angle[tileID] = hit_angle_tmp
        return hit_angle

    #-----------------------------------------
    def calcAngleTruth_byZ(self, angle="theta"):
        hit_angle = self.calcAngleTruth_byId(angle)

        z_arr     = []
        angle_arr = []

        for tileID in hit_angle:
            for hit in hit_angle[tileID]:
                z_arr.append(self.tile[tileID].pos[2])
                angle_arr.append(hit)

        self.result_z     = np.array(z_arr)
        self.result_angle = np.array(angle_arr)
        return self.result_z, self.result_angle

    #-----------------------------------------
    def getBinned(self):
        return np.histogram2d(self.result_z, self.result_angle, bins=[220,180])
