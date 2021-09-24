#---------------------------------------------------------------------
#  TILE CLASS
#       - pos (position)
#       - dir (normal vector to tile surface)
#       - id
#---------------------------------------------------------------------
import numpy as np

class Tile():
    def __init__ (self, tile_pos, tile_dir, tile_id):
        self.pos  = tile_pos
        self.dir  = tile_dir
        self.id   = tile_id

        self.hits  = []
        self.angle = []


class TileDetector():
    def __init__ (self, tiles):
        self.tile   = tiles

    def rateId(self):
        hit_rate = {}
        for tileID in self.tile:
            hit_rate[tileID] = 0
            for hit in self.tile[tileID].hits:
                hit_rate[tileID] += 1
        return hit_rate

    def rateZ(self):
        dict_rate = self.rateId()

        z_arr   = []
        hit_arr = []

        for tileID in dict_rate:
            z_arr.append(self.tile[tileID].pos[2])
            hit_arr.append(dict_rate[tileID])

        return z_arr, hit_arr

    def addHit(self, tile, hit):
        self.tile[tile].hits.append(hit)

    def getPos(self, tileID):
        return self.tile[tileID].pos
