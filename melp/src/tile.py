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

        self.hits = []


class TileDet():
    def __init__ (self, tiles):
        self.tile = tiles

    def rate(self):
        hit_rate = np.zeros(len(self.tile))
        index = 0
        for tileID in self.tile:
            for hit in self.tile[tileID].hits:
                hit_rate[index] += 1
            index += 1
        return hit_rate

    def addHit(self, tile, hit):
        self.tile[tile].hits.append(hit)
