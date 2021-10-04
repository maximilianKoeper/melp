# ---------------------------------------------------------------------
#  TILE CLASS
#       - pos (position)
#       - dir (normal vector to tile surface)
#       - id
# ---------------------------------------------------------------------

class Tile:
    def __init__(self, tile_pos, tile_dir, tile_id):
        self.pos = tile_pos
        self.dir = tile_dir
        self.id = tile_id

        self.hits = []
        # self.impact_vec = []


class TileDetector:
    def __init__(self, tiles: dict):
        self.tile = tiles

    # -----------------------------------------
    #  public functions
    # -----------------------------------------

    def addHit(self, tile, hit):
        self.tile[tile].hits.append(hit)

    # -----------------------------------------

    def getPos(self, tileID):
        return self.tile[tileID].pos

    # -----------------------------------------

    def addRateResult(self, hitrate):
        self.hitrate = hitrate

    # -----------------------------------------

    def addAngleResult(self, hitangle):
        self.hitangle = hitangle
    # -----------------------------------------
