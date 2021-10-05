# ---------------------------------------------------------------------
#  TILE CLASS
#       - pos (position)
#       - dir (normal vector to tile surface)
#       - id
# ---------------------------------------------------------------------
import dataclasses


@dataclasses.dataclass
class Tile:
    id: int
    pos: list
    dir: list

    hits: list = dataclasses.field(default_factory=list)


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
