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
    dt: float = 0.

    hits: list = dataclasses.field(default_factory=list)

    # -----------------------------------------
    #  public functions
    # -----------------------------------------

    def info(self):
        print("------------------------------")
        print("Tile information\n")
        print("  - Tile ID: ", self.id)
        print("  - Position: ", self.pos)
        print("  - Direction: ", self.dir)
        print("  - Total Hits: ", len(self.hits))
        print("  - Truth Time Misal: ", self.dt)
        print("------------------------------")


class TileDetector:
    def __init__(self, tiles: dict, misal=False):
        self.tile = tiles
        self.AddedRuns = []
        self.hitrate = []
        self.hitangle = []
        self.tilemisal = misal

    # -----------------------------------------
    #  public functions
    # -----------------------------------------

    def addHit(self, tile, hit):
        self.tile[tile].hits.append(hit)

    # -----------------------------------------

    def getPos(self, tileID):
        return self.tile[tileID].pos

    # -----------------------------------------

    def addRateResult(self, hitrate_tmp):
        self.hitrate.append(hitrate_tmp)

    # -----------------------------------------

    def addAngleResult(self, hitangle):
        self.hitangle.append(hitangle)
    # -----------------------------------------

    def addDT(self, tile: int, dt: float):
        self.tile[tile].dt = dt
