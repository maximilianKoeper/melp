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

    # -----------------------------------------
    # returns neighbour id for given tile
    # return False if it doesnt exist
    #
    # left: -56
    # right: +56
    # up: +1 (exception: end of ring)
    # down: -1 (exception beginning of ring)
    #
    def getNeighbour(self, tileid: int, position: str):

        # left
        if position == "l" or position == "left":
            neighbour_id = tileid - 56
            if neighbour_id not in self.tile.keys():
                return False
            else:
                return neighbour_id

        # right
        elif position == "r" or position == "right":
            neighbour_id = tileid + 56
            if neighbour_id not in self.tile.keys():
                return False
            else:
                return neighbour_id

        # down
        elif position == "d" or position == "down":
            # -----------------
            # checking for 56th tile (56 +1)th tile is not in the same ring anymore !!!
            tmp_hittile = tileid - 200000
            if tmp_hittile >= 100000:
                tmp_hittile -= 100000
            if tmp_hittile % 56 == 0 or tmp_hittile == 0:
                neighbour_id = tileid + 55
            else:
                neighbour_id = tileid - 1
            # -----------------
            if neighbour_id not in self.tile.keys():
                return False

            return neighbour_id

        # up
        elif position == "u" or position == "up":
            # -----------------
            # checking for 56th tile (56 +1)th tile is not in the same ring anymore !!!
            tmp_hittile = tileid - 200000
            if tmp_hittile >= 100000:
                tmp_hittile -= 100000
            if (tmp_hittile + 1) % 56 == 0:
                neighbour_id = tileid - 55
            else:
                neighbour_id = tileid + 1
            # -----------------
            if neighbour_id not in self.tile.keys():
                return False

            return neighbour_id

        else:
            raise ValueError(f"getNeighbour(tileid, position=(l,r,d,u)) expected, got: {position}")
