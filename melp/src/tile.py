# ---------------------------------------------------------------------
#  TILES and TILE DETECTOR
# ---------------------------------------------------------------------
import dataclasses
import warnings
import math
import random
from functools import lru_cache


# ---------------------------------------------------------------------
#  TILE CLASS
#       - pos (position)
#       - dir (normal vector to tile surface)
#       - id
#       - dt (time misal truth)
#       - dt_cal
# ---------------------------------------------------------------------
@dataclasses.dataclass
class Tile:
    id: int
    pos: list
    dir: list
    dt_truth: float = 0.
    dt_cal: float = 0.
    dt_truth_abs = 0.
    dt_cal_abs = 0.

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
        print("  - Truth Time Misal: ", self.dt_truth)
        print("  - Calibrated Time Misal: ", self.dt_cal)
        print("------------------------------")

    def row(self) -> int:
        tmp_id = self.id - 200000
        if self.station() == 2:
            tmp_id -= 100000

        res = tmp_id % 56
        return res

    def column(self) -> int:
        tmp_id = self.id - 200000
        if self.station() == 2:
            tmp_id -= 100000

        tmp = math.floor(tmp_id / 56)
        return tmp

    def station(self) -> int:
        if self.id >= 300000:
            return 2
        elif self.id < 300000:
            return 1

    def get_truth_offset(self) -> float:
        return self.dt_truth

    def get_calibrated_offset(self) -> float:
        return self.dt_cal

    def get_offset(self) -> float:
        return self.dt_truth - self.dt_cal

    def update_calibration(self, offset):
        self.dt_cal += offset

    def __eq__(self, other):
        return self.id == other.id


# ---------------------------------------------------------------------
#  TILE Detector CLASS
# ---------------------------------------------------------------------
class TileDetector:
    def __init__(self, tiles: dict, misal=False):
        self.tile = tiles
        self.AddedRuns = []
        self.hitrate = []
        self.hitangle = []
        self.tilemisal = misal
        self.calibrated = False

    def __str__(self):
        return f'Loaded Tiles: {len(self.tile)}'
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
        self.tile[tile].dt_truth = dt

    # -----------------------------------------
    # returns neighbour id for given tile
    # return False if it doesnt exist
    #
    # left: -56
    # right: +56
    # up: +1 (exception: end of ring)
    # down: -1 (exception beginning of ring)
    #
    @lru_cache
    def getNeighbour(self, tileid: int, position: str) -> int or bool:

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

    # -----------------------------------------
    # z - direction
    @lru_cache
    def row_ids(self, phi: int, station_offset: int) -> list:
        if 0 <= phi <= 55:
            tile_ids = []
            for tile in range(0, 52):
                tile_id = station_offset + phi + tile * 56
                tile_ids.append(tile_id)

            if set(tile_ids).issubset(set(self.tile.keys())) is not True:
                warnings.warn("Tile IDs do not match loaded geometry")

            return tile_ids
        else:
            raise ValueError("Row must be between 0-55")

    # -----------------------------------------
    # phi - direction
    @lru_cache
    def column_ids(self, z: int, station_offset: int) -> list:
        if 0 <= z <= 51:
            tile_ids = []
            for tile in range(0, 56):
                tile_id = station_offset + tile + z * 56
                tile_ids.append(tile_id)

            if set(tile_ids).issubset(set(self.tile.keys())) is not True:
                warnings.warn("Tile IDs do not match loaded geometry")

            return tile_ids
        else:
            raise ValueError("Column must be between 0-51")

    # -----------------------------------------
    # get tile id from row and column
    @lru_cache
    def id_from_row_col(self, row: int, column: int, station_offset: int) -> int:
        tile_id = station_offset + row + column * 56

        if tile_id not in self.tile.keys():
            raise ValueError(f"No tile at given coordinates {row}, {column}")

        if row > 55 or column > 51:
            raise ValueError("expected row < 56 and column < 52")

        return tile_id

    # -----------------------------------------
    # generate Time Misalignment
    def generate_misal(self, **kwargs):
        for tileid in self.tile:
            self.tile[tileid].dt_truth = random.gauss(mu=0.0, sigma=kwargs["sigma"])

        self.tilemisal = True
