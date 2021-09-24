#---------------------------------------------------------------------
#  DETECTOR CLASS
#---------------------------------------------------------------------

from melp.src.tile import Tile

class Detector():
    def __init__ (self, tile_id_pos, tile_id_dir, sensor_id_pos):
        self.Tiles   = {}
        self.Sensors = {}
        for id in tile_id_pos:
            self.Tiles[id] = Tile(tile_id_pos[id], tile_id_dir[id], id)

        for id in sensor_id_pos:
            self.Sensors[id] = Sensor(sensor_id_pos[id], id)
