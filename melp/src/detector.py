#---------------------------------------------------------------------
#  DETECTOR CLASS
#       - creates a detecor with tiles and pixel sensors
#  TODO:
#       - add fibre / mmpcs to detecor
#---------------------------------------------------------------------
import ROOT

import pickle

from melp.src.sensor import Sensor
from melp.src.tile import Tile
from melp.src.hit import Hit
from melp.src.tile import TileDet

class Detector():
    def __init__ (self, tiles, sensors, hits = {}):
        #self.Tiles   = tiles
        self.Sensors = sensors

        self.Tiles = TileDet(tiles)
    #-----------------------------------------
    #  Load Detector geometry from Root File
    #-----------------------------------------
    @classmethod
    def initFromROOT (cls, filename):
        file            = ROOT.TFile(filename)
        ttree_sensor    = file.Get("alignment/sensors")
        ttree_tiles     = file.Get("alignment/tiles")

        tile_id_pos      = {}
        tile_id_dir      = {}

        sensor_id_pos    = {}

        for i in range(ttree_tiles.GetEntries()):
            ttree_tiles.GetEntry(i)
            # direction
            xyz = []
            xyz.append(ttree_tiles.dirx)
            xyz.append(ttree_tiles.diry)
            xyz.append(ttree_tiles.dirz)

            tile_id_dir[ttree_tiles.sensor] = xyz

            # position
            tile_xyz = []
            tile_xyz.append(ttree_tiles.posx)
            tile_xyz.append(ttree_tiles.posy)
            tile_xyz.append(ttree_tiles.posz)
            tile_id_pos[ttree_tiles.sensor] = tile_xyz

        Tiles = {}
        for id in tile_id_pos:
            Tiles[id] = Tile(tile_id_pos[id], tile_id_dir[id], id)

        for id in sensor_id_pos:
            #self.Sensors[id] = Sensor(sensor_id_pos[id], id)
            pass

        return cls(Tiles, [0,1])

    #-----------------------------------------
    #  Load Detector geometry from Save File
    #-----------------------------------------
    @classmethod
    def initFromSave (cls, filename):
        data = []
        with open(filename, "rb") as f:
            for i in pickle.load(f):
                data.append(i)

        return cls(data[0],data[1])
    #-----------------------------------------
    #  private functions
    #-----------------------------------------

    #-----------------------------------------
    #  public functions
    #-----------------------------------------

    #-----------------------------------------

    def save(self, filename):
        data = [self.Tiles, self.Sensors]

        with open(filename, "wb") as f:
            pickle.dump(data, f)


    def addTileHits(self, filename):
        file          = ROOT.TFile(filename)
        ttree_mu3e    = file.Get("mu3e")
        #ttree_mu3e_mc = self.file.Get("mu3e_mchits")
        for frame in range(ttree_mu3e.GetEntries()):
            ttree_mu3e.GetEntry(frame)
            for i in range(len(ttree_mu3e.tilehit_tile)):
                tile = ttree_mu3e.tilehit_tile[i]
                edep = ttree_mu3e.tilehit_edep[i]
                mc_i = ttree_mu3e.tilehit_mc_i[i]


                tilehit = Hit(edep = edep, mc_i = mc_i)
                self.Tiles.addHit(tile, tilehit)
