# ---------------------------------------------------------------------
#  DETECTOR CLASS
#       - creates a detector with tiles and pixel sensors
#  TODO:
#       - add fibre / mmpcs to detector
# ---------------------------------------------------------------------
import ROOT

import pickle
import numpy as np

from melp.src.sensor import Sensor
from melp.src.sensor import SensorModule
from melp.src.tile import Tile
from melp.src.tile import TileDetector


class Detector:
    def __init__(self, tiles, sensors, runs=None):
        # self.Tiles   = tiles
        self.SensorsModules = sensors
        self.TileDetector = tiles

        self.info()

    # -----------------------------------------
    #  Load Detector geometry from Root File
    # -----------------------------------------
    @classmethod
    def initFromROOT(cls, filename: str):
        file = ROOT.TFile(filename)
        ttree_sensor = file.Get("alignment/sensors")
        ttree_tiles = file.Get("alignment/tiles")
        tilemisal = False
        if file.Get("alignment/tiles_mis"):
            ttree_tiles_misal = file.Get("alignment/tiles_mis")
            tilemisal = True

        # ---------------
        # TILES
        Tiles = {}
        for i in range(ttree_tiles.GetEntries()):
            kwargs = {}
            ttree_tiles.GetEntry(i)

            kwargs["id"] = ttree_tiles.sensor
            kwargs["pos"] = [ttree_tiles.posx, ttree_tiles.posy, ttree_tiles.posz]
            kwargs["dir"] = [ttree_tiles.dirx, ttree_tiles.diry, ttree_tiles.dirz]
            if tilemisal:
                kwargs["dt"] = ttree_tiles_misal.dt

            Tiles[ttree_tiles.sensor] = Tile(**kwargs)

        # ---------------
        # PIXEL
        Sensors = {}

        for i in range(ttree_sensor.GetEntries()):
            ttree_sensor.GetEntry(i)
            sensor_pos = np.array([ttree_sensor.vx, ttree_sensor.vy, ttree_sensor.vz])
            sensor_row = np.array([ttree_sensor.rowx, ttree_sensor.rowy, ttree_sensor.rowz])
            sensor_col = np.array([ttree_sensor.colx, ttree_sensor.coly, ttree_sensor.colz])
            Sensors[ttree_sensor.sensor] = Sensor(sensor_pos, sensor_row, sensor_col, ttree_sensor.sensor)
            pass

        return cls(TileDetector(Tiles, tilemisal), SensorModule(Sensors))

    def __str__(self):
        return f'Detector(TileDetector={self.TileDetector}, SensorModules={self.SensorsModules}))'

    # -----------------------------------------
    #  Load Detector geometry from Save File
    # -----------------------------------------
    @classmethod
    def initFromSave(cls, filename: str):
        data = []
        with open(filename, "rb") as f:
            for i in pickle.load(f):
                data.append(i)

        return cls(data[0], data[1])

    # -----------------------------------------
    #  private functions
    # -----------------------------------------

    # -----------------------------------------
    #  public functions
    # -----------------------------------------

    def info(self):
        print("------------------------------")
        print("Detector information\n")
        print("Stats:")
        print("  - Tiles: ", len(self.TileDetector.tile))
        print("    -> misal: ", self.TileDetector.tilemisal)
        print("  - Pixel Modules: ", len(self.SensorsModules.sensor))
        print("  - Loaded Runs (Tiles): ", self.TileDetector.AddedRuns)
        print("  - Loaded Runs (Pixel): ", self.SensorsModules.AddedRuns)
        print("------------------------------")

    # -----------------------------------------

    def save(self, filename: str):
        data = [self.TileDetector, self.SensorsModules]

        with open(filename, "wb") as f:
            pickle.dump(data, f)
