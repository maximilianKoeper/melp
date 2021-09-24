#---------------------------------------------------------------------
#  SENSOR CLASS
#       - pos (position)
#       - id
#       - layer
#       - section
#---------------------------------------------------------------------



import ROOT

import numpy as np

class Sensor():
    def __init__ (self, sensor_pos, sensor_row, sensor_col, sensor_id):
        if sensor_id > 16500:
            sensor_id = sensor_id >> 16
        self.pos    = sensor_pos
        self.row    = sensor_row
        self.col    = sensor_col
        self.id     = sensor_id

        if (sensor_id >= 0 and sensor_id < 1024):
            self.layer   = 0
            self.section = 0
        elif (sensor_id >= 1024 and sensor_id < 2000):
            self.layer   = 1
            self.section = 0
        elif (sensor_id >= 2000 and sensor_id < 3000):
            self.layer   = 2
            self.section = 0
        elif (sensor_id >= 10000 and sensor_id < 11500):
            self.layer   = 2
            self.section = -1
        elif (sensor_id >= 14000 and sensor_id < 15200):
            self.layer   = 2
            self.section = 1
        elif (sensor_id >= 3000 and sensor_id < 4000):
            self.layer   = 3
            self.section = 0
        elif (sensor_id >= 11500 and sensor_id < 12500):
            self.layer   = 3
            self.section = -1
        elif (sensor_id >= 15200 and sensor_id < 16500):
            self.layer   = 3
            self.section = 1


class SensorModul():
    def __init__ (self, sensors):
        self.sensor   = sensors

    #-----------------------------------------
    #  public functions
    #-----------------------------------------

    #-----------------------------------------
    def getPixelPos (self, pixelid_row_col):
        pixel_id = pixelid_row_col >> 16

        row_param = pixelid_row_col & 0xFF
        col_param = (pixelid_row_col >> 8) & 0xFF

        sensor_pos_vxyz = self.sensor[pixel_id].pos
        sensor_pos_col  = self.sensor[pixel_id].col
        sensor_pos_row  = self.sensor[pixel_id].row

        pos = np.array(sensor_pos_vxyz) + (col_param+0.5)*np.array(sensor_pos_col) + (row_param+0.5)*np.array(sensor_pos_row)
        return pos
