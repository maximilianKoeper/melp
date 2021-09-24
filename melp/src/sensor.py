#---------------------------------------------------------------------
#  SENSOR CLASS
#       - pos (position)
#       - id  (is automaticly shifted if bigger than 16500)
#       - layer
#       - section
#---------------------------------------------------------------------



import ROOT

import numpy as np

class Sensor():
    def __init__ (self, sensor_pos, sensor_id):
        if sensor_id > 16500:
            sensor_id = sensor_id >> 16
        self.pos    = sensor_pos
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

    @classmethod
    def initRaw (self, pixelid_row_col, ttree_alignment_sensors):
        pixel_id = pixelid_row_col >> 16

        sensor_id_to_index  = {}
        for i in range(ttree_alignment_sensors.GetEntries()):
            ttree_alignment_sensors.GetEntry(i)
            sensor_id_to_index[ttree_alignment_sensors.sensor] = i

        pixel_index = sensor_id_to_index[pixel_id]
        ttree_alignment_sensors.GetEntry(pixel_index)

        row_param = pixelid_row_col & 0xFF
        col_param = (pixelid_row_col >> 8) & 0xFF

        sensor_pos_vxyz = []
        sensor_pos_vxyz.append(ttree_alignment_sensors.vx)
        sensor_pos_vxyz.append(ttree_alignment_sensors.vy)
        sensor_pos_vxyz.append(ttree_alignment_sensors.vz)

        sensor_pos_col = []
        sensor_pos_col.append(ttree_alignment_sensors.colx)
        sensor_pos_col.append(ttree_alignment_sensors.coly)
        sensor_pos_col.append(ttree_alignment_sensors.colz)

        sensor_pos_row = []
        sensor_pos_row.append(ttree_alignment_sensors.rowx)
        sensor_pos_row.append(ttree_alignment_sensors.rowy)
        sensor_pos_row.append(ttree_alignment_sensors.rowz)

        pos = np.array(sensor_pos_vxyz) + (col_param+0.5)*np.array(sensor_pos_col) + (row_param+0.5)*np.array(sensor_pos_row)
        return cls(pos, pixel_id)

    @classmethod
    def initFromAlignment (self, pixel_id, alignment_sensors):
        if pixel_id > 16500:
            pixel_id = pixel_id >> 16

        pos = alignment_sensors[pixel_id]

        return cls(pos, pixel_id)
