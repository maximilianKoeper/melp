#---------------------------------------------------------------------
#  TRAJECTORY CLASS
#       - id
#       - v_pos         (vertex origin)
#       - v_dir         (vertex momentum)
#       - sensor_hits   (saves hits along the trajectory)
#       - tile_hits     (saves hits along the trajectory)
#       - fibre_hits    (saves hits along the trajectory)
#---------------------------------------------------------------------

class Trajectory():
    def __init__ (self, traj_id):
        self.id          = traj_id

        self.sensor_hits = []
        self.tile_hits   = []
        slef.fibre_hits  = []

        self.v_pos       = np.zeros(3)
        slef.v_dir       = np.zeros(3)

    #-----------------------------------------
    #  private functions
    #-----------------------------------------

    #-----------------------------------------
    #  public functions
    #-----------------------------------------

    def add sensorHit(self, sensor):
        self.sensor_hits.append(sensor)

    #-----------------------------------------

    def add tileHit(tile):
        self.tile_hits.append(self, tile)

    #-----------------------------------------

    def add fibreHit(fibre):
        self.fibre_hits.append(self, fibre)
