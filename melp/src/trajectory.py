# ---------------------------------------------------------------------
#  TRAJECTORY CLASS
# ---------------------------------------------------------------------

class Trajectory:
    def __init__(self, traj_id: int, v_pos, v_dir):
        self.id = traj_id

        self.v_pos = v_pos
        self.v_dir = v_dir
