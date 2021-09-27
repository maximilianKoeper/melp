#---------------------------------------------------------------------
#  HIT CLASS
#       - pos (position)
#       - id
#       - time
#       - edep
#---------------------------------------------------------------------


class Hit():
    def __init__ (self, edep=0, angle=0, mc_i = 0, traj_id = -1, run_id = -1, hid=0):
        self.edep       = edep
        self.mc_i       = mc_i
        self.traj_id    = traj_id

        self.run_id     = run_id

        self.hid        = hid

        self.angle      = angle
