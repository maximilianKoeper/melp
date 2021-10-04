# ---------------------------------------------------------------------
#  HIT CLASS
# ---------------------------------------------------------------------


class Hit:
    def __init__(self, edep=0, mc_i=0, traj_id=-1, run_id=-1, hid=0, impact_vec=None, trajectory=None):
        self.edep = edep
        self.mc_i = mc_i
        self.traj_id = traj_id

        self.run_id = run_id

        self.hid = hid

        self.impact_vec = impact_vec

        self.trajectory = trajectory
