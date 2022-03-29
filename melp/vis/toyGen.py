import numpy as np
from .src.trajectory import Trajectory


class ToyEventGenerator:
    def __init__(self, pt: float, pz: float, particle_type: int):
        self.id = 0
        self.pt = pt
        self.pz = pz
        self.particle_type = particle_type

    def new_event(self):
        vx = 0.
        vy = 0.
        vz = 0.

        theta = np.random.uniform(0, 2*np.pi)
        px = self.pt * np.cos(theta)
        py = self.pt * np.sin(theta)

        pid = self.particle_type
        self.id += 1
        return Trajectory(id=self.id,
                          vx=vx,
                          vy=vy,
                          vz=vz,
                          px=px,
                          py=py,
                          pz=self.pz,
                          pid=pid,
                          mother=-1,
                          tile_hit_ids=[])

    def _get_charge_(self):
        if self.particle_type == 0:  # photon
            charge_q = np.NINF
        elif self.particle_type == 1:  # positron
            charge_q = +1
        elif self.particle_type == 2:  # electron
            charge_q = -1
        else:
            charge_q = np.NaN

        return charge_q
