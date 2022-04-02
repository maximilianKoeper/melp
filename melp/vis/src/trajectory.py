import dataclasses

import numpy as np

from melp.libs import mathfunctions as mf


@dataclasses.dataclass
class Trajectory:
    id: int

    vx: float
    vy: float
    vz: float

    px: float
    py: float
    pz: float

    pid: int

    mother: int

    tile_hit_ids: list

    def get_particle_type(self):
        return int(repr(self.pid)[-1])

    def get_pt(self):
        return np.sqrt(self.px ** 2 + self.py ** 2)

    def is_toy_event(self):
        if self.mother == -1:
            return true
        else:
            return false

    def get_helix_path(self):
        x, y, z = [], [], []
        if self.get_particle_type() in [1, 2, 3, 4]:
            x, y, z = self._get_e_helix_path_()
        elif self.get_particle_type() == 0:
            x, y, z = self._get_photon_path_()
        return x, y, z

    # -------------------------
    # "PRIVATE" METHODS
    # -------------------------
    def _get_e_helix_path_(self):
        # calculate radius for helix
        radius = self._get_helix_radius_()

        # calculate theta for helix (angle towards z direction)
        theta = abs(np.arctan2(self.pz, self.get_pt()))
        if self.pz < 0:
            theta *= -1

        # calculate middle point of helices
        vec = np.cross([self.px, self.py, 0], [0, 0, 1])
        vev_normed = vec / np.sqrt(np.sum(vec ** 2))
        hx = self.vx + radius * vev_normed[0]
        hy = self.vy + radius * vev_normed[1]
        hz = self.vz

        # calculating phi offsets
        phi = mf.angle_between_phi(np.array([hx - self.vx, hy - self.vy]), np.array([0, 1]))
        phi_offset = np.deg2rad(-phi + 270)

        # calculating helix path
        if self.get_particle_type() in [3, 4]:  # MU+-
            x, y, z = self._get_helix_(hx, hy, hz, radius, theta, phi_offset, 3)
        elif self.get_particle_type() in [1, 2]:  # e+-
            x, y, z = self._get_helix_(hx, hy, hz, radius, theta, phi_offset, 2)
        return x, y, z

    def _get_photon_path_(self):
        x, y, z = [self.vx], [self.vy], [self.vz]

        x.append(self.vx + 50 * self.px)
        y.append(self.vy + 50 * self.py)
        z.append(self.vz + 50 * self.pz)

        return x, y, z

    def _get_helix_radius_(self, magnetic_field: float = -1.):

        # transversal momentum
        p_t = self.get_pt()

        if self.get_particle_type() == 0:  # photon
            return np.NINF
        elif self.get_particle_type() == 1:  # positron
            charge_q = +1
        elif self.get_particle_type() == 2:  # electron
            charge_q = -1
        elif self.get_particle_type() == 3:  # Mu+
            charge_q = +1
        elif self.get_particle_type() == 4:  # Mu-
            charge_q = -1
        else:
            return np.NaN

        r = p_t / (0.3 * magnetic_field * charge_q)
        return r

    # -------------------------
    # "STATIC" METHODS
    # -------------------------
    @staticmethod
    def _get_helix_(hx: float, hy: float, hz: float, r: float, theta_i: float, phi_offset: float, l: float):
        theta_max = l * np.pi
        theta = np.linspace(0, theta_max, 100)

        x = abs(r) * np.cos(-np.sign(r) * theta + phi_offset) + hx
        y = abs(r) * np.sin(-np.sign(r) * theta + phi_offset) + hy
        z = abs(r) * np.tan(theta_i) * theta + hz

        return x, y, z
