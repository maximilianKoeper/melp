import ROOT

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

from melp import Detector
from melp.libs import mathfunctions as mf


# -------------------------
#  VISUALIZER CLASS
# -------------------------
class Visualizer:

    def __init__(self, filename, frame_id=0, trajectories=[]):
        self.frame = frame_id
        self.trajectories = trajectories
        # self.fig, self.ax_arr = plt.subplots(1, 2, figsize=(20, 10))

        self.mu3e_detector = Detector.initFromROOT(filename)

        # -------------------------
        # INIT ROOT File and Branches
        # -------------------------
        self.file = ROOT.TFile(filename)
        self.ttree_mu3e = self.file.Get("mu3e")
        self.ttree_mu3e_mc = self.file.Get("mu3e_mchits")
        self.__update__()

    def __str__(self):
        return "Visualizer info:  Frame:" + self.frame + " Trajectories: " + self.trajectories

    # -------------------------
    # UPDATES Traj information
    # -------------------------
    def __update__(self):
        print("\n -> Fetching Trajectories for new Frame ... ", end='')

        if self.frame >= self.ttree_mu3e.GetEntries():
            raise ValueError('Frame ID > Length of ROOT File')

        self.ttree_mu3e.GetEntry(self.frame)

        # traj ids
        self.list_traj = list(self.ttree_mu3e.traj_ID)
        self.list_traj_mother = list(self.ttree_mu3e.traj_mother)

        # vertex position
        self.list_traj_vx = list(self.ttree_mu3e.traj_vx)
        self.list_traj_vy = list(self.ttree_mu3e.traj_vy)
        self.list_traj_vz = list(self.ttree_mu3e.traj_vz)

        # momentum information
        self.list_traj_px = list(self.ttree_mu3e.traj_px)
        self.list_traj_py = list(self.ttree_mu3e.traj_py)
        self.list_traj_pz = list(self.ttree_mu3e.traj_pz)

        # particle types
        self.list_traj_pid = list(self.ttree_mu3e.traj_type)
        self.list_traj_particle_type = []
        for i in range(len(self.list_traj)):
            self.list_traj_particle_type.append(int(repr(self.list_traj_pid[i])[-1]))

        # Tile Hit information
        self.list_tilehit_tile = list(self.ttree_mu3e.tilehit_tile)
        self.list_tilehit_mci = list(self.ttree_mu3e.tilehit_mc_i)
        self.list_tilehit_tid = []

        for index in self.list_tilehit_mci:
            self.ttree_mu3e_mc.GetEntry(index)
            self.list_tilehit_tid.append(self.ttree_mu3e_mc.tid)

        print("done | #TRAJ: ", len(self.list_traj))

    # -------------------------
    # "CALLABLE" METHODS
    # -------------------------
    def set_frame_id(self, frame_id: int):
        self.frame = frame_id
        self.__update__()

    def set_trajectories(self, trajectories: list):
        self.trajectories = trajectories

    def show(self):
        x, y, z = self._calculate_helix_path_()
        x_t, y_t, z_t = self._get_tile_hit_positions_()

        x_T, y_T, z_T = self._get_tile_detector_positions_()

        fig, ax_arr = plt.subplots(1, 2, figsize=(30, 10))

        ax = ax_arr[0]
        for index_current in self.trajectories:
            ax.scatter(self.list_traj_vx[index_current], self.list_traj_vy[index_current])
            ax.scatter(np.array(x_t), np.array(y_t), marker="o", color="r", linewidths=3)
            ax.arrow(self.list_traj_vx[index_current], self.list_traj_vy[index_current], self.list_traj_px[index_current],
                     self.list_traj_py[index_current], length_includes_head=True, head_width=5, head_length=5)
        for i in range(len(x)):
            ax.plot(x[i], y[i])
        ax.scatter(x_T, y_T, marker=".", alpha=0.1, color="g")

        circle2 = plt.Circle((0, 0), 19, color='black', alpha=0.1)
        ax.add_patch(circle2)
        ax.axis('equal')

        ax = ax_arr[1]
        for index_current in self.trajectories:
            ax.scatter(self.list_traj_vz[index_current], self.list_traj_vx[index_current])
            ax.scatter(np.array(z_t), np.array(x_t), marker="o", color="r", linewidths=3)
            ax.arrow(self.list_traj_vz[index_current], self.list_traj_vx[index_current], self.list_traj_pz[index_current],
                     self.list_traj_px[index_current], length_includes_head=True, head_width=5, head_length=5)
        for i in range(len(x)):
            ax.plot(z[i], x[i])
        ax.scatter(z_T, x_T, marker=".", alpha=0.1, color="g")

        ax.axis('equal')

        plt.show()

    def show_3d(self):
        #matplotlib.use('WebAgg')
        x, y, z = self._calculate_helix_path_()
        x_t, y_t, z_t = self._get_tile_hit_positions_()
        x_T, y_T, z_T = self._get_tile_detector_positions_()

        fig = plt.figure(figsize=plt.figaspect(1))
        ax = fig.add_subplot(111, projection='3d')

        circle2 = plt.Circle((0, 0), 19, color='black', alpha=0.3)
        ax.add_patch(circle2)
        art3d.pathpatch_2d_to_3d(circle2, z=0, zdir="x")
        for i in range(len(x)):
            ax.plot(z[i], x[i], y[i])
        ax.scatter(np.array(z_t), np.array(x_t), np.array(y_t), marker="o", color="r", linewidths=3)
        ax.scatter(np.array(z_T), np.array(x_T), np.array(y_T), marker=".", alpha=0.1, color="g")
        ax.set_xlim(-600, 600)
        ax.set_ylim(-600, 600)
        ax.set_zlim(-600, 600)
        plt.show(block=True)

    def list_traj_info(self):
        for i in range(len(self.list_traj)):
            str_print = ""
            str_print += str(i) + " | Traj_ID: " + str(self.list_traj[i])
            str_print += " | p_z: " + str(np.round(self.list_traj_pz[i], 2))
            str_print += " | particle: " + str(self.list_traj_particle_type[i])
            str_print += " | mother: " + str(self.list_traj_mother[i]) + "\n"
            print(str_print)

    # -------------------------
    # "PRIVATE" METHODS
    # -------------------------
    def _get_tile_detector_positions_(self):
        x, y, z = [], [], []
        for tileid in self.mu3e_detector.TileDetector.tile:
            x.append(self.mu3e_detector.TileDetector.tile[tileid].pos[0])
            y.append(self.mu3e_detector.TileDetector.tile[tileid].pos[1])
            z.append(self.mu3e_detector.TileDetector.tile[tileid].pos[2])
        return x, y, z

    def _get_tile_hit_positions_(self):
        x_t = []
        y_t = []
        z_t = []
        for index_current in self.trajectories:
            for index in range(len(self.list_tilehit_tid)):
                if self.list_tilehit_tid[index] == self.list_traj[index_current]:
                    x_t.append(self.mu3e_detector.TileDetector.tile[self.list_tilehit_tile[index]].pos[0])
                    y_t.append(self.mu3e_detector.TileDetector.tile[self.list_tilehit_tile[index]].pos[1])
                    z_t.append(self.mu3e_detector.TileDetector.tile[self.list_tilehit_tile[index]].pos[2])

        return x_t, y_t, z_t

    def _calculate_helix_path_(self):
        # calculate radius for helices
        list_traj_r = []
        for i in self.trajectories:
            list_traj_r.append(self._get_helix_radius_(self.list_traj_px[i], self.list_traj_py[i], self.list_traj_particle_type[i], -1.))

        # calculate theta for helices
        list_traj_theta = []
        for i in self.trajectories:
            px = self.list_traj_px[i]
            py = self.list_traj_py[i]
            pz = self.list_traj_pz[i]

            sign = 1
            if self.list_traj_particle_type[i] == 2:
                sign *= -1
            list_traj_theta.append(sign*mf.angle_between(np.array([px, py, pz]), np.array([0, 0, 1]))-90)

        # calculate middle point of helices
        hx, hy, hz = [], [], []
        for i in range(len(self.trajectories)):
            index_traj = self.trajectories[i]
            v = np.cross([self.list_traj_px[index_traj], self.list_traj_py[index_traj], 0], [0, 0, 1])

            v_n = v / np.sqrt(np.sum(v ** 2))
            r = list_traj_r[i]

            hx.append(self.list_traj_vx[index_traj] + r * v_n[0])
            hy.append(self.list_traj_vy[index_traj] + r * v_n[1])
            hz.append(self.list_traj_vz[index_traj])

        # calculating phi offsets
        phi_offsets = []
        for i in range(len(self.trajectories)):
            index_traj = self.trajectories[i]
            phi = mf.angle_between_phi(np.array([hx[i]-self.list_traj_vx[index_traj], hy[i]-self.list_traj_vy[index_traj]]), np.array([0, 1]))
            phi_offsets.append(np.deg2rad(-phi+270))

        # calculating helix path
        x, y, z = {}, {}, {}
        for i in range(len(self.trajectories)):
            x_tmp, y_tmp, z_tmp = self._get_helix_(hx[i], hy[i],  hz[i], list_traj_r[i], list_traj_theta[i], phi_offsets[i], 2)
            x[i] = x_tmp
            y[i] = y_tmp
            z[i] = z_tmp
        return x, y, z

    # -------------------------
    # "STATIC" METHODS
    # -------------------------
    @staticmethod
    def _get_helix_radius_(px: float, py: float, particle_type: int, magnetic_field: float):

        # transversal momentum
        p_t = np.sqrt(px ** 2 + py ** 2)

        if particle_type == 0:  # photon
            return np.NINF
        elif particle_type == 1:  # positron
            charge_q = +1
        elif particle_type == 2:  # electron
            charge_q = -1
        else:
            return np.NaN

        r = p_t / (0.3 * magnetic_field * charge_q)
        return r

    @staticmethod
    def _get_helix_(hx: float, hy: float, hz: float, r: float, theta_i: float, phi_offset: float, l: float):
        theta_max = l * np.pi
        theta = np.linspace(0, theta_max, 100)

        x = abs(r) * np.cos(-np.sign(r)*theta + phi_offset) + hx
        y = abs(r) * np.sin(-np.sign(r)*theta + phi_offset) + hy
        z = 2 * np.pi * r * np.tan(np.pi * theta_i / 180) * theta / (2 * np.pi) + hz

        return x, y, z
