import ROOT

import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

from melp import Detector
from melp.libs import mathfunctions as mf
from .src.trajectory import Trajectory


# -------------------------
#  VISUALIZER CLASS
# -------------------------
class Visualizer:

    def __init__(self, filename, frame_id=0, trajectories=[]):
        self.trajectories = trajectories
        self.traj_objects = []
        self.frame = frame_id

        self.mu3e_detector = Detector.initFromROOT(filename)

        # -------------------------
        # INIT ROOT File and Branches
        # -------------------------
        self.file = ROOT.TFile(filename)
        self.ttree_mu3e = self.file.Get("mu3e")
        self.ttree_mu3e_mc = self.file.Get("mu3e_mchits")
        self._update_()

        # -------------------------
        # Set default configuration for visuals
        # -------------------------
        self.color_tile = "gray"
        self.color_electrons = "blue"
        self.color_positrons = "red"
        self.f_size = 25

    def __str__(self):
        return "Visualizer info:  Frame:" + self.frame + " Trajectories: " + self.trajectories

    # -------------------------
    # UPDATES Traj information
    # -------------------------
    def _update_(self):
        print("\n -> Fetching Trajectories for new Frame ... ", end='')

        if self.frame >= self.ttree_mu3e.GetEntries():
            raise ValueError('Frame ID > Length of ROOT File')

        self.ttree_mu3e.GetEntry(self.frame)

        # traj ids
        list_traj = list(self.ttree_mu3e.traj_ID)
        list_traj_mother = list(self.ttree_mu3e.traj_mother)

        # vertex position
        list_traj_vx = list(self.ttree_mu3e.traj_vx)
        list_traj_vy = list(self.ttree_mu3e.traj_vy)
        list_traj_vz = list(self.ttree_mu3e.traj_vz)

        # momentum information
        list_traj_px = list(self.ttree_mu3e.traj_px)
        list_traj_py = list(self.ttree_mu3e.traj_py)
        list_traj_pz = list(self.ttree_mu3e.traj_pz)

        # particle types
        list_traj_pid = list(self.ttree_mu3e.traj_type)

        # Tile Hit information
        self.list_tilehit_tile = list(self.ttree_mu3e.tilehit_tile)
        list_tilehit_mci = list(self.ttree_mu3e.tilehit_mc_i)
        self.list_tilehit_tid = []

        for index in list_tilehit_mci:
            self.ttree_mu3e_mc.GetEntry(index)
            self.list_tilehit_tid.append(self.ttree_mu3e_mc.tid)

        self.traj_objects = []
        for i in range(len(list_traj)):
            self.traj_objects.append(Trajectory(id=list_traj[i],
                                                vx=list_traj_vx[i],
                                                vy=list_traj_vy[i],
                                                vz=list_traj_vz[i],
                                                px=list_traj_px[i],
                                                py=list_traj_py[i],
                                                pz=list_traj_pz[i],
                                                pid=list_traj_pid[i],
                                                mother=list_traj_mother[i],
                                                tile_hit_ids=[]))

        self.select_all_trajectories()
        print("done | #TRAJ: ", len(self.traj_objects))

    # -------------------------
    # "CALLABLE" METHODS
    # -------------------------
    def set_frame_id(self, frame_id: int):
        self.frame = frame_id
        self._update_()

    def reset_frame(self):
        self._update_()

    def set_trajectories(self, trajectories: list):
        self.trajectories = trajectories

    def add_toy_event(self, toy_traj: Trajectory):
        self.traj_objects.append(toy_traj)

    def select_all_trajectories(self):
        self.trajectories = []
        for traj in self.traj_objects:
            self.trajectories.append(traj.id)

    def select_all_toy_trajectories(self):
        self.trajectories = []
        for traj in self.traj_objects:
            if traj.mother == -1:
                self.trajectories.append(traj.id)

    def show(self):
        warnings.filterwarnings(action='ignore')
        plt.rcParams.update({'font.size': self.f_size})
        x_t, y_t, z_t = self._get_tile_hit_positions_()

        x_t_geom, y_t_geom, z_t_geom = self._get_tile_detector_positions_()

        fig, ax_arr = plt.subplots(1, 2, figsize=(30, 10))

        ax = ax_arr[0]
        for index_current in self.trajectories:
            #ax.scatter(self.list_traj_vx[index_current], self.list_traj_vy[index_current])
            ax.scatter(np.array(x_t), np.array(y_t), marker="o", color="r", linewidths=3)
            #ax.arrow(self.list_traj_vx[index_current],
            #         self.list_traj_vy[index_current],
            #         self.list_traj_px[index_current],
            #         self.list_traj_py[index_current], length_includes_head=True, head_width=5, head_length=5)
        for traj in self.traj_objects:
            if traj.id not in self.trajectories:
                continue
            x, y, z = traj.get_helix_path()
            color = "gray"
            if traj.get_particle_type() == 1:
                color = "red"
            elif traj.get_particle_type() == 2:
                color = "blue"
            elif traj.get_particle_type() == 0:
                color = "green"
            ax.plot(x, y, color=color)

        ax.scatter(x_t_geom, y_t_geom, marker=".", alpha=0.1, color=self.color_tile)

        circle2 = plt.Circle((0, 0), 19, color='black', alpha=0.1)
        ax.add_patch(circle2)

        ax.axis('equal')
        ax.set_title("xy-plane", fontsize=self.f_size)
        ax.set_ylabel("y", fontsize=self.f_size)
        ax.set_xlabel("x", fontsize=self.f_size)

        ax = ax_arr[1]
        for index_current in self.trajectories:
            #ax.scatter(self.list_traj_vz[index_current], self.list_traj_vx[index_current])
            ax.scatter(np.array(z_t), np.array(x_t), marker="o", color="r", linewidths=3)
            #ax.arrow(self.list_traj_vz[index_current],
            #         self.list_traj_vx[index_current],
            #         self.list_traj_pz[index_current],
            #         self.list_traj_px[index_current], length_includes_head=True, head_width=5, head_length=5)
        for traj in self.traj_objects:
            if traj.id not in self.trajectories:
                continue
            x, y, z = traj.get_helix_path()
            color = "gray"
            if traj.get_particle_type() == 1:
                color = "red"
            elif traj.get_particle_type() == 2:
                color = "blue"
            elif traj.get_particle_type() == 0:
                color = "green"
            ax.plot(z, x, color=color)
        ax.scatter(z_t_geom, x_t_geom, marker=".", alpha=0.1, color=self.color_tile)

        ax.axis('equal')
        ax.set_title("zx-plane", fontsize=self.f_size)
        ax.set_ylabel("y", fontsize=self.f_size)
        ax.set_xlabel("z", fontsize=self.f_size)

        plt.show()

    def show_3d(self):
        x_t, y_t, z_t = self._get_tile_hit_positions_()
        x_t_geom, y_t_geom, z_t_geom = self._get_tile_detector_positions_()

        fig = plt.figure(figsize=plt.figaspect(1) * 3)
        ax = fig.add_subplot(111, projection='3d')

        # plot 2d target (3d not yet working)
        circle2 = plt.Circle((0, 0), 19, color='black', alpha=0.3)
        ax.add_patch(circle2)
        art3d.pathpatch_2d_to_3d(circle2, z=0, zdir="x")

        # add trajectories
        for traj in self.traj_objects:
            if traj.id not in self.trajectories:
                continue
            x, y, z = traj.get_helix_path()
            color = "gray"
            if traj.get_particle_type() == 1:
                color = "red"
            elif traj.get_particle_type() == 2:
                color = "blue"
            elif traj.get_particle_type() == 0:
                color = "green"
            ax.plot(z, x, y, color=color)

        # add tile hits
        ax.scatter(np.array(z_t), np.array(x_t), np.array(y_t), marker="o", color="r", linewidths=3)

        # add tile detector geometry
        ax.scatter(np.array(z_t_geom), np.array(x_t_geom), np.array(y_t_geom), marker=".", alpha=0.1, color=self.color_tile)

        ax.set_xlim(-600, 600)
        ax.set_ylim(-600, 600)
        ax.set_zlim(-600, 600)
        plt.show(block=True)

    def list_traj_info(self):
        for traj in self.traj_objects:
            str_print = ""
            str_print += " | Traj_ID: " + str(traj.id)
            str_print += " | p_z: " + str(np.round(traj.pz, 2))
            str_print += " | p_t: " + str(np.round(traj.get_pt(), 2))
            str_print += " | particle: " + str(traj.get_particle_type())
            str_print += " | mother: " + str(traj.mother) + "\n"
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
        for traj in self.traj_objects:
            if traj.id not in self.trajectories:
                continue

            for index in range(len(self.list_tilehit_tid)):
                if self.list_tilehit_tid[index] == traj.id:
                    x_t.append(self.mu3e_detector.TileDetector.tile[self.list_tilehit_tile[index]].pos[0])
                    y_t.append(self.mu3e_detector.TileDetector.tile[self.list_tilehit_tile[index]].pos[1])
                    z_t.append(self.mu3e_detector.TileDetector.tile[self.list_tilehit_tile[index]].pos[2])

        return x_t, y_t, z_t
