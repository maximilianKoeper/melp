import numpy as np
from scipy.optimize import minimize, brute

from melp.libs import mathfunctions as mf


# DEBUG: Delete after testing is finished
# import matplotlib.pyplot as plt

class Helices:
    def __init__(self, vx, vy, vz, px, py, pz, htype, tile_pos):

        self.bfield = -1

        self.z0 = vz
        self.type = htype

        self.tile_pos = tile_pos

        pt = np.hypot(px, py)

        self.r = pt / (0.3 * self.bfield)
        self.theta = np.arctan2(pz, pt)
        self.phi = np.arctan2(py, px)

        if self.type == 1:
            self.phi += np.pi / 2
        elif self.type == 2:
            self.phi -= np.pi / 2
        else:
            raise ValueError('Helices: init: type not supported')

        self.dz = 2 * np.pi * self.r * np.tan(self.theta)
        self.xc = self.r * np.cos(self.phi) + vx
        self.yc = self.r * np.sin(self.phi) + vy

        # self.xy_vec     = 0.

        # self.test1      = 0.

    #####################
    # private functions #
    #####################

    def __Helix(self, alpha):
        xyz = np.zeros(3)
        xyz[0] = self.xc - self.r * np.cos(alpha)
        xyz[1] = self.yc - self.r * np.sin(alpha)
        xyz[2] = self.z0 + self.dz * ((alpha - self.phi) / (2 * np.pi))

        return xyz

    # ------------------------------------

    def __Minimize_Func_Angle(self, alpha):
        xyz_circ = self.__Helix(alpha)
        xyz_tile = np.zeros(3)
        xyz_tile[0] = self.tile_pos[0]
        xyz_tile[1] = self.tile_pos[1]
        xyz_tile[2] = self.tile_pos[2]
        # distance    = mf.distance_between_2d(xyz_tile[0:2], xyz_circ[0:2])

        distance = mf.distance_between_3d(xyz_tile, xyz_circ)

        return distance

    # ------------------------------------

    def __Get_Phi(self):
        tmp_min = brute(self.__Minimize_Func_Angle, ranges=((-10 * np.pi, +10 * np.pi),), Ns=100)[0]
        tmp_min = minimize(self.__Minimize_Func_Angle, tmp_min).x

        return tmp_min

    # ------------------------------------

    def __Get_Primary_Tile_Hit_Vector(self):
        temp_phi = self.__Get_Phi()
        v1_tmp = self.__Helix(temp_phi)

        if self.type == 2:
            offset = - 0.1
        else:
            offset = + 0.1

        v2_tmp = self.__Helix(temp_phi + offset)
        xy_hit_vector = np.array(v1_tmp) - np.array(v2_tmp)
        return xy_hit_vector

    # ------------------------------------

    def __Get_Primary_Tile_Hit_Angle(self, tile_norm_vec, angle):
        if angle == "phi":
            norm_vec = -np.array(tile_norm_vec)[0:2]
            temp_vec = mf.angle_between_phi(self.__Get_Primary_Tile_Hit_Vector()[0:2], norm_vec)
            return temp_vec
        elif angle == "theta":
            norm_vec = -np.array([0, 0, 1])
            temp_vec = mf.angle_between(self.__Get_Primary_Tile_Hit_Vector(), norm_vec)
            return temp_vec
        elif angle == "norm":
            norm_vec = tile_norm_vec
            temp_vec = mf.angle_between(self.__Get_Primary_Tile_Hit_Vector(), norm_vec)
            return temp_vec
        else:
            raise ValueError("angle != [phi/theta/norm]")

    #####################
    # public  functions #
    #####################

    def hitAngle(self, tile_norm_vec, angle="phi"):
        return self.__Get_Primary_Tile_Hit_Angle(tile_norm_vec, angle)

    #####################
    # TESTING functions #
    #####################


"""
    def test (self):
        #self.test1 = minimize(self.__Minimize_Func_Angle, 0).x
        self.test1 = brute(self.__Minimize_Func_Angle, ranges=((-4*np.pi,+4*np.pi),), Ns=100)[0]
        self.test1 = minimize(self.__Minimize_Func_Angle, self.test1).x
        print(self.test1)

    def plottest (self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        #startpoint
        ax.scatter(self.__Helix(self.phi)[0],self.__Helix(self.phi)[1], self.__Helix(self.phi)[2], color="r")
        i=self.phi
        while i > -20:
            #trajectory
            ax.scatter(self.__Helix(i)[0],self.__Helix(i)[1],self.__Helix(i)[2], color="b", alpha = 0.1)
            i-=0.1
        #tilepos
        ax.scatter(self.tile_pos[0],self.tile_pos[1],self.tile_pos[2], color="g")
        #hitpos
        ax.scatter(self.__Helix(self.test1)[0], self.__Helix(self.test1)[1], self.__Helix(self.test1)[2],color="y")

        #ax.scatter(self.__Helix(self.test1+self.phi)[0], self.__Helix(self.test1+self.phi)[1], self.__Helix(self.test1+self.phi)[2],color="y")

        plt.show()
"""
