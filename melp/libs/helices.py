import ROOT
import numpy as np
from scipy.optimize import minimize

# DEBUG: Delete after testing is finished
import matplotlib.pyplot as plt

from melp.libs import mathfunctions as mf

class Helices ():
    # WARNING: NOT TESTED
    def __init__ (self, vx, vy, vz, px, py, pz, type, tile_pos):

        self.bfield     = -1

        self.z0         = vz
        self.type       = type

        self.tile_pos   = tile_pos

        pt              = np.hypot(px, py)

        self.r          = pt / (0.3 * self.bfield)
        self.theta      = np.arctan2(pz, pt)
        self.phi        = np.arctan2(py, px)

        if self.type == 2:
            self.phi += np.pi / 2
        elif self.type == 3:
            self.phi -= np.pi / 2
        else:
            raise ValueError('Helices: init: type not supported')


        self.dz         = 2*np.pi * self.r * np.tan(self.theta);
        self.xc         = self.r*np.cos(self.phi) + vx;
        self.yc         = self.r*np.sin(self.phi) + vy;

        self.xy_vec     = 0.

        # DEBUG: FOR TESTING only
        self.test1 = 0


    #####################
    # private functions #
    #####################

    def __Circle (self, alpha):
        xy    = np.zeros(2)
        xy[0] = self.xc - self.r * np.cos(alpha)
        xy[1] = self.yc - self.r * np.sin(alpha)

        return xy

    def __Minimize_Func_Angle (self, alpha):
        xy_circ    = self.__Circle(alpha)
        xy_tile    = np.zeros(2)
        xy_tile[0] = self.tile_pos[0]
        xy_tile[1] = self.tile_pos[1]
        distance   = mf.distance_between_2d(xy_tile, xy_circ)

        return distance

    def __Get_XY_Plane_Phi (self):
        return minimize(self.__Minimize_Func_Angle, 0).x


    def __Get_Primary_Tile_Hit_Vector_XY (self):
        temp_phi    = self.__Get_XY_Plane_Phi()
        v1_tmp      = self.__Circle(temp_phi)

        if self.type == 3:
            offset = - 0.1
        else:
            offset = + 0.1

        v2_tmp        = self.__Circle(temp_phi + offset)
        xy_hit_vector = np.array(v1_tmp) - np.array(v2_tmp)
        return xy_hit_vector


    def __Get_Primary_Tile_Hit_Angle (self, tile_norm_vec , angle="phi"):
        """
            TODO:
                - add theta
                - add total
        """
        norm_vec = np.array(tile_norm_vec[0:2])
        temp_vec = mf.angle_between(self.__Get_Primary_Tile_Hit_Vector_XY(), norm_vec)

        return temp_vec


    #####################
    # public  functions #
    #####################

    def hitAngle(self, tile_norm_vec , angle="phi"):
        return self.__Get_Primary_Tile_Hit_Angle(tile_norm_vec , angle)


    #####################
    # TESTING functions #
    #####################

    def test (self):
        self.test1 = minimize(self.__Minimize_Func_Angle, 0).x

    def plottest (self):

        plt.scatter(self.__Circle(self.phi)[0],self.__Circle(self.phi)[1], color="r")
        i=self.phi
        while i > -5:
            plt.scatter(self.__Circle(i)[0],self.__Circle(i)[1], color="b", alpha = 0.1)
            i-=0.1
        plt.scatter(self.tile_pos[0],self.tile_pos[1])
        plt.scatter(self.__Circle(self.test1)[0],self.__Circle(self.test1)[1])

        plt.show()
