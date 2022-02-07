import ROOT

import numpy as np

from melp.src.hit import Hit
from melp.src.trajectory import Trajectory

from melp.libs import mathfunctions as mf
from melp.libs import helices as hl
from melp.libs.misc import *
import melp

# ---------------------------------------------------------------------
#  Define global variables and functions to select Detector
# ---------------------------------------------------------------------

__detector__: melp.Detector = None


def select(selection: melp.Detector):
    global __detector__
    __detector__ = selection


def info():
    global __detector__
    print(__detector__)


# ---------------------------------------------------------------------
#  addTileHit information
# ---------------------------------------------------------------------
def addTileHits(filename, traj=True, truth=False):
    global __detector__
    if __detector__ is None:
        print("Error: Detector not selected!")
        return

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e_mc = file.Get("mu3e_mchits")
    ttree_mu3e.GetEntry(0)
    if ttree_mu3e.run in __detector__.TileDetector.AddedRuns:
        print("ERROR: RUN ", ttree_mu3e.run, " already loaded into Detector")
        return
    else:
        __detector__.TileDetector.AddedRuns.append(ttree_mu3e.run)
        pass

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        for i in range(len(ttree_mu3e.tilehit_tile)):
            kwargs = {}
            tile = ttree_mu3e.tilehit_tile[i]
            edep = ttree_mu3e.tilehit_edep[i]
            mc_i = ttree_mu3e.tilehit_mc_i[i]
            ttree_mu3e_mc.GetEntry(mc_i)
            hid = ttree_mu3e_mc.hid
            frame_id = frame

            kwargs["edep"] = edep
            kwargs["mc_i"] = mc_i
            kwargs["hid"] = hid
            kwargs["frame_id"] = frame_id

            # ---------------------
            # Add Truth information
            # ---------------------
            if truth:
                px = ttree_mu3e_mc.p_in_x
                py = ttree_mu3e_mc.p_in_y
                pz = ttree_mu3e_mc.p_in_z
                p_xyz = [px, py, pz]
                kwargs["impact_vec"] = p_xyz

            # ---------------------
            # Add Traj information
            # ---------------------
            if traj:
                tid = ttree_mu3e_mc.tid
                traj_ids = list(ttree_mu3e.traj_ID)
                index = index_finder(traj_ids, tid)
                try:
                    index = int(*index)

                    vx = ttree_mu3e.traj_vx[index]
                    vy = ttree_mu3e.traj_vy[index]
                    vz = ttree_mu3e.traj_vz[index]
                    px = ttree_mu3e.traj_px[index]
                    py = ttree_mu3e.traj_py[index]
                    pz = ttree_mu3e.traj_pz[index]
                    traj_type = ttree_mu3e.traj_type[index]

                    v_xyz = [vx, vy, vz]
                    p_xyz = [px, py, pz]

                    trajectory = Trajectory(id=tid, v_pos=v_xyz, v_dir=p_xyz, traj_type=traj_type)

                    kwargs["trajectory"] = trajectory
                except:
                    pass

            tilehit = Hit(**kwargs)
            __detector__.TileDetector.addHit(tile, tilehit)


# ---------------------------------------------------------------------
#  addSensorHit information
# ---------------------------------------------------------------------
def addSensorHits(filename, traj=True):
    global __detector__
    if __detector__ is None:
        print("Error: Detector not selected!")
        return

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    ttree_mu3e_mc = file.Get("mu3e_mchits")
    ttree_mu3e.GetEntry(0)
    if ttree_mu3e.run in __detector__.SensorsModules.AddedRuns:
        print("ERROR: RUN ", ttree_mu3e.run, " already loaded into Detector")
        return
    else:
        __detector__.SensorsModules.AddedRuns.append(ttree_mu3e.run)
        pass

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        for i in range(ttree_mu3e.Nhit):
            kwargs = {}
            frame_id = frame

            sensorid = ttree_mu3e.hit_pixelid[i]
            mc_i = ttree_mu3e.hit_mc_i[i]
            ttree_mu3e_mc.GetEntry(mc_i)

            pos = __detector__.SensorsModules.getPixelPos(sensorid)

            kwargs["pos"] = pos
            kwargs["frame_id"] = frame_id
            # ---------------------
            # Add Traj information
            # ---------------------

            if traj:
                tid = ttree_mu3e_mc.tid
                kwargs["tid"] = tid

            sensorHit = Hit(**kwargs)
            __detector__.SensorsModules.addHit(sensorid, sensorHit)


# ---------------------------------------------------------------------
# getHitRate
#
# returns: (z_pos, total_hits, primary_hits, secondary_hits, tertiary_hits, edep)
#
# ---------------------------------------------------------------------
def getHitRate(tileID=-1) -> list:
    global __detector__
    if __detector__ is None:
        print("Error: Detector not selected!")
        return

    hitrate = [[0], [0], [0], [0], [0], [0]]
    if tileID == -1:
        for tile in __detector__.TileDetector.tile:
            hitrate[0].append(__detector__.TileDetector.getPos(tile)[2])
            hitrate[1].append(0)
            hitrate[2].append(0)
            hitrate[3].append(0)
            hitrate[4].append(0)
            hitrate[5].append(0)

            for hit in __detector__.TileDetector.tile[tile].hits:

                edep = hit.edep
                hid = hit.hid

                hitrate[1][-1] += 1  # add a hit
                hitrate[5][-1] += edep  # add edep
                if abs(hid) == 1:
                    hitrate[2][-1] += 1  # add a primary hit
                elif abs(hid) == 2:
                    hitrate[3][-1] += 1  # add secondary hit
                elif abs(hid) == 3:
                    hitrate[4][-1] += 1  # add tertiary hit


    else:
        hitrate[0][0] = __detector__.TileDetector.tile[tileID].pos[2]  # z position

        for hit in __detector__.TileDetector.tile[tileID].hits:

            edep = hit.edep
            hid = hit.hid

            hitrate[1][0] += 1  # add a hit
            hitrate[5][0] += edep  # add edep
            if abs(hid) == 1:
                hitrate[2][0] += 1  # add a primary hit
            elif abs(hid) == 2:
                hitrate[3][0] += 1  # add secondary hit
            elif abs(hid) == 3:
                hitrate[4][0] += 1  # add tertiary hit

    __detector__.TileDetector.addRateResult(hitrate)
    return hitrate


# ---------------------------------------------------------------------
#  HIT ANGLE
#
# TODO:
#   - PDG Check
# ---------------------------------------------------------------------
def getHitAngle(tileID=-1, rec_type="Truth", hit_type="primary", angle="phi", particle_type="all") -> list:
    global __detector__
    if __detector__ is None:
        print("Error: Detector not selected!")
        return

    hitangle = [[0.], [0.], rec_type, hit_type, angle, particle_type]

    for tileID in __detector__.TileDetector.tile:
        for hit in __detector__.TileDetector.tile[tileID].hits:

            #############
            # HID CHECK #
            #############
            if hit_type == "primary":
                # only primary hit gets analyzed
                if abs(hit.hid) != 1:
                    continue
            elif hit_type == "secondary":
                if abs(hit.hid) != 2:
                    continue
            elif hit_type == "all":
                pass
            else:
                raise ValueError("hit_type: not supported")

            ##############
            # CALC ANGLE #
            ##############
            if rec_type == "Truth":
                hitangle[1].append(__truth_angle__(angle, hit, tileID))
                hitangle[0].append(__detector__.TileDetector.tile[tileID].pos[2])
            elif rec_type == "Helix":
                tmp_angle = __helix_angle__(angle, hit, tileID)
                if tmp_angle is not None:
                    hitangle[1].append(tmp_angle)
                    hitangle[0].append(__detector__.TileDetector.tile[tileID].pos[2])
            elif rec_type == "TilePixel":
                tmp_angle = __tile_pixel_rec__(angle, hit, tileID)
                if tmp_angle is not None:
                    hitangle[1].append(tmp_angle)
                    hitangle[0].append(__detector__.TileDetector.tile[tileID].pos[2])
            else:
                raise ValueError("hit_type: not supported")

    __detector__.TileDetector.addAngleResult(hitangle)
    return hitangle


# ---------------------------------------------------------------------
def __truth_angle__(angle: str, hit, tileID: int) -> float:
    result: float

    if angle == "norm":
        result = mf.angle_between(hit.impact_vec, __detector__.TileDetector.tile[tileID].dir)
    elif angle == "theta":
        result = mf.angle_between(hit.impact_vec, np.array([0, 0, -1]))
    elif angle == "phi":
        vector = -np.array(__detector__.TileDetector.tile[tileID].dir)
        result = -mf.angle_between_phi(hit.impact_vec[0:2], vector[0:2])
    else:
        raise ValueError('ERROR: angle != [norm, theta, phi]')

    return result


# ---------------------------------------------------------------------
def __helix_angle__(angle: str, hit, tileID: int) -> float:
    result: float = 0
    if hit.trajectory is not None:
        type_1 = abs(int(repr(hit.trajectory.traj_type)[-1]))
        if type_1 == 1 or type_1 == 2:
            v_xyz = hit.trajectory.v_pos
            p_xyz = hit.trajectory.v_dir
            helix = hl.Helices(v_xyz[0], v_xyz[1], v_xyz[2], p_xyz[0], p_xyz[1], p_xyz[2], type_1,
                               __detector__.TileDetector.tile[tileID].pos)
            result = helix.hitAngle(__detector__.TileDetector.tile[tileID].dir, angle)
            del helix
        else:
            return None
    return result


# ---------------------------------------------------------------------
def __tile_pixel_rec__(angle: str, hit, tileID: int) -> float:
    result: float = 0

    return result
