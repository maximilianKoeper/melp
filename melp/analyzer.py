import ROOT

import numpy as np

from melp.src.hit import Hit

from melp.libs import mathfunctions as mf

#---------------------------------------------------------------------
#  HIT RATE
#---------------------------------------------------------------------

class TileAnalyzer ():
    def __init__ (self, detecor, *args):
        self.Detecor = detecor


    #---------------------------------------------------------------------
    #  addTileHits
    #---------------------------------------------------------------------
    def addTileHits(self, filename):
        file          = ROOT.TFile(filename)
        ttree_mu3e    = file.Get("mu3e")
        ttree_mu3e_mc = file.Get("mu3e_mchits")
        ttree_mu3e.GetEntry(0)
        if ttree_mu3e.run in self.Detecor.AddedRuns:
            print("ERROR: RUN ", ttree_mu3e.run, " already loaded into Detecor")
            return
        else:
            self.Detecor.AddedRuns.append(ttree_mu3e.run)
            pass

        for frame in range(ttree_mu3e.GetEntries()):
            ttree_mu3e.GetEntry(frame)
            for i in range(len(ttree_mu3e.tilehit_tile)):
                tile = ttree_mu3e.tilehit_tile[i]
                edep = ttree_mu3e.tilehit_edep[i]
                mc_i = ttree_mu3e.tilehit_mc_i[i]
                ttree_mu3e_mc.GetEntry(mc_i)
                hid  = ttree_mu3e_mc.hid

                # p_in_xyz is not always in Root file
                #try:
                px   = ttree_mu3e_mc.p_in_x
                py   = ttree_mu3e_mc.p_in_y
                pz   = ttree_mu3e_mc.p_in_z

                p_xyz    = [0,0,0]
                p_xyz[0] = px
                p_xyz[1] = py
                p_xyz[2] = pz

                tilehit = Hit(edep = edep, mc_i = mc_i, hid = hid, impact_vec = p_xyz)
                #except:
                #    print("error")
                #    tilehit = Hit(edep = edep, mc_i = mc_i, hid = hid)

                #tilehit = Hit(edep = edep, mc_i = mc_i, hid = hid, impact_vec = p_xyz)
                self.Detecor.TileDetector.addHit(tile, tilehit)

    #---------------------------------------------------------------------
    #  addTileHits
    #
    # returns: (z_pos, total_hits, primary_hits, secondary_hits, tertiary_hits, edep)
    #
    #---------------------------------------------------------------------
    def getHitRate(self, tileID = -1):
        hitrate = [[0],[0],[0],[0],[0],[0]]
        if tileID == -1:
            for tile in self.Detecor.TileDetector.tile:
                hitrate[0].append(self.Detecor.TileDetector.getPos(tile)[2])
                hitrate[1].append(0)
                hitrate[2].append(0)
                hitrate[3].append(0)
                hitrate[4].append(0)
                hitrate[5].append(0)

                for hit in self.Detecor.TileDetector.tile[tile].hits:

                    edep = hit.edep
                    hid  = hit.hid

                    hitrate[1][-1] += 1 # add a hit
                    hitrate[5][-1] += edep # add edep
                    if abs(hid) == 1:
                        hitrate[2][-1] += 1 # add a primary hit
                    elif abs(hid) == 2:
                        hitrate[3][-1] += 1 # add secondary hit
                    elif abs(hid) == 3:
                        hitrate[4][-1] += 1 # add tertiary hit


        else:
            hitrate[0][0] = self.Detecor.TileDetector.tile[tileID].pos[2] # z position

            for hit in self.Detecor.TileDetector.tile[tileID].hits:

                edep = hit.edep
                hid  = hit.hid

                hitrate[1][0] += 1 # add a hit
                hitrate[5][0] += edep # add edep
                if abs(hid) == 1:
                    hitrate[2][0] += 1 # add a primary hit
                elif abs(hid) == 2:
                    hitrate[3][0] += 1 # add secondary hit
                elif abs(hid) == 3:
                    hitrate[4][0] += 1 # add tertiary hit

        self.Detecor.TileDetector.addRateResult(hitrate)
        return hitrate

    #---------------------------------------------------------------------
    #  HIT ANGLE
    #
    # TODO:
    #   - PDG Check
    #---------------------------------------------------------------------
    def getHitAngle(self, tileID = -1, rec_type="Truth", hit_type="primary", angle="phi"):

        hitangle = [[0],[0],rec_type,hit_type,angle]

        for tileID in self.Detecor.TileDetector.tile:
            for hit in self.Detecor.TileDetector.tile[tileID].hits:

                #############
                # HID CHECK #
                #############
                if hit_type == "primary":
                    # only primary hit gets analyzed
                    if hit.hid != 1:
                        continue
                elif hit_type == "secondary":
                    if hit.hid != 1:
                        continue
                elif hit_type == "all":
                    pass
                else:
                    raise ValueError("hit_type: not supported")


                ##############
                # CALC ANGLE #
                ##############
                if  angle == "norm":
                    hitangle[1].append(mf.angle_between(hit.impact_vec, self.Detecor.TileDetector.tile[tileID].dir))
                elif angle == "theta":
                    hitangle[1].append(mf.angle_between(hit.impact_vec, np.array([0,0,-1])))
                elif angle == "phi":
                    vector = -np.array(self.Detecor.TileDetector.tile[tileID].dir)
                    hitangle[1].append(-mf.angle_between_phi(hit.impact_vec[0:2], vector[0:2]))
                else:
                    raise ValueError('ERROR: angle != [norm, theta, phi]')

                hitangle[0].append(self.Detecor.TileDetector.tile[tileID].pos[2])

        self.Detecor.TileDetector.addAngleResult(hitangle)
        return hitangle
