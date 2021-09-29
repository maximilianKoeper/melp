import ROOT

from melp.src.hit import Hit

#---------------------------------------------------------------------
#  HIT RATE
#---------------------------------------------------------------------

class HitRate ():
    def __init__ (self, detecor, *args):
        self.Detecor = detecor

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

                tilehit = Hit(edep = edep, mc_i = mc_i, hid = hid)
                self.Detecor.TileDetector.addHit(tile, tilehit)

    # hitrate (z_pos, total_hits, primary_hits, secondary_hits, tertiary_hits, edep)
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

            return hitrate


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

            return hitrate


class HitAngle ():
    def __init__ (self, detecor, *args):
        pass
