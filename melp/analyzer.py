import ROOT

from melp.src.hit import Hit

class HitRate ():
    def __init__ (self, detecor, *args):
        self.Detecor = detecor

    def addTileHits(self, filename):
        file          = ROOT.TFile(filename)
        ttree_mu3e    = file.Get("mu3e")
        ttree_mu3e_mc = file.Get("mu3e_mchits")
        ttree_mu3e.GetEntry(0)
        if ttree_mu3e.run in self.Detecor.AddedRuns:
            print("ERROR: RUN already loaded into Detecor")
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

    def getHitRate(self):
        for tile in self.Detecor.TileDetector.tile:
            pass


class HitAngle ():
    def __init__ (self, detecor, *args):
        pass
