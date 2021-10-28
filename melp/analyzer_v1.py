import ROOT

import numpy as np

from melp.libs import mathfunctions as mf
import melp.analyzer_v2 as v2


def getHitAngleRec(filename, rec_type="pixel", hit_type="primary", angle="phi", particle_type="all") -> list:
    if v2.__detector__ is None:
        print("Error: Detector not selected!")
        return

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    #ttree_mu3e.GetEntry(0)

    hitangle = [[0.], [0.], rec_type, hit_type, angle, particle_type]

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)
        for tileid in range(ttree_mu3e.tilehit_tile):
            pass