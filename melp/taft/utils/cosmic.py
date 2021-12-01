import ROOT
import numpy as np


def find_cosmic_events(filename: str, **kwargs):
    root_file = ROOT.TFile.Open(filename, "READ")
    ttree_mu3e = root_file.Get(kwargs["ttree_loc"])

    frame_with_hit_counter = 0
    time = []
    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        if len(ttree_mu3e.tilehit_tile) > 0:
            frame_with_hit_counter += 1
            time.append(sum(ttree_mu3e.tilehit_time)/len(ttree_mu3e.tilehit_time) - min(ttree_mu3e.tilehit_time))

    print(frame_with_hit_counter, " von ", ttree_mu3e.GetEntries())
    return np.array(time)
