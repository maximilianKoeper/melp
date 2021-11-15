from melp import Detector
import ROOT

def hittimes_in_file (filename):
    hittimes = {}
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
            hittime = ttree_mu3e.tilehit_timestamp[hit_tile_index]
            hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

            if hit_tile in hittimes.keys():
                hittimes[hit_tile].append(hittime)
            else:
                hittimes[hit_tile] = [hittime]

    return hittimes


def hittimes_in_frame (filename, frame_id):
    hittimes = {}
    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")
    
    ttree_mu3e.GetEntry(frame_id)

    for hit_tile_index in range(len(ttree_mu3e.tilehit_tile)):
        hittime = ttree_mu3e.tilehit_timestamp[hit_tile_index]
        hit_tile = ttree_mu3e.tilehit_tile[hit_tile_index]

        if hit_tile in hittimes.keys():
            hittimes[hit_tile].append(hittime)
        else:
            hittimes[hit_tile] = [hittime]

    return hittimes



