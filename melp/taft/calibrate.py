import ROOT

import numpy as np
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

def calibrate(filename: str, **kwargs):
    global __detector__
    if __detector__ is None:
        print("ERROR: Detector not selected")
        return

    file = ROOT.TFile(filename)
    ttree_mu3e = file.Get("mu3e")


    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        # Analyzing frame
        for hit_tile in ttree_mu3e.tilehit_tile:

            # Look for clusters in z-dir
            if (hit_tile + 56) in ttree_mu3e.tilehit_tile:
                print("Found cluster for z")

            # Look for clusters in phi-dir
            if (hit_tile + 1) in ttree_mu3e.tilehit_tile:
                print("Found cluster for phi")