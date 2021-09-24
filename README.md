# Mu3eHelper for tile analysis (for python) aka melp  

# MELP V2

## Example:
```
import matplotlib.pyplot as plt
from melp import Detector

mu3e_detector = Detector.initFromROOT("run.root")

mu3e_detector.addTileHits("run.root")

z_arr, hit_arr = mu3e_detector.Tiles.rateZ()
plt.plot(z_arr,hit_arr)
plt.show()

mu3e_detector.calcImpactVec("sorted_truth.root")

mu3e_detector.Tiles.calcAngleTruth_byZ()

binned, xedges, yedges = mu3e_detector.Tiles.getBinned()

import numpy as np
import matplotlib as mpl

fig = plt.figure(figsize=(10, 8))

ax = fig.add_subplot(111, title='theta distribution')
X, Y = np.meshgrid(xedges, yedges)
im = ax.pcolormesh(X, Y, binned.T, cmap="turbo", norm = mpl.colors.LogNorm())
plt.ylabel("theta")
plt.xlabel("z")
plt.colorbar(im)
plt.show()
```
