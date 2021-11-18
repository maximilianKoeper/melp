[![CodeQL](https://github.com/maximilianKoeper/melp/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/maximilianKoeper/melp/actions/workflows/codeql-analysis.yml)

# Mu3eHelper for tile analysis (for python) aka melp  
___
# MELP V2
- TAFT - Tile alignment for tiles
___
# IN DEVELOPMENT: UNSTABLE API
___
## Improvements to V1:

- detector class can handle everything concerning the geometry of tiles and pixels.
- the detector class can be saved and loaded with all imported run data.
- helix reconstruction is a lot faster and need less memory.
- modular structure makes changes easy
- usage of dataclasses makes it easy to debug functions

## Example:
```
import matplotlib.pyplot as plt
from melp import Detector
import melp
import melp.taft as taft
```

```
mu3e_detector = Detector.initFromROOT("run.root")
```
```
taft.selct(mu3e_detector)

taft.calibrate("time_misal.root")
# Not ready to use!
```
```
melp.select(mu3e_detector)
melp.info()

melp.addTileHits("run.root", truth=True, traj=True)
melp.addSensorHits("sorted_truth.root", traj=True)

mu3e_detector.info()
mu3e_detector.save("detector_savefile")

hitangle = melp.getHitAngle(rec_type="Helix")
```

```
import numpy as np

binned_data, xedges, yedges = np.histogram2d(hitangle[0], hitangle[1], bins=[220, 180])

fig = plt.figure(figsize=(10, 8))
import matplotlib as mpl
ax = fig.add_subplot(111, title='theta distribution')
X, Y = np.meshgrid(xedges, yedges)
im = ax.pcolormesh(X, Y, binned_data.T, cmap="turbo", norm = mpl.colors.LogNorm())
plt.ylabel("theta")
plt.xlabel("z")
plt.colorbar(im)
plt.show()
```
