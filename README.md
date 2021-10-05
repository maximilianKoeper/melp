# Mu3eHelper for tile analysis (for python) aka melp  

# MELP V2

# IN DEVELOPMENT: UNSTABLE API

## Example:
```
import matplotlib.pyplot as plt
from melp import Detector
import melp
```

```
mu3e_detector = Detector.initFromROOT("run.root")

melp.select(mu3e_detector)
melp.info()

melp.addTileHits("run.root", truth=True, traj=True)

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
