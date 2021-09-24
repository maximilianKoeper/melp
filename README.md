# Mu3eHelper for tile analysis (for python) aka melp  

# MELP V2

## Example:
```
import matplotlib.pyplot as plt
from melp.src.detecor import Detector

mu3e_detector = Detector.initFromROOT("run.root")

mu3e_detector.addTileHits("run.root")

z_arr, hit_arr = mu3e_detector.Tiles.rateZ()
plt.plot(z_arr,hit_arr)
plt.show()
```
