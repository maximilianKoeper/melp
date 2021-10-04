# Mu3eHelper for tile analysis (for python) aka melp  

# MELP V2

# IN DEVELOPMENT

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

melp.addTileHits("run.root")

print(melp.getHitRate())
```
