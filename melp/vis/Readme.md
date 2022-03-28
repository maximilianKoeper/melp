# MELP.VIS - Visualizer for single Tracks
___

## Example:
```
from melp.vis import Visualizer
import matplotlib
```

```
vis = Visualizer("run42_20000_sorted_test.root")
vis.set_frame_id(1199)
```

```
vis.list_traj_info()
```

```
vis.set_trajectories([36,39, 60])
```

```
vis.show()
vis.show_3d()
```
