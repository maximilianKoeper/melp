# MELP.VIS - Visualizer for single Tracks
___

## Visualizer
- import Visualizer and import matplotlib
```
from melp.vis import Visualizer
import matplotlib
```

- initialize Visualizer (sorted root file from mu3e (v.4.6 tested))
```
vis = Visualizer("run42_20000_sorted_test.root")
```

### Functions:
- set_frame_id(int)
- set_trajectories(list)
- add_toy_event(Trajectory (returned by ToyEventGenerator.new_event()))
- select_all_trajectories()
- select_all_toy_trajectories()
- reset_frame()
- show()
- show_3d()
- list_traj_info()

___
## ToyEventGenerator
- import ToyEventVisualizer
```
from melp.vis import ToyEventGenerator
```

- initialize ToyEventGenerator
```
TEG = ToyEventGenerator(pt=20, pz=-15, particle_type=1)
```
(electron = 2 / positron = 1)
### Functions:
- new_event() -> Trajectory object

#### example:
```
vis.add_toy_event(TEG.new_event())
vis.select_all_toy_trajectories()
vis.show()
```