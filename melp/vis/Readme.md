# MELP.VIS - Visualizer for single Tracks
___
<img src="../../pictures/vis_2d.png" alt="vis_2d" style="width:300px;"/>
<img src="../../pictures/vis_3d.png" alt="vis_3d" style="width:300px;"/>

### Example Jupyter Notebook:
[Jupyter Notebook](../../jupyter-notbooks/melpVis/melpVis.ipynb)

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
- select_only_particle_types(list (e.g. \[1,2\]))
- reset_frame()
- set_plt_options_2d(**kwargs)  -> xlim, ylim, zlim, dpi, fontsize 
- show()
- show_3d()
- list_traj_info()

___

- 0: photon
- 1: e+
- 2: e+
- 3: mu+
- 4: mu-

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
### Functions:
- new_event() -> Trajectory object

#### example:
```
vis.add_toy_event(TEG.new_event())
vis.select_all_toy_trajectories()
vis.show()
```