# ---------------------------------------------------------------------
#  TRAJECTORY CLASS
# ---------------------------------------------------------------------

import dataclasses


@dataclasses.dataclass
class Trajectory:
    id: int
    v_pos: list
    v_dir: list
    traj_type: int
