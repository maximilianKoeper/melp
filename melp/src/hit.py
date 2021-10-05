# ---------------------------------------------------------------------
#  HIT CLASS
# ---------------------------------------------------------------------

import dataclasses


@dataclasses.dataclass
class Hit:
    edep: float = 0
    mc_i: int = 0
    traj_id: int = -1
    run_id: int = -1
    hid: int = 0
    impact_vec: list = dataclasses.field(default_factory=list)
    trajectory: object = None
