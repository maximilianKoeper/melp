# ---------------------------------------------------------------------
#  HIT CLASS
# ---------------------------------------------------------------------

import dataclasses


@dataclasses.dataclass
class Hit:
    edep: float = 0
    mc_i: int = 0
    tid: int = -1
    frame_id: int = -1
    run_id: int = -1
    hid: int = 0
    impact_vec: list = None
    trajectory: object = None
    pos: list = None
