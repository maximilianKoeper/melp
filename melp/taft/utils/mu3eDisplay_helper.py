import dataclasses


@dataclasses.dataclass
class Trajectory:
    tile1_pos: list
    tile2_pos: list


def generate_txt_event_file(trajectories: list, nevents: int) -> bool:
    file = open("Display_cosmics.txt", "w")
    file.write("Event\n")
    index = 0
    for traj in trajectories:
        if index == nevents:
            file.write("Event\n")
            index = 0
        pos1 = traj.tile1_pos
        pos2 = traj.tile2_pos
        str_tmp = f'Line {round(pos1[0])} {round(pos1[1])} {round(pos1[2])} {round(pos2[0])} {round(pos2[1])} {round(pos2[2])}\n'
        file.write(str_tmp)
        index += 1

    file.close()
    return True
