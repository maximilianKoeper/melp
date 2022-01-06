import warnings

import ROOT
import numpy as np

from melp.libs.misc import index_finder
from melp.taft.utils.mu3eDisplay_helper import Trajectory, generate_txt_event_file


# ---------------------------------------------
def cosmic_correction_z(__detector__, **kwargs):
    print("analyzing cosmics file #1")
    kwargs["station"] = 1
    hist_z_1 = cosmic_linear_correction(kwargs["cosmic_file"], __detector__, **kwargs)
    correction_1 = np.median(hist_z_1)
    print("analyzing cosmics file #2")
    kwargs["station"] = 2
    hist_z_2 = cosmic_linear_correction(kwargs["cosmic_file"], __detector__, **kwargs)
    correction_2 = np.median(hist_z_2)
    print("done")
    phi = __detector__.TileDetector.column_ids(0, 200000)

    for p in range(len(phi)):
        corr = []
        row = __detector__.TileDetector.row_ids(p, 200000)
        for i in range(len(row)):
            corr.append(correction_1 * i)
        for i in range(len(row)):
            __detector__.TileDetector.tile[row[i]].dt_cal -= corr[i]

    for p in range(len(phi)):
        corr = []
        row = __detector__.TileDetector.row_ids(p, 300000)
        for i in range(len(row)):
            corr.append(correction_2 * i)
        for i in range(len(row)):
            __detector__.TileDetector.tile[row[i]].dt_cal -= corr[i]


# --------------------------------------------------
#
def cosmic_linear_correction(filename: str, detector, **kwargs):
    root_file = ROOT.TFile.Open(filename, "READ")
    ttree_mu3e = root_file.Get(kwargs["ttree_loc"])
    time_dist_z = []

    trajectories = []

    # it -> iterator (frame_id). -1 if EOF
    it = find_next_cosmic_event(ttree_mu3e, it=0, station=kwargs["station"])
    while it != -1:

        if it % 100000 == 0:
            print(round(it / ttree_mu3e.GetEntries() * 100), " % | Total Frames: ", ttree_mu3e.GetEntries(),
                  end='\r')

        if kwargs["mc_primary"] is False:
            test_dict = check_cosmic_events(ttree_mu3e)
        else:
            test_dict = check_cosmic_events_mc(ttree_mu3e)

        for key in test_dict:
            tmp_ids = test_dict[key][0]

            if kwargs["station"] == 1:
                if any(y >= 300000 for y in tmp_ids):
                    continue
            elif kwargs["station"] == 2:
                if any(y < 300000 for y in tmp_ids):
                    continue

            tmp_time_1 = min(test_dict[key][1])
            tmp_tile_id_1 = test_dict[key][0][int(*index_finder(list(test_dict[key][1]), tmp_time_1))]
            tmp_time_1 += detector.TileDetector.tile[tmp_tile_id_1].get_offset()

            tmp_time_2 = max(test_dict[key][1])
            tmp_tile_id_2 = test_dict[key][0][int(*index_finder(list(test_dict[key][1]), tmp_time_2))]
            tmp_time_2 += detector.TileDetector.tile[tmp_tile_id_2].get_offset()
            # calculating tof
            pos1 = detector.TileDetector.tile[tmp_tile_id_1].pos
            pos2 = detector.TileDetector.tile[tmp_tile_id_2].pos
            dist = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)  # mm
            dist *= 0.001  # m
            tof = (dist / 299792458) * (10 ** 9)  # ns

            z_dist = (detector.TileDetector.tile[tmp_tile_id_1].column() - detector.TileDetector.tile[
                tmp_tile_id_2].column())

            if abs(z_dist) > kwargs["cosmic_threshold"]:
                time_dist_z.append((abs(tmp_time_1 - tmp_time_2) - tof) / z_dist)

            trajectories.append(Trajectory(tile1_pos=pos1, tile2_pos=pos2))

        it += 1
        it = find_next_cosmic_event(ttree_mu3e, it, 1)

    generate_txt_event_file(trajectories, 500)
    return time_dist_z


# --------------------------------------
def station_station_timing(filename: str, detector, **kwargs):
    trajectories = []

    root_file = ROOT.TFile.Open(filename, "READ")
    ttree_mu3e = root_file.Get(kwargs["ttree_loc"])

    time_dist_betw_stations = []

    it = find_next_cosmic_event(ttree_mu3e, it=0, station=1)
    while it != -1:

        if it % 100000 == 0:
            print(round(it / ttree_mu3e.GetEntries() * 100), " % | Total Frames: ", ttree_mu3e.GetEntries(),
                  end='\r')

        if kwargs["mc_primary"] is False:
            test_dict = check_cosmic_events(ttree_mu3e)
        else:
            test_dict = check_cosmic_events_mc(ttree_mu3e)

        for key in test_dict:
            tmp_ids = test_dict[key][0]
            if any(y > 300000 for y in tmp_ids) and any(y < 300000 for y in tmp_ids):
                # print(test_dict[key][0])
                tmp_time_arr_2 = []
                tmp_id_arr_2 = []
                tmp_time_arr_1 = []
                tmp_id_arr_1 = []
                for hit_index in range(len(tmp_ids)):
                    tile_id = tmp_ids[hit_index]
                    if tile_id >= 300000:
                        tmp_time_arr_2.append(test_dict[key][1][hit_index])
                        tmp_id_arr_2.append(tile_id)
                    else:
                        tmp_time_arr_1.append(test_dict[key][1][hit_index])
                        tmp_id_arr_1.append(tile_id)

                # tmp_time_1 = sum(tmp_time_arr_1) / len(tmp_time_arr_1)
                # tmp_time_2 = sum(tmp_time_arr_2) / len(tmp_time_arr_2)

                # first and last hit used for timing information
                tilehit_times_2, tilehit_ids_2 = (list(t) for t in zip(*sorted(zip(tmp_time_arr_2, tmp_id_arr_2))))
                tmp_time_2 = tilehit_times_2[-1]
                tilehit_times_1, tilehit_ids_1 = (list(t) for t in zip(*sorted(zip(tmp_time_arr_1, tmp_id_arr_1))))
                tmp_time_1 = tilehit_times_1[0]

                # tof = 0.
                pos1 = detector.TileDetector.tile[tilehit_ids_1[0]].pos
                pos2 = detector.TileDetector.tile[tilehit_ids_2[-1]].pos
                # pos1 = detector.TileDetector.tile[max(tmp_ids)].pos
                # pos2 = detector.TileDetector.tile[min(tmp_ids)].pos
                dist = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)  # mm
                dist *= 0.001  # m
                tof = (dist / 299792458) * (10 ** 9)

                if kwargs["tof"]:
                    if tmp_time_2 > tmp_time_1:
                        time_dist_betw_stations.append((tmp_time_1 - tmp_time_2) + tof)
                    else:
                        time_dist_betw_stations.append((tmp_time_1 - tmp_time_2) - tof)
                else:
                    time_dist_betw_stations.append((tmp_time_1 - tmp_time_2))

                trajectories.append(Trajectory(tile1_pos=pos1, tile2_pos=pos2))

        it += 1
        it = find_next_cosmic_event(ttree_mu3e, it, 1)

    generate_txt_event_file(trajectories, 500)
    return time_dist_betw_stations


# --------------------------------------------------
#
def get_cosmic_data_from_file(filename: str, detector, **kwargs):
    root_file = ROOT.TFile.Open(filename, "READ")
    ttree_mu3e = root_file.Get(kwargs["ttree_loc"])
    time_offset_between_hits = []
    position_hits_hit1_column = []
    position_hits_hit1_row = []
    position_hits_hit2_column = []
    position_hits_hit2_row = []

    # it -> iterator (frame_id). -1 if EOF
    it = find_next_cosmic_event(ttree_mu3e, it=0, station=kwargs["station"])
    while it != -1:
        # TODO: just for debugging
        if it >= 1500000:
            it = 100000000000001
        if it % 100000 == 0:
            print(round(it / ttree_mu3e.GetEntries() * 100), " % | Total Frames: ", ttree_mu3e.GetEntries(),
                  end='\r')

        if kwargs["mc_primary"] is False:
            test_dict = check_cosmic_events(ttree_mu3e)
        else:
            test_dict = check_cosmic_events_mc(ttree_mu3e)

        for key in test_dict:
            tmp_ids = test_dict[key][0]

            if kwargs["station"] == 1:
                if any(y >= 300000 for y in tmp_ids):
                    continue
            elif kwargs["station"] == 2:
                if any(y < 300000 for y in tmp_ids):
                    continue

            tmp_time_1 = min(test_dict[key][1])
            tmp_tile_id_1 = test_dict[key][0][int(*index_finder(list(test_dict[key][1]), tmp_time_1))]
            tmp_time_1 += detector.TileDetector.tile[tmp_tile_id_1].get_offset()

            tmp_time_2 = max(test_dict[key][1])
            tmp_tile_id_2 = test_dict[key][0][int(*index_finder(list(test_dict[key][1]), tmp_time_2))]
            tmp_time_2 += detector.TileDetector.tile[tmp_tile_id_2].get_offset()

            # calculating tof
            pos1 = detector.TileDetector.tile[tmp_tile_id_1].pos
            pos2 = detector.TileDetector.tile[tmp_tile_id_2].pos
            dist = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)  # mm
            dist *= 0.001  # m
            tof = (dist / 299792458) * (10 ** 9)  # ns

            if abs(dist) >= 0.05:
                time_offset_between_hits.append((abs(tmp_time_1 - tmp_time_2) - tof))
                position_hits_hit1_column.append(detector.TileDetector.tile[tmp_tile_id_1].column())
                position_hits_hit1_row.append(detector.TileDetector.tile[tmp_tile_id_1].row())
                position_hits_hit2_column.append(detector.TileDetector.tile[tmp_tile_id_2].column())
                position_hits_hit2_row.append(detector.TileDetector.tile[tmp_tile_id_2].row())

        it += 1
        it = find_next_cosmic_event(ttree_mu3e, it, 1)

    position_hits = (position_hits_hit1_column, position_hits_hit1_row, position_hits_hit2_column, position_hits_hit2_row)

    return time_offset_between_hits, position_hits


# --------------------------------------
# returns index of next frame with cosmic event
# it: iterator -> frame_id
def find_next_cosmic_event(ttree_mu3e, it: int, station, threshold=1) -> int:
    # EOF
    if it >= ttree_mu3e.GetEntries():
        return -1

    station_offset = 200000
    if station == 2:
        station_offset += 100000

    # searches next frame with at least the number of frames defined in threshold
    # (default 1) in the station defined by station
    start = it
    for it in range(start, ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(it)
        if ttree_mu3e.Ntilehit >= threshold:
            if station == 3:
                return it
            for hit_index in range(len(ttree_mu3e.tilehit_tile)):
                tile_id = ttree_mu3e.tilehit_tile[hit_index]
                if station_offset <= tile_id < station_offset + 100000:
                    return it

    return -1


# --------------------------------------
# checks event for useful data
def check_cosmic_events(ttree_mu3e):
    # sort hits in time
    tilehit_times, tilehit_ids = (list(t) for t in
                                  zip(*sorted(zip(list(ttree_mu3e.tilehit_time), list(ttree_mu3e.tilehit_tile)))))
    # print(tilehit_ids)
    # print(np.array(tilehit_times) - min(tilehit_times))

    sorted_tracks = get_single_tracks_time_cut(tilehit_ids, tilehit_times)

    return sorted_tracks


# --------------------------------------
# checks event for useful data
# with primary id information
def check_cosmic_events_mc(ttree_mu3e):
    # sort hits for primary id
    indices = np.argsort(list(ttree_mu3e.tilehit_primary))
    tilehit_times = np.asarray(list(ttree_mu3e.tilehit_time))[indices]
    tilehit_ids = np.asarray(list(ttree_mu3e.tilehit_tile))[indices]
    tilehit_primaries = np.asarray(list(ttree_mu3e.tilehit_primary))[indices]

    # print(tilehit_ids)
    # print(np.array(tilehit_times) - min(tilehit_times))

    sorted_tracks = get_single_tracks_primary_mc(tilehit_ids, tilehit_times, tilehit_primaries)

    return sorted_tracks


# --------------------------------------
# split list into chunks for each event (time cut)
# returns dictionary with data (key is arbitrary but unique)
def get_single_tracks_time_cut(tilehit_ids: list, tilehit_times: list, threshold: float = 8.) -> dict:
    single_events = {}

    tmp_time_reference = tilehit_times[0]
    index_start_track = 0
    for index in range(len(tilehit_times)):
        if abs(tmp_time_reference - tilehit_times[index]) > threshold:
            single_events[index] = [tilehit_ids[index_start_track:index], tilehit_times[index_start_track:index]]
            index_start_track = index
            tmp_time_reference = tilehit_times[index]

    # fill up remaining event
    if index_start_track != len(tilehit_times):
        single_events[len(tilehit_times)] = [tilehit_ids[index_start_track:],
                                             tilehit_times[index_start_track:]]

    # delete entries with only one hit
    keys_to_delete = []
    for key in single_events:
        if len(single_events[key][0]) == 1:
            keys_to_delete.append(key)

    for key in keys_to_delete:
        del single_events[key]

    return single_events


# --------------------------------------
# split list into chunks for each event (mc event)
# returns dictionary with data (key is arbitrary but unique)
def get_single_tracks_primary_mc(tilehit_ids: list, tilehit_times: list, tilehit_primaries: list) -> dict:
    single_events = {}

    tmp_primary_reference = tilehit_primaries[0]
    index_start_track = 0
    for index in range(len(tilehit_times)):
        if tilehit_primaries[index] != tmp_primary_reference:
            single_events[index] = [tilehit_ids[index_start_track:index], tilehit_times[index_start_track:index]]
            index_start_track = index
            tmp_primary_reference = tilehit_primaries[index]

    # fill up remaining event
    if index_start_track != len(tilehit_times):
        single_events[len(tilehit_times)] = [tilehit_ids[index_start_track:],
                                             tilehit_times[index_start_track:]]

    # delete entries with only one hit
    keys_to_delete = []
    for key in single_events:
        if len(single_events[key][0]) == 1:
            keys_to_delete.append(key)

    for key in keys_to_delete:
        del single_events[key]

    return single_events
