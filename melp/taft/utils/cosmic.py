import ROOT
import numpy as np


# --------------------------------------------------
# searches all frames for hits
# returns time dist from these hits
def find_cosmic_events(filename: str, **kwargs):
    root_file = ROOT.TFile.Open(filename, "READ")
    ttree_mu3e = root_file.Get(kwargs["ttree_loc"])

    frame_with_hit_counter = 0
    time_total = []
    time_dist_betw_stations = []
    event_size_per_station = []

    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        if frame % 100000 == 0:
            print(round(frame / ttree_mu3e.GetEntries() * 100), " % | Total Frames: ", ttree_mu3e.GetEntries(),
                  end='\r')

        if ttree_mu3e.Ntilehit > 0:
            frame_with_hit_counter += 1
            time_total.append(
                sum(ttree_mu3e.tilehit_time) / len(ttree_mu3e.tilehit_time) - min(ttree_mu3e.tilehit_time))

            tmp_arr_1 = []
            tmp_arr_2 = []
            for hit_index in range(len(ttree_mu3e.tilehit_tile)):
                tile_id = ttree_mu3e.tilehit_tile[hit_index]
                if tile_id > 300000:
                    tmp_arr_2.append(ttree_mu3e.tilehit_time[hit_index])
                else:
                    tmp_arr_1.append(ttree_mu3e.tilehit_time[hit_index] + 3)
            if len(tmp_arr_1) > 0 and len(tmp_arr_2) > 0:
                event_size_per_station.append(len(tmp_arr_2))
                event_size_per_station.append(len(tmp_arr_1))
                # global_time = min(ttree_mu3e.tilehit_time)
                tmp_time_1 = sum(tmp_arr_1) / len(tmp_arr_1)
                tmp_time_2 = sum(tmp_arr_2) / len(tmp_arr_2)

                time_dist_betw_stations.append(tmp_time_1 - tmp_time_2)

    print("100 % | Total Frames: ", ttree_mu3e.GetEntries())
    print("Frames with hits: ", frame_with_hit_counter)
    return np.array(time_total), np.array(time_dist_betw_stations), np.array(event_size_per_station)


def station_station_timing(filename: str, detector, **kwargs):
    root_file = ROOT.TFile.Open(filename, "READ")
    ttree_mu3e = root_file.Get(kwargs["ttree_loc"])

    time_dist_betw_stations = []

    it = find_next_cosmic_event(ttree_mu3e, 0, 1)
    while it != -1:

        if it % 100000 == 0:
            print(round(it / ttree_mu3e.GetEntries() * 100), " % | Total Frames: ", ttree_mu3e.GetEntries(),
                  end='\r')

        if kwargs["mc_primary"] is False:
            test_dict = check_cosmic_events(ttree_mu3e, None, None)
        else:
            test_dict = check_cosmic_events_mc(ttree_mu3e, None, None)

        for key in test_dict:
            tmp_ids = test_dict[key][0]
            if any(y > 300000 for y in tmp_ids) and any(y < 300000 for y in tmp_ids):
                # print(test_dict[key][0])
                tmp_arr_2 = []
                tmp_arr_1 = []
                for hit_index in range(len(tmp_ids)):
                    tile_id = tmp_ids[hit_index]
                    if tile_id >= 300000:
                        tmp_arr_2.append(test_dict[key][1][hit_index])
                    else:
                        tmp_arr_1.append(test_dict[key][1][hit_index] + 3)

                tmp_time_1 = sum(tmp_arr_1) / len(tmp_arr_1)
                tmp_time_2 = sum(tmp_arr_2) / len(tmp_arr_2)

                tof = 0.
                pos1 = detector.TileDetector.tile[max(tmp_ids)].pos
                pos2 = detector.TileDetector.tile[min(tmp_ids)].pos
                dist = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)  # mm
                dist *= 0.001  # m
                tof = (dist / 299792458) * (10 ** 9)

                if kwargs["tof"]:
                    if tmp_time_2 > tmp_time_1:
                        time_dist_betw_stations.append((tmp_time_1 - tmp_time_2) + tof)
                    else:
                        time_dist_betw_stations.append((tmp_time_1 - tmp_time_2) - tof)
                else:
                    time_dist_betw_stations.append((tmp_time_1 - tmp_time_2))  # - tof)
        it += 1
        it = find_next_cosmic_event(ttree_mu3e, it, 1)

    return time_dist_betw_stations


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
            for hit_index in range(len(ttree_mu3e.tilehit_tile)):
                tile_id = ttree_mu3e.tilehit_tile[hit_index]
                if station_offset <= tile_id < station_offset + 100000:
                    return it

    return -1


# --------------------------------------
# checks event for useful data
def check_cosmic_events(ttree_mu3e, it: int, station: int):
    # sort hits in time
    tilehit_times, tilehit_ids = (list(t) for t in
                                  zip(*sorted(zip(list(ttree_mu3e.tilehit_time), list(ttree_mu3e.tilehit_tile)))))
    # print(tilehit_ids)
    # print(np.array(tilehit_times) - min(tilehit_times))

    sorted_tracks = get_single_tracks_time_cut(tilehit_ids, tilehit_times)

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
# checks event for useful data
def check_cosmic_events_mc(ttree_mu3e, it: int, station: int):
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
