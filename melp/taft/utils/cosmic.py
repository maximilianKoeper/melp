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
    for frame in range(ttree_mu3e.GetEntries()):
        ttree_mu3e.GetEntry(frame)

        if frame % 100000 == 0:
            print(round(frame/ttree_mu3e.GetEntries() * 100), " % | Total Frames: ", ttree_mu3e.GetEntries(), end='\r')

        if len(ttree_mu3e.tilehit_tile) > 0:
            frame_with_hit_counter += 1
            time_total.append(sum(ttree_mu3e.tilehit_time) / len(ttree_mu3e.tilehit_time) - min(ttree_mu3e.tilehit_time))

            tmp_arr_1 = []
            tmp_arr_2 = []
            for hit_index in range(len(ttree_mu3e.tilehit_tile)):
                tile_id = ttree_mu3e.tilehit_tile[hit_index]
                if tile_id > 300000:
                    tmp_arr_2.append(ttree_mu3e.tilehit_time[hit_index])
                else:
                    tmp_arr_1.append(ttree_mu3e.tilehit_time[hit_index]+3)
            if len(tmp_arr_1) > 0 and len(tmp_arr_2) > 0:
                # global_time = min(ttree_mu3e.tilehit_time)
                tmp_time_1 = sum(tmp_arr_1)/len(tmp_arr_1)
                tmp_time_2 = sum(tmp_arr_2)/len(tmp_arr_2)

                time_dist_betw_stations.append(tmp_time_1-tmp_time_2)

    print("100 % | Total Frames: ", ttree_mu3e.GetEntries())
    print("Frames with hits: ", frame_with_hit_counter)
    return np.array(time_total), np.array(time_dist_betw_stations)
