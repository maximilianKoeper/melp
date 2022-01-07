import matplotlib.pyplot as plt
import numpy as np


# --------------------------------------------------------
def plot_station_calibration(mu3e_detector, station: int):
    f_size = 15

    if station == 1:
        station_offset = 200000
    elif station == 2:
        station_offset = 300000
    else:
        raise ValueError("station=(1 or 2) expected")

    # ---------------------------------------
    # get data
    # ---------------------------------------
    grid_relative = np.zeros((52, 56))
    grid_truth = np.zeros((52, 56))
    grid_calibrated = np.zeros((52, 56))

    for tile_id in mu3e_detector.TileDetector.tile:
        tile = mu3e_detector.TileDetector.tile[tile_id]
        if tile.id >= 300000 and station == 1:
            continue
        if tile.id < 300000 and station == 2:
            continue
        y = tile.row()
        x = tile.column()
        # print(tile_id, x, y)

        grid_relative[x][y] = float(tile.get_offset() - mu3e_detector.TileDetector.tile[station_offset].dt_truth)
        grid_truth[x][y] = tile.get_truth_offset()
        grid_calibrated[x][y] = tile.get_calibrated_offset() + mu3e_detector.TileDetector.tile[station_offset].dt_truth

    # ---------------------------------------
    # Plotting
    # ---------------------------------------
    fig, ax_arr = plt.subplots(2, 2, figsize=(20, 20))

    min = np.min(grid_truth)
    max = np.max(grid_truth)

    # plot truth offsets
    ax = ax_arr[0][0]
    heatplot_truth = ax.imshow(grid_truth.T, cmap='magma', vmin=min, vmax=max)
    ax.set_title("truth offset", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)
    fig.colorbar(heatplot_truth, ax=ax)

    # plot calibrated offsets
    ax = ax_arr[0][1]
    heatplot_calibrated = ax.imshow(grid_calibrated.T, cmap='magma', vmin=min, vmax=max)
    ax.set_title("calibrated offset", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)
    fig.colorbar(heatplot_calibrated, ax=ax)

    # plot deviation
    ax = ax_arr[1][0]
    heatplot_relative = ax.imshow(grid_relative.T, cmap='magma')
    ax.set_title("relative error to truth", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)
    fig.colorbar(heatplot_relative, ax=ax)

    # plot deviation with same range as offsets
    ax = ax_arr[1][1]
    heatplot_calibrated = ax.imshow(grid_relative.T, cmap='magma', vmin=min, vmax=max)
    ax.set_title("relative error from truth (same scaling as offsets)", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)
    fig.colorbar(heatplot_calibrated, ax=ax)

    plt.show()


# --------------------------------------------------------
def plot_calibration(mu3e_detector):
    f_size = 15
    # ---------------------------------------
    # get data
    # ---------------------------------------
    grid_relative_1 = np.zeros((52, 56))
    grid_relative_2 = np.zeros((52, 56))

    for tile_id in mu3e_detector.TileDetector.tile:
        tile = mu3e_detector.TileDetector.tile[tile_id]
        y = tile.row()
        x = tile.column()

        if tile.id >= 300000:
            grid_relative_2[x][y] = float(tile.get_offset() - mu3e_detector.TileDetector.tile[300000].dt_truth)
        if tile.id < 300000:
            grid_relative_1[x][y] = float(tile.get_offset() - mu3e_detector.TileDetector.tile[200000].dt_truth)

    # ---------------------------------------
    # Plotting
    # ---------------------------------------
    fig, ax_arr = plt.subplots(1, 2, figsize=(20, 10))

    # plot truth offsets
    ax = ax_arr[0]
    heatplot_truth = ax.imshow(grid_relative_1.T, cmap='magma')
    ax.set_title("error station 1", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)
    fig.colorbar(heatplot_truth, ax=ax)

    # plot calibrated offsets
    ax = ax_arr[1]
    heatplot_calibrated = ax.imshow(grid_relative_2.T, cmap='magma')
    ax.set_title("error station 2", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)
    fig.colorbar(heatplot_calibrated, ax=ax)

    plt.show()

    min1 = np.min(grid_relative_1)
    max1 = np.max(grid_relative_1)
    print("Station 1: max error: ", np.round(max1, 5), " min error: ", np.round(min1, 5))

    min2 = np.min(grid_relative_2)
    max2 = np.max(grid_relative_2)
    print("Station 2: max error: ", np.round(max2, 5), " min error: ", np.round(min2, 5))


def plot_error_dist(mu3e_detector):
    f_size = 15
    # ---------------------------------------
    # get data
    # ---------------------------------------
    relative_1 = []
    relative_2 = []

    for tile_id in mu3e_detector.TileDetector.tile:
        tile = mu3e_detector.TileDetector.tile[tile_id]

        if tile.id >= 300000:
            relative_2.append(float(tile.get_offset() - mu3e_detector.TileDetector.tile[300000].dt_truth))
        if tile.id < 300000:
            relative_1.append(float(tile.get_offset() - mu3e_detector.TileDetector.tile[200000].dt_truth))

    # ---------------------------------------
    # Plotting
    # ---------------------------------------
    fig, ax_arr = plt.subplots(1, 2, figsize=(20, 10))

    # plot truth offsets
    ax = ax_arr[0]
    ax.hist(relative_1, bins=50)
    ax.set_title("error station 1", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)

    # plot calibrated offsets
    ax = ax_arr[1]
    ax.hist(relative_2, bins=50)
    ax.set_title("error station 2", fontsize=f_size)
    ax.set_ylabel("phi (column)", fontsize=f_size)
    ax.set_xlabel("z", fontsize=f_size)

    plt.show()

    min1 = np.min(relative_1)
    max1 = np.max(relative_1)
    print("Station 1: max error: ", np.round(max1, 5), " min error: ", np.round(min1, 5))

    min2 = np.min(relative_2)
    max2 = np.max(relative_2)
    print("Station 2: max error: ", np.round(max2, 5), " min error: ", np.round(min2, 5))

def plot_correction_function(cal_function, *popt):
    x1 = np.linspace(0, 51, 52)
    y1 = np.linspace(0, 55, 56)
    x2, y2 = np.meshgrid(x1, y1)

    plt.imshow(cal_function((x2, y2), *popt) - cal_function((0, 0), *popt), cmap="magma")
    plt.colorbar()
    plt.show()