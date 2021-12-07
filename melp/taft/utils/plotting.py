import matplotlib.pyplot as plt
import numpy as np


# --------------------------------------------------------
def plot_station_calibration(mu3e_detector, station: int):
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

        grid_relative[x][y] = float(
            tile.dt_truth - tile.dt_cal - mu3e_detector.TileDetector.tile[station_offset].dt_truth)
        grid_truth[x][y] = tile.dt_truth
        grid_calibrated[x][y] = tile.dt_cal + mu3e_detector.TileDetector.tile[station_offset].dt_truth

    #grid_relative[0][0] = 0
    #grid_truth[0][0] = mu3e_detector.TileDetector.tile[station_offset].dt_truth

    # ---------------------------------------
    # Plotting
    # ---------------------------------------
    fig, ax_arr = plt.subplots(2, 2, figsize=(20, 20))

    min = np.min(grid_truth)
    max = np.max(grid_truth)

    # plot truth offsets
    ax = ax_arr[0][0]
    heatplot_truth = ax.imshow(grid_truth.T, cmap='magma', vmin=min, vmax=max)
    ax.set_title("truth offset")
    fig.colorbar(heatplot_truth, ax=ax)

    # plot calibrated offsets
    ax = ax_arr[0][1]
    heatplot_calibrated = ax.imshow(grid_calibrated.T, cmap='magma', vmin=min, vmax=max)
    ax.set_title("calibrated offset")
    fig.colorbar(heatplot_calibrated, ax=ax)

    # plot deviation
    ax = ax_arr[1][0]
    heatplot_relative = ax.imshow(grid_relative.T, cmap='magma')
    ax.set_title("relative error to truth")
    fig.colorbar(heatplot_relative, ax=ax)

    # plot deviation with same range as offsets
    ax = ax_arr[1][1]
    heatplot_calibrated = ax.imshow(grid_relative.T, cmap='magma', vmin=min, vmax=max)
    ax.set_title("relative error from truth (same scaling as offsets)")
    fig.colorbar(heatplot_calibrated, ax=ax)

    plt.show()


# --------------------------------------------------------
def plot_calibration(mu3e_detector):
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
            grid_relative_2[x][y] = float(
                tile.dt_truth - tile.dt_cal - mu3e_detector.TileDetector.tile[300000].dt_truth)
        if tile.id < 300000:
            grid_relative_1[x][y] = float(
                tile.dt_truth - tile.dt_cal - mu3e_detector.TileDetector.tile[200000].dt_truth)

    # ---------------------------------------
    # Plotting
    # ---------------------------------------
    fig, ax_arr = plt.subplots(1, 2, figsize=(20, 10))

    # plot truth offsets
    ax = ax_arr[0]
    heatplot_truth = ax.imshow(grid_relative_1.T, cmap='magma')
    ax.set_title("truth offset")
    fig.colorbar(heatplot_truth, ax=ax)

    # plot calibrated offsets
    ax = ax_arr[1]
    heatplot_calibrated = ax.imshow(grid_relative_2.T, cmap='magma')
    ax.set_title("calibrated offset")
    fig.colorbar(heatplot_calibrated, ax=ax)

    plt.show()
