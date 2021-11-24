# ---------------------------------------
def loop_correction_phi(detector, dt_phi_rel: dict, station: int):
    for z in range(len(detector.TileDetector.row_ids(0, station))):
        tile_ids = detector.TileDetector.column_ids(z, station)

        sum_dt_column = 0.
        number_unfilled_dt = 0
        for id_index in tile_ids:
            try:
                sum_dt_column += dt_phi_rel[id_index]
            except KeyError:
                number_unfilled_dt += 1
                continue

        sum_dt_column /= (len(tile_ids) - number_unfilled_dt)

        for id_index in tile_ids:
            try:
                dt_phi_rel[id_index] -= sum_dt_column
            except KeyError:
                continue

    return dt_phi_rel


# ---------------------------------------
def loop_correction_z(detector, dt_phi_rel: dict, dt_z_rel: dict, station: int):
    for phi in range(len(detector.TileDetector.column_ids(0, station)) - 1):
        # tile_ids_row_1 = detector.TileDetector.row_ids(z, station)
        # tile_ids_row_2 = detector.TileDetector.row_ids(z + 1, station)

        sum_dt_row = 0.
        number_unfilled_dt = 0
        counter = 0
        for z in range(len(detector.TileDetector.row_ids(0, station)) - 1):
            tile_id_1 = detector.TileDetector.id_from_row_col(row=phi, column=z, station_offset=station)
            tile_id_2 = detector.TileDetector.id_from_row_col(row=phi, column=z, station_offset=station)
            try:
                sum_dt_row += dt_z_rel[tile_id_1] - dt_z_rel[tile_id_2]
            except KeyError:
                number_unfilled_dt += 1
                continue
            counter += 1
        tile_id_1 = detector.TileDetector.id_from_row_col(row=phi, column=0, station_offset=station)
        tile_id_2 = detector.TileDetector.id_from_row_col(row=phi, column=51, station_offset=station)
        sum_dt_row += dt_phi_rel[tile_id_1] - dt_phi_rel[tile_id_2]

        counter += 1
        sum_dt_row /= (2 * counter - 2 * number_unfilled_dt)

        # correcting
        for z in range(len(detector.TileDetector.row_ids(0, station)) - 1):
            tile_id_1 = detector.TileDetector.id_from_row_col(row=phi, column=z, station_offset=station)
            tile_id_2 = detector.TileDetector.id_from_row_col(row=phi, column=z, station_offset=station)
            try:
                dt_z_rel[tile_id_1] += sum_dt_row
                dt_z_rel[tile_id_2] -= sum_dt_row
            except KeyError:
                number_unfilled_dt += 1
                continue
            counter += 1
        tile_id_1 = detector.TileDetector.id_from_row_col(row=phi, column=0, station_offset=station)
        tile_id_2 = detector.TileDetector.id_from_row_col(row=phi, column=51, station_offset=station)
        dt_phi_rel[tile_id_1] += sum_dt_row
        dt_phi_rel[tile_id_2] -= sum_dt_row

    return dt_phi_rel, dt_z_rel
