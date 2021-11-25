# ---------------------------------------
def loop_correction_phi(__detector__, dt_phi_rel: dict, station: int):
    for z in range(len(__detector__.TileDetector.row_ids(0, station))):
        tile_ids = __detector__.TileDetector.column_ids(z, station)

        sum_dt_row = 0.
        number_unfilled_dt = 0
        for id_index in tile_ids:
            try:
                sum_dt_row += dt_phi_rel[id_index]
            except KeyError:
                number_unfilled_dt += 1
                continue

        sum_dt_row /= (len(tile_ids) - number_unfilled_dt)

        for id_index in tile_ids:
            try:
                dt_phi_rel[id_index] -= sum_dt_row
            except KeyError:
                continue

    return dt_phi_rel
