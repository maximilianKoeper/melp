# ---------------------------------------
def loop_correction_phi(detector, dt_phi_rel: dict, station: int):
    for z in range(len(detector.TileDetector.row_ids(0, station))):
        tile_ids = detector.TileDetector.column_ids(z, station)

        sum_dt_column = 0.
        for id_index in tile_ids:
            sum_dt_column += dt_phi_rel[id_index]

        sum_dt_column /= len(tile_ids)

        for id_index in tile_ids:
            dt_phi_rel[id_index] -= sum_dt_column

    return dt_phi_rel

# ---------------------------------------
# OLD VERSION
# def loop_correction_phi(detector, dt_phi_rel: dict, station: int):
#    for z in range(len(detector.TileDetector.row_ids(0, station))):
#        tile_ids = detector.TileDetector.column_ids(z, station)
#
#        sum_dt_column = 0.
#        number_unfilled_dt = 0
#        for id_index in tile_ids:
#            try:
#                sum_dt_column += dt_phi_rel[id_index]
#            except KeyError:
#                number_unfilled_dt += 1
#                continue
#
#        sum_dt_column /= (len(tile_ids) - number_unfilled_dt)
#
#        for id_index in tile_ids:
#            try:
#                dt_phi_rel[id_index] -= sum_dt_column
#            except KeyError:
#                continue
#
#    return dt_phi_rel
