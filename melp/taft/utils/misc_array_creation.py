import numpy as np


def create_array_from_dict(detector, dt_dict_z: dict, dt_dict_phi: dict, station: int) -> np.array:
    number_columns = len(detector.TileDetector.row_ids(0, 200000))
    number_rows = len(detector.TileDetector.column_ids(0, 200000))

    column_len_dict = number_columns
    row_len_dict = number_rows

    print(column_len_dict, row_len_dict)

    result = np.zeros((2, number_columns, number_rows))
    for column_index in range(column_len_dict):
        for row_index in range(row_len_dict):
            print(column_index, row_index)
            tile_id = detector.TileDetector.id_from_row_col(column=column_index, row=row_index, station_offset=station)
            result[0][column_index][row_index] = dt_dict_phi[tile_id]
            if column_index < 51:
                result[1][column_index][row_index] = dt_dict_z[tile_id]

    return result
