#---------------------------------------------------------------------
#  TILE CLASS
#       - pos (position)
#       - dir (normal vector to tile surface)
#       - id
#---------------------------------------------------------------------


class Tile():
    def __init__ (self, tile_pos, tile_dir, tile_id):
        self.pos = tile_pos
        self.dir = tile_dir
        self.id  = tile_id
