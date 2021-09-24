#---------------------------------------------------------------------
#  HIT CLASS
#       - pos (position)
#       - id
#       - time
#       - edep
#---------------------------------------------------------------------


class Hit():
    def __init__ (self, time=0, edep=0, timestamp=0, mc_i=0, mc_n=0, primary=0):
        self.time       = time
        self.edep       = edep
        self.timestamp  = timestamp
        self.mc_i       = mc_i
        self.mc_n       = mc_n
        self.primary    = primary
