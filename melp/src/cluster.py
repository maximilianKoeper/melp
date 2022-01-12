import dataclasses


#cluster hit class
@dataclasses.dataclass
class ClusterHit:
    tile_id: int
    #edep: float = 0
    mc_i: int = 0
    tid: int = -1
    frame_id: int = -1
    primary: int = 0
    time: float = -1
    #run_id: int = -1
    #hid: int = 0
    #impact_vec: list = None
    #trajectory: object = None
    #pos: list = None

    # -----------------------------------------
    #  public functions
    # -----------------------------------------

    #returns True if hits are the same
    #def __eq__(self, other):
    #    return self.mc_i == other.mc_i

    #returns 


#cluster class
@dataclasses.dataclass
class Cluster:
    id: int
    frame_id: int
    master_id: int
    master_primary: int = 0
    hits: list = dataclasses.field(default_factory=list)
    
    # -----------------------------------------
    #  public functions
    # -----------------------------------------

    #returns True if the two clusters are the same
    def __eq__(self, other):
        return self.id == other.id

    #returns the cluster id and the number of hits in the cluster
    #def __str__(self):
    #    return (f'Cluster: ID={self.id}, # Hits={len(self.hits)}')

    #returns the number of hits in the cluster
    def __len__(self):
        return len(self.hits)
    
    #returns all tile ids in cluster
    def get_tile_ids(self):
        tile_ids = []
        for hit in self.hits:
            tile_ids.append(hit.tile_id)
        return tile_ids

    #returns all primary ids in cluster
    def get_primaries(self):
        primary_ids = []
        for hit in self.hits:
            primary_ids.append(hit.primary)
        return primary_ids
    
    #returns all times in cluster
    def get_times(self):
        times = []
        for hit in self.hits:
            times.append(hit.time)
        return times

    #return all frame ids in cluster
    def get_frame_ids(self):
        frame_ids = []
        for hit in self.hits:
            frame_ids.append(hit.frame_id)
        return frame_ids

    #return the value and index of hit with smallest time in cluster
    def get_min_time(self):
        times = []
        for hit in self.hits:
            times.append(hit.time) 
        val, idx = min((val, idx) for (idx, val) in enumerate(times))
        return self.hits[idx], idx
    


