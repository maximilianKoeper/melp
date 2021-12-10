import dataclasses


@dataclasses.dataclass
class Cluster:
    id: int
    hits: list = dataclasses.field(default_factory=list)
    time: float = None

    def __eq__(self, other):
        return self.id == other.id

    def __str__(self):
        print(f'Cluster: ID={self.id}, # Hits={len(self.hits)}')

    def __len__(self):
        return len(self.hits)
    