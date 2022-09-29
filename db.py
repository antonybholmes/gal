from typing import Mapping, Union
from . import genomic


class DbVersion(genomic.Annotation):
    def __init__(self, genome: str = "Human", assembly: str = "GRCh38", track: str = "RefSeq"):
        super().__init__('db-version')
        self._genome = genome
        self._assembly = assembly
        self._track = track
        self._db = "/".join([genome, assembly, track])

    def get_names(self):
        return ["Data Source"]

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        return [self._db]
