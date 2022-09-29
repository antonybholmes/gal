# -*- coding: utf-8 -*-
"""
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.
This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with 
this program. If not, see <https://www.gnu.org/licenses/>. 

Copyright (C) 2022 Antony Holmes.
"""

from typing import Any, Mapping, Union
from . import genomic
from . import text


class Centromeres(genomic.ClassifyRegion):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self, centromeres: genomic.SearchGenomicBedFeatures, classification: str = "centromeric"):
        self._overlaps = genomic.GenomicFeaturesOverlap(centromeres)
        self._classification = classification

    def get_classification(self, location: genomic.Location):
        overlap = self._overlaps.get_max_overlap(location)

        if overlap is not None:
            p = overlap.length / location.length
            return self._classification
        else:
            return text.NA


class Telomeres(Centromeres):
    def __init__(self, telomeres: genomic.SearchGenomicBedFeatures):
        super().__init__(telomeres, "telomeric")


class Repetitive(genomic.Annotation):
    """
    Determine whether a location overlaps a repetitive region
    """

    def __init__(self, centromeres:genomic.SearchGenomicBedFeatures, telomeres:genomic.SearchGenomicBedFeatures):
        super().__init__('repetitive')
        self._centromeres = Centromeres(centromeres)
        self._telomeres = Telomeres(telomeres)

    def get_classification(self, location: genomic.Location):
        classifications = set()

        classifications.add(self._centromeres.get_classification(location))
        classifications.add(self._telomeres.get_classification(location))

        # If there are multiple classifications, get rid of the n/a
        if len(classifications) > 1:
            classifications.remove(text.NA)

        return sorted(classifications)

    def get_names(self):
        return ["Centromere/Telomere"]

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        ret = ",".join(sorted(self.get_classification(location)))

        return [ret]


class SimpleTandemRepeats(genomic.Annotation):
    def __init__(self, trf:genomic.SearchGenomicBedFeatures):
        super().__init__('simple-tandem-repeats')
        self._trf_overlaps = genomic.GenomicFeaturesOverlap(trf)

    def get_names(self):
        return ["Simple Tandem Repeats"]

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "tandem_repeat"
        else:
            classification = text.NA

        return [classification]