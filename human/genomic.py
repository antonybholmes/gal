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
from ... import text
from ... import genomic

PATH = __path__
REFSEQ_FILE = f"{__path__}/assets/human/grch38/ucsc_refseq_exons_entrez_grch38.v20220921.tsv"
TAD_FILE = f"{PATH}/assets/human/tads.gencode.v38lift37.genes.approved.tsv"
TANDEM_REPEATS_FILE = f"{PATH}/assets/human/grch38/simple_tandem_repeats.bed"
CHROM_SIZE_FILE = f"{PATH}/assets/human/grch38/chrom.sizes.txt"
ENCODE_BLACKLIST_FILE = f"{PATH}/assets/human/hg19/encode_blacklist.bed"
GIULIA_BLACKLIST_FILE = f"{PATH}/assets/human/hg19/rdf_giulia_blacklist_hg19.bed"
UCSC_CENTROMERES_FILE = f"{PATH}/assets/human/grch38/ucsc_centromeres.bed"
UCSC_TELOMERES_FILE = f"{PATH}/assets/human/grch38/ucsc_telomeres.bed"
RDF_PERICENTROMERES_FILE = f"{PATH}/assets/human/hg19/rdf_pericentromeres_hg19.bed"



class HumanChromosomes(genomic.Chromosomes):
    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        print("Initializing Human Chromosomes")
        super().__init__(CHROM_SIZE_FILE)
        

class HumanCentromeres(genomic.SearchGenomicBedFeatures):
    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        print("Initializing Human Centromeres")
        super().__init__(UCSC_CENTROMERES_FILE)


# class HumanPericentromeres(genomic.SearchGenomicBedFeatures):
#     _instance = None

#     def __new__(cls):
#         if cls._instance is None:
#             cls._instance = super().__new__(cls)
#         return cls._instance

#     def __init__(self):
#         print("Initializing Human Pericentromeres")
#         super().__init__(human.RDF_PERICENTROMERES_FILE)
        

class HumanTelomeres(genomic.SearchGenomicBedFeatures):
    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        print("Initializing Human Telomeres")
        super().__init__(UCSC_TELOMERES_FILE)


# class Telomeres(genomic.ClassifyRegion):
#     """
#     Determine whether a location overlaps a centromere
#     """

#     _instance = None

#     def __new__(cls):
#         if cls._instance is None:
#             cls._instance = super().__new__(cls)
#         return cls._instance

#     def __init__(self):
#         self._chromosomes = HumanChromosomes()

#     def get_classification(self, location: genomic.Location):
#         end = self._chromosomes.get_size(
#             location.chr) - genomic.TELOMERE_SIZE

#         if location.start <= genomic.TELOMERE_SIZE:
#             classification = "peri_telomeric"
#         elif location.end >= end:
#             classification = "peri_telomeric"
#         else:
#             classification = text.NA

#         return classification

class Telomeres(genomic.ClassifyRegion):
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        self._overlaps = genomic.GenomicFeaturesOverlap(HumanTelomeres())

    def get_classification(self, location: genomic.Location):
        classification = text.NA

        overlap = self._overlaps.get_max_overlap(location)

        if overlap is not None:
            p = overlap.length / location.length
            classification = "telomeric"

        return classification


class Centromeres(genomic.ClassifyRegion):
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        self.cen_overlaps = genomic.GenomicFeaturesOverlap(HumanCentromeres())
        #self.p_cen_overlaps = genomic.GenomicFeaturesOverlap(HumanPericentromeres())

    def get_classification(self, location: genomic.Location):
        classification = text.NA

        # Are we in a centromere
        in_centromere = False

        overlap = self.cen_overlaps.get_max_overlap(location)

        p = 0

        if overlap is not None:
            p = overlap.length / location.length

            classification = "centromeric"
            in_centromere = True

        #
        # Are we in a pericentromere
        #

        # if not in_centromere:
        #     overlap = self.p_cen_overlaps.get_max_overlap(location)

        #     if overlap is not None:
        #         p = overlap.length / location.length
        #         classification = "peri_centromeric"

        return classification


class Repetitive(genomic.Annotation):
    """
    Determine whether a location overlaps a repetitive region
    """

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__('repetitive')
        self._centromeres = Centromeres()
        self._telomeres = Telomeres()

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
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__('simple-tandem-repeats')
        trf = genomic.SearchGenomicBedFeatures(
            TANDEM_REPEATS_FILE)
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


class EncodeBlacklist(genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__('encode-black-list')
        self._trf_overlaps = genomic.GenomicFeaturesOverlap(genomic.SearchGenomicBedFeatures(ENCODE_BLACKLIST_FILE))

    def get_names(self):
        return ["ENCODE blacklist"]

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "encode_blacklist"
        else:
            classification = text.NA

        return [classification]


class GiuliaBlacklist(genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__('giulia-black-list')
        self._trf_overlaps = genomic.GenomicFeaturesOverlap(genomic.SearchGenomicBedFeatures(GIULIA_BLACKLIST_FILE))

    def get_names(self):
        return ["Giulia blacklist"]

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "giulia_blacklist"
        else:
            classification = text.NA

        return [classification]


# class Nnnn(genomic.Annotation):
#     """
#     Determine whether a location overlaps a centromere
#     """

#     def __init__(self):
#         nnnn = genomic.SearchGenomicBedFeatures(
#             "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/assembly/hg19/nnnn_hg19.bed")

#         self._nnnn_overlaps = genomic.GenomicFeaturesOverlap(nnnn)

#     def get_names(self):
#         return ["NNNNs"]

#     def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
#         overlap = self._nnnn_overlaps.get_max_overlap(location)

#         if overlap is not None:
#             classification = "nnnn"
#         else:
#             classification = text.NA

#         return [classification]