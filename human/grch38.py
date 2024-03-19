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
from typing import Mapping, Union
from .. import genes
from .. import mir
from .. import genomic
from .. import centromeres
from .. import peaks
from .. import tss
from .. import core
from .. import db
from .. import tad
from .. import bed
from .. import text
import os

REFSEQ_FILE = os.path.join(
    core.PATH, "assets/human/grch38/ucsc_refseq_exons_entrez_grch38.v20220921.tsv"
)
MIR_BED_FILE = os.path.join(core.PATH, "assets/human/grch38/mirbase/v22/mir.bed")
TAD_FILE = os.path.join(core.PATH, "assets/human/grch38/gcb_tads_approved_20221109.tsv")
TANDEM_REPEATS_FILE = os.path.join(
    core.PATH, "assets/human/grch38/simple_tandem_repeats.bed"
)
CHROM_SIZE_FILE = os.path.join(core.PATH, "assets/human/grch38/chrom.sizes.txt")
UCSC_CENTROMERES_FILE = os.path.join(
    core.PATH, "assets/human/grch38/ucsc_centromeres.bed"
)
UCSC_TELOMERES_FILE = os.path.join(core.PATH, "assets/human/grch38/ucsc_telomeres.bed")
ENCODE_BLACKLIST_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/encode/blacklists/grch38/GRCh38_unified_blacklist.bed"


class HumanRefSeqGenes(genes.RefSeqGenes):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__(REFSEQ_FILE)


REFSEQ_GENES = HumanRefSeqGenes()

class EncodeBlacklist(genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        super().__init__('encode-black-list')
        self._trf_overlaps = genomic.GenomicFeaturesOverlap(
            bed.SearchGenomicBedFeatures(ENCODE_BLACKLIST_FILE))

    def get_names(self):
        return ["ENCODE blacklist GRCh38"]

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "encode_blacklist"
        else:
            classification = text.NA

        return [classification]

class HumanAnnotateGene(peaks.AnnotateGene):
    def __init__(
        self, promoter_type: str, prom_ext_5p: int, prom_ext_3p: int, bin_size: int
    ):
        super().__init__(
            promoter_type,
            tss.RefSeqAnnotationFactory.getInstance(
                REFSEQ_FILE, prom_ext_5p, prom_ext_3p, bin_size
            ),
            REFSEQ_GENES,
        )


class HumanMirs(mir.SearchMirs):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        print("Initializing mIRs")
        super().__init__(MIR_BED_FILE)


class HumanChromosomes(genomic.Chromosomes):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        print("Initializing Human Chromosomes")
        super().__init__(CHROM_SIZE_FILE)


class HumanRepetitive(centromeres.Repetitive):
    """
    Determine whether a location overlaps a repetitive region
    """

    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__(
            bed.SearchGenomicBedFeatures(UCSC_CENTROMERES_FILE),
            bed.SearchGenomicBedFeatures(UCSC_TELOMERES_FILE),
        )


class HumanSimpleTandemRepeats(centromeres.SimpleTandemRepeats):
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__(bed.SearchGenomicBedFeatures(TANDEM_REPEATS_FILE))


class HumanOverlapTss(tss.OverlapTss):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self, block_size=100):
        super().__init__(REFSEQ_FILE, block_size)


class GencodeTADAnnotation(tad.TADAnnotation):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self, bin_size=100):
        super().__init__(TAD_FILE, bin_size=bin_size)


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


class GRCh38Version(db.DbVersion):
    def __init__(self):
        super().__init__("Human", "GRCh38", "RefSeq")


class AnnotatePeak(peaks.AnnotatePeak):
    """
    Core annotation for annotating peaks/regions
    """

    def __init__(
        self,
        type,
        prom_ext_5p: int = 5000,
        prom_ext_3p: int = 4000,
        bin_size: int = 10000,
        n_closest: int = 5,
    ):
        super().__init__(type)
        # annotations specific to human

        # skip version module
        # self.add_module(GRCh38Version())

        self.add_module(
            HumanAnnotateGene(
                genomic.promoter_type(prom_ext_5p, prom_ext_3p),
                prom_ext_5p,
                prom_ext_3p,
                bin_size,
            )
        )

        # default closest gene
        refseq_start = tss.RefSeqStart(
            REFSEQ_FILE, REFSEQ_GENES, prom_ext_5p, prom_ext_3p
        )

        self.add_module(refseq_start)

        self.add_module(
            tss.RefSeqEnd(REFSEQ_FILE, REFSEQ_GENES, prom_ext_5p, prom_ext_3p)
        )

        # since we already add the closest, start at 2 and
        self.add_module(
            peaks.NClosestGenes(REFSEQ_GENES, refseq_start, start=2, n=n_closest - 1)
        )

        # human.REFSEQ_GENES, prom_ext_5p, bin_size))
        self.add_module(mir.MirAnnotation(HumanMirs()))
        self.add_module(HumanRepetitive())
        self.add_module(HumanSimpleTandemRepeats())
        self.add_module(EncodeBlacklist())
        # self.add_module(hg.GiuliaBlacklist())
        self.add_module(GencodeTADAnnotation())
        self.add_module(tad.IsTADAnnotation())
        self.add_module(tad.IsClosestTADAnnotation())
        self.add_module(HumanOverlapTss())
