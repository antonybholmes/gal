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
from .. import genes
from .. import mir
from .. import text
from .. import genomic
from .. import centromeres
from .. import peaks
from .. import tss
from .. import core
from .. import tad
from .. import db
from .. import bed
from .. import species
import os

# REFSEQ_FILE = os.path.join(
#    core.PATH, "assets/human/hg19/ucsc_refseq_exons_entrez_hg19.v20220921.tsv")
REFSEQ_FILE = os.path.join(
    core.PATH, "assets/mouse/mm10/ucsc_refseq_exons_entrez_mm10.v20231129.txt"
)
MIR_BED_FILE = os.path.join(core.PATH, "assets/human/hg19/mirbase/v20/mir.bed")
TAD_FILE = os.path.join(
    core.PATH, "assets/human/hg19/gcb_tads.gencode.v38lift37.genes.approved.tsv"
)
TANDEM_REPEATS_FILE = os.path.join(
    core.PATH, "assets/human/hg19/simple_tandem_repeats_hg19.bed"
)
CHROM_SIZE_FILE = os.path.join(core.PATH, "assets/mouse/mm10/chrom.sizes.txt")
UCSC_CENTROMERES_FILE = os.path.join(
    core.PATH, "assets/human/hg19/ucsc_centromeres_hg19.bed"
)
UCSC_TELOMERES_FILE = os.path.join(core.PATH, "assets/human/hg19/ucsc_telomeres.bed")

ENCODE_BLACKLIST_FILE = os.path.join(
    core.PATH, "assets/human/hg19/encode_blacklist.bed"
)

GIULIA_BLACKLIST_FILE = os.path.join(
    core.PATH, "assets/human/hg19/rdf_giulia_blacklist_hg19.bed"
)


class MouseRefSeqGenes(genes.RefSeqGenes):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__(REFSEQ_FILE)


REFSEQ_GENES = MouseRefSeqGenes()


class MouseAnnotateGene(peaks.AnnotateGene):
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


class MouseChromosomes(genomic.Chromosomes):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        print("Initializing Mouse Chromosomes")
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


class MouseOverlapTss(tss.OverlapTss):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self, block_size=100):
        super().__init__(REFSEQ_FILE, block_size)


class EncodeBlacklist(genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        super().__init__("encode-black-list")
        self._trf_overlaps = genomic.GenomicFeaturesOverlap(
            bed.SearchGenomicBedFeatures(ENCODE_BLACKLIST_FILE)
        )

    def get_names(self):
        return ["ENCODE blacklist"]

    def update_row(
        self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]
    ):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "encode_blacklist"
        else:
            classification = text.NA

        return [classification]


class GencodeTADAnnotation(tad.TADAnnotation):
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self, bin_size=100):
        super().__init__(TAD_FILE, bin_size=bin_size)


class MM10Version(db.DbVersion):
    def __init__(self):
        super().__init__("Mouse", "mm10", "RefSeq")


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

        self.add_module(MM10Version())

        self.add_module(
            MouseAnnotateGene(
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

        self.add_module(
            peaks.NClosestGenes(
                REFSEQ_GENES, refseq_start, n=n_closest, species=species.Species.MOUSE
            )
        )

        # human.REFSEQ_GENES, prom_ext_5p, bin_size))
        # self.add_module(mir.MirAnnotation(HumanMirs()))
        # self.add_module(HumanRepetitive())
        # self.add_module(HumanSimpleTandemRepeats())
        # self.add_module(EncodeBlacklist())
        # self.add_module(GiuliaBlacklist())
        # self.add_module(GencodeTADAnnotation())
        # self.add_module(tad.IsTADAnnotation())
        # self.add_module(tad.IsClosestTADAnnotation())
        self.add_module(MouseOverlapTss())
