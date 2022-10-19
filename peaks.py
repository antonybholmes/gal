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

from abc import ABC
import json
import sys
import collections
import re
from typing import Any, Mapping, Optional, Union
import numpy as np
import pandas as pd

from . import text
from . import genomic
from . import headings
from . import genes
from . import tss
from . import species

BIN_SIZE = 10000
PEAKS_PER_GENE_HEADING = 'Peaks Per Gene'
LOCATION_HEADING = 'Genomic Location'
PVALUE_HEADING = 'P-value (ChIPseeqer)'
SCORE_HEADING = 'Score (ChIPseeqer)'
WIDTH_HEADING = 'Peak Width'


def parse_peaks(file):
    peaks = collections.defaultdict(
        lambda: collections.defaultdict(lambda: (int, int)))

    f = open(file, 'r')

    # skip header
    # f.readline()

    for line in f:
        tokens = line.split('\t')

        chr = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])

        id = ':'.join([chr, '-'.join([str(start), str(end)])])

        #sys.stderr.write(id + '\n')

        peaks[chr][id] = (start, end)

    f.close()

    return peaks





def duplicate_peaks(type, file):
    f = open(file, 'r')

    header = f.readline().strip().split('\t')

    entrez_column = text.find_index(
        header, headings.ENTREZ_ID)
    refseq_column = text.find_index(
        header, headings.REFSEQ_ID)
    symbol_column = text.find_index(
        header, headings.GENE_SYMBOL)
    overlap_type_column = text.find_index(
        header, 'Relative To Gene')
    tss_column = text.find_index(header, 'TSS Distance')

    #mir_column = text.find_index(header, 'miR Symbol')
    #mir_type_column = text.find_index(
    #    header, 'Relative To miR')
    #mss_column = text.find_index(
    #    header, 'miR Start Distance')

    print('\t'.join(header))

    for line in f:
        ls = line.strip()

        if len(ls) == 0:
            continue

        tokens = ls.split('\t')

        entrezes = tokens[entrez_column].split(';')
        refseqs = tokens[refseq_column].split(';')
        overlap_types = tokens[overlap_type_column].split(';')
        symbols = tokens[symbol_column].split(';')
        tsses = tokens[tss_column].split(';')

        split = False

        for i in range(0, len(refseqs)):
            refseq = refseqs[i]

            if refseq == text.NA:
                continue

            entrez = entrezes[i]
            overlap_type = overlap_types[i]
            symbol = symbols[i]
            tss = tsses[i]

            # clone
            new_tokens = tokens[:]

            new_tokens[overlap_type_column] = overlap_type
            new_tokens[entrez_column] = entrez
            new_tokens[refseq_column] = refseq
            new_tokens[symbol_column] = symbol
            new_tokens[tss_column] = tss

            #new_tokens[mir_column] = text.NA
            #new_tokens[mir_type_column] = text.NA

            print('\t'.join(new_tokens))

            split = True

        # Only process the mirs if they exist

        # if mir_column != -1:
        #     mirs = tokens[mir_column].split(';')
        #     mir_types = tokens[mir_type_column].split(';')
        #     mir_mss = tokens[mss_column].split(';')

        #     for i in range(0, len(mirs)):
        #         mir = mirs[i]

        #         if mir == text.NA:
        #             continue

        #         mir_type = mir_types[i]
        #         mss = mir_mss[i]

        #         new_tokens = tokens[:]

        #         new_tokens[overlap_type_column] = text.NA
        #         new_tokens[entrez_column] = text.NA
        #         new_tokens[refseq_column] = text.NA
        #         new_tokens[symbol_column] = text.NA
        #         new_tokens[tss_column] = text.NA

        #         new_tokens[mir_column] = mir
        #         new_tokens[mir_type_column] = mir_type
        #         new_tokens[mss_column] = mss

        #         print('\t'.join(new_tokens))

        #         split = True

        if split == False:
            # no lines were split so print the row as it was originally
            print('\t'.join(tokens))

    f.close()


def filter_peaks(file, max_tss_5p_dist, max_tss_3p_dist):
    """
    Exclude peaks that are too far from a tss
    """

    f = open(file, 'r')

    header = f.readline().strip().split('\t')

    tss_column = text.find_index(
        header, headings.TSS_DISTANCE)

    print('\t'.join(header))

    for line in f:
        line = line.strip()

        if len(line) == 0:
            continue

        tokens = line.split('\t')

        tss = tokens[tss_column]

        if tss == text.NA:
            continue

        d = int(tss)

        if d < max_tss_5p_dist or d > max_tss_3p_dist:
            continue

        print('\t'.join(tokens))

    f.close()


class NearestPeak:
    """
    Finds the closest peak to another set of peaks.
    """

    def __init__(self, file=None):
        self._peak_start_map = collections.defaultdict(
            lambda: collections.defaultdict(genomic.Location))
        self._starts = collections.defaultdict(list)

        if file is not None:
            self.load(file)

    def load(self, file):
        print('Loading peaks from {file}...', file=sys.stderr)

        f = open(file, 'r')

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')

            location = genomic.parse_location(tokens[0])

            mid_point = genomic.mid_point(location)

            self._peak_start_map[location.chr][mid_point] = location

        f.close()

        for chr in self._peak_start_map:
            self.starts[chr] = sorted(self.peak_starts[chr])

        print('Finished.', file=sys.stderr)

    def load_from_cols(self, file):
        print('Loading peaks from {file}...', file=sys.stderr)

        f = open(file, 'r')

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')
            location = genomic.parse_location_cols(tokens)
            mid_point = genomic.mid_point(location)
            self._peak_start_map[location.chr][mid_point] = location

        f.close()

        for chr in self._peak_start_map:
            self.starts[chr] = sorted(self._peak_start_map[chr])

        print('Finished.', file=sys.stderr)

    def get_nearest_peak(self, location: genomic.Location):
        chr = location.chr

        if len(self.starts[chr]) == 0:
            return None

        mid_point = genomic.mid_point(location)

        # use a binary style search to find closest peak

        ps = 0
        pe = len(self._starts[chr]) - 1

        while pe - ps > 1:
            pm = int((pe + ps) / 2)
            test_start = self._starts[chr][pm]

            #sys.stderr.write(str(test_start) + ' ' + str(mid_point) + ' ' + str(ps) + ' ' + str(pe) + ' ' + str(pm) + '\n')

            # perfect match
            if test_start > mid_point:
                pe = pm  # - 1
            elif test_start < mid_point:
                ps = pm  # + 1
            else:
                return self._peak_start_map[chr][mid_point]

        # point lies between two midpoints, so pick the closest

        sd = mid_point - self._starts[chr][ps]
        ed = mid_point - self._starts[chr][pe]

        abss = abs(sd)
        abse = abs(ed)

        if (abss <= abse):
            start = self._starts[chr][ps]
        else:
            start = self._starts[chr][pe]

        return self._peak_start_map[chr][start]

    def get_nearest_abs_dist(self, location):
        return abs(self.get_nearest_dist(location))

    def get_nearest_dist(self, location):
        """
        Find the closest distance to a feature
        """

        chr = location.chr

        if len(self._starts[chr]) == 0:
            return sys.maxint

        mid_point = genomic.mid_point(location)

        # use a binary style search to find closest peak

        ps = 0
        pe = len(self._starts[chr]) - 1

        while pe - ps > 1:
            pm = int((pe + ps) / 2)
            test_start = self._starts[chr][pm]

            #sys.stderr.write(str(test_start) + ' ' + str(mid_point) + ' ' + str(ps) + ' ' + str(pe) + ' ' + str(pm) + '\n')

            # perfect match
            if test_start > mid_point:
                pe = pm  # - 1
            elif test_start < mid_point:
                ps = pm  # + 1
            else:
                return 0

        # point lies between two midpoints, so pick the closest

        # sys.stderr.write('hmm ' + chr + ' ' + str(ps) + '\n')

        sd = mid_point - self._starts[chr][ps]
        ed = mid_point - self._starts[chr][pe]

        abss = abs(sd)
        abse = abs(ed)

        if (abss <= abse):
            return sd
        else:
            return ed


class NClosestGenes(genomic.Annotation):
    """
    Find the n closest (usually 5) closest genes to a location
    and add annotation columns
    """

    def __init__(self,
                 refseq_genes: tss.RefSeqAnnotation,
                 refseq_start: tss.RefSeqTss,
                 species: species.Species = species.Species.HUMAN,
                 n: int = 5):
        super().__init__('n-closest-genes', species=species)
        self._refseq_genes = refseq_genes
        self._refseq_start = refseq_start
        self._n: int = n
        self._ncols = n * 4
        self._header = []
        for i in range(1, self._n + 1):
            self._header.append(f'#{i} Closest {headings.REFSEQ_ID}')
            self._header.append(f'#{i} Closest {headings.ENTREZ_ID}')
            self._header.append(
                f'#{i} Closest {headings.GENE_SYMBOL}')
            self._header.append(f'#{i} TSS Closest Distance')

    def get_names(self):
        return self._header

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
        # find a large panel of closest genes so we can hopefully find
        # n unique genes within this
        tss_genes = self._refseq_start.get_n_closest_genes(location, n=50)

        used = set()
        n = 1

        # find n unique entrez ids since there can be multiple refseqs with
        # same gene id

        ret = []

        for tss_gene in tss_genes:
            gene_annotation = self._refseq_genes.get_gene(tss_gene.id)

            if gene_annotation.entrez in used:
                continue

            used.add(gene_annotation.entrez)

            ret.append(gene_annotation.refseq)
            ret.append(gene_annotation.entrez)
            ret.append(gene_annotation.symbol)
            ret.append(str(tss_gene.d))

            # keep track of how many unique genes we've found

            if isinstance(row_map, dict):
                row_map[f'#{n} Closest {headings.REFSEQ_ID}'] = gene_annotation.refseq
                row_map[f'#{n} Closest {headings.ENTREZ_ID}'] = gene_annotation.entrez
                row_map[f'#{n} {headings.CLOSEST_GENE_SYMBOL}'] = gene_annotation.symbol
                row_map[f'#{n} TSS Closest Distance'] = tss_gene.d

            if n == self._n:
                break

            n += 1

        while len(ret) < self._ncols:
            ret.append(text.NA)

        #print(ret, len(ret), len(self._header))
        return ret


# class ClosestGene(genomic.Annotation):
#     """
#     Closest gene to a location
#     """

#     def __init__(self,
#                  promoter_type: str,
#                  refseq_genes: tss.RefSeqAnnotation,
#                  refseq_tss: tss.RefSeqTss):
#         super().__init__('closest-gene')
#         self._promoter_type = promoter_type
#         self._refseq_genes: tss.RefSeqAnnotation = refseq_genes
#         self._refseq_tss: tss.RefSeqTss = refseq_tss

#     def get_names(self):
#         ret = []
#         ret.append(f'Closest {headings.REFSEQ_ID}')
#         ret.append(f'Closest {headings.ENTREZ_ID}')
#         ret.append(headings.CLOSEST_GENE_SYMBOL)
#         ret.append(
#             f'Relative To Closest Gene ({self._promoter_type})')
#         ret.append(f'TSS Closest Distance')
#         return ret

#     def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
#         # find a large panel of closest genes so we can hopefully find
#         # n unique genes within this
#         ret: list[str] = []

#         tss_gene = self._refseq_tss.get_closest_gene(location)

#         gene_annotation = self._refseq_genes.get_gene(tss_gene.id)

#         ret.append(gene_annotation.refseq)
#         ret.append(gene_annotation.entrez)
#         ret.append(gene_annotation.symbol)
#         ret.append(reversed(sorted(tss_gene.types)))
#         ret.append(str(tss_gene.d))

#         row_map[f'Closest {headings.REFSEQ_ID}'] = gene_annotation.refseq
#         row_map[f'Closest {headings.ENTREZ_ID}'] = gene_annotation.entrez
#         row_map[headings.CLOSEST_GENE_SYMBOL] = gene_annotation.symbol
#         row_map[f'Relative To Closest Gene ({self._promoter_type})'] = reversed(
#             sorted(tss_gene.types))
#         row_map[f'TSS Closest Distance'] = tss_gene.d

#         return ret


class AnnotateGene(genomic.Annotation):
    """
    Basic gene info
    """

    def __init__(self,
                 promoter_type: str,
                 refseq_annotation: tss.RefSeqAnnotation,
                 refseq_genes: genes.RefSeqGenes):
        super().__init__('annotate-gene')
        self._promoter_type = promoter_type
        self._refseq_annotation = refseq_annotation
        self._refseq_genes = refseq_genes

        self._header = []
        self._header.append(headings.REFSEQ_ID)
        self._header.append(headings.ENTREZ_ID)
        self._header.append(headings.GENE_SYMBOL)
        self._header.append(f'Relative To Gene ({self._promoter_type})')
        self._header.append(headings.TSS_DISTANCE)
        self._empty_row = [text.NA] * 5

        self._metadata.update(refseq_annotation.metadata)
        self._metadata['promoter'] = promoter_type

    def get_names(self):
        return self._header

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
        genes = self._refseq_annotation.annotate_location(location)
        row_map['genes'] = genes

        refseqs = []
        symbols = []
        entrezes = []
        chrs = []
        strands = []
        starts = []
        ends = []
        types = []
        tss_list = []

        for gene in genes:
            gene_annotation = self._refseq_genes.get_gene(gene.id)

            chrs.append(gene_annotation.chr)
            starts.append(gene_annotation.start)
            strands.append(gene_annotation.strand)
            ends.append(gene_annotation.end)
            entrezes.append(gene_annotation.entrez)
            refseqs.append(gene_annotation.refseq)
            symbols.append(gene_annotation.symbol)

            # include tss distance for all peaks regardless
            tss_list.append(str(gene.d))

            types.append(','.join(reversed(sorted(gene.types))))

        row = []

        if len(entrezes) > 0:
            row.append(';'.join(refseqs))
            row.append(';'.join(entrezes))
            row.append(';'.join(symbols))
            row.append(';'.join(types))
            row.append(';'.join(tss_list))
            row_map[headings.GENE_SYMBOL] = ';'.join(symbols)
        else:
            row.extend(self._empty_row)
            row_map[headings.GENE_SYMBOL] = text.NA

        return row


class TableAnnotation(ABC):
    def annotate(self, header: list[str], df: pd.DataFrame) -> pd.DataFrame:
        """
        Should return the name of the column or an array of columns
        """
        ...


class PeaksPerGeneAnnotation(TableAnnotation):
    """
    Flag when a gene has more than 3 peaks
    """

    def annotate(self, header: list[str], df: pd.DataFrame) -> pd.DataFrame:
        header = np.array(header)

        loc_col = text.find_index(header, LOCATION_HEADING)
        # np.where(header == headings.GENE_SYMBOL)[0][0]
        gene_col = text.find_index(header, headings.GENE_SYMBOL)

        peaks_genes = collections.defaultdict(set)
        genes_peaks = collections.defaultdict(set)

        for i in range(0, df.shape[0]):
            loc = df.iloc[i, loc_col]
            genes = df.iloc[i, gene_col].split(';')

            for gene in genes:
                if gene != text.NA:
                    peaks_genes[loc].add(gene)
                    genes_peaks[gene].add(loc)

        # allocate array to save using append
        counts = np.zeros(df.shape[0], dtype=int)

        for i in range(0, df.shape[0]):
            loc = df.iloc[i, loc_col]

            c = 0

            # see which genes we are on
            if loc in peaks_genes:
                genes = peaks_genes[loc]

                # see how many peaks are on this gene
                for gene in genes:
                    if gene in genes_peaks:
                        c = max(c, len(genes_peaks[gene]))

            # peaks per gene
            counts[i] = c

        df[PEAKS_PER_GENE_HEADING] = counts

        return df


class AnnotatePeak:
    """
    Core annotation for annotating peaks/regions
    """

    def __init__(self,
                 type: str):

        self._header: list[str] = []

        self._type = type
        #self.prom_ext_5p = prom_ext_5p
        #self.prom_ext_3p = prom_ext_3p
        #self.bin_size = 10000

        #self.mir_max_gap = prom_ext_5p
        # pychipseq.human.tss.RefSeqAnnotation(prom_ext_5p, prom_ext_3p, self.bin_size)

        #self._refseq_annotation = refseq_annotation
        # self._refseq_genes = refseq_genes  # RefSeqGenes()

        # pychipseq.human.tss.RefSeqTss(prom_ext_5p, prom_ext_3p)
        #self._refseq_start = refseq_start
        # pychipseq.human.tss.RefSeqEnd(prom_ext_5p, prom_ext_3p)
        #self._refseq_end = refseq_end

        # For annotating if in centromere or not
        #self.repetitive = genomic.Repetitive()

        #
        # Create a mir database
        #
        #self.mir_annotation = pychipseq.human.mir.MirAnnotation(bin_size)
        #self.mir_transcript_annotation = pychipseq.human.mir.MirTranscriptAnnotation()

        # self._promoter_type = genomic.promoter_heading f'(prom=-{prom_ext_5p / 1000}/+{prom_ext_3p / 1000}kb)'

        #
        # Modules are self contained annotation blocks to make it easier
        # to update them
        #

        self._annotation_modules: list[genomic.Annotation] = []
        self._table_annotation_modules: list[TableAnnotation] = []

        # self._annotation_modules.append(AnnotateGene(
        #    self._promoter_type, refseq_annotation, refseq_genes))

        # self._annotation_modules.append(pychipseq.human.genomic.SimpleTandemRepeats())
        # self._annotation_modules.append(pychipseq.human.genomic.Nnnn())

        # self._annotation_modules.append(ClosestGene(self._type,
        #     self._promoter_type,
        #     refseq_genes,
        #     refseq_tss))

        # self._annotation_modules.append(
        #    NClosestGenes(refseq_genes, refseq_start, n=n_closest))

        self._table_annotation_modules.append(PeaksPerGeneAnnotation())

        # if annotating multiple files, ensures the header construction
        # only happens once
        self._update_header = True

    def add_module(self, module: genomic.Annotation):
        self._annotation_modules.append(module)

    def print_header(self) -> list[str]:
        ret = self._header
        sys.stdout.write('\t'.join(ret))
        return ret

    @property
    def metadata(self) -> str:
        return {'modules': [module.metadata for module in self._annotation_modules]}

    def write_metadata(self):
        with open('metadata.json', 'w') as f:
            json.dump(self.metadata, f, indent=2)

    def print_metadata(self):
        print(json.dumps(self.metadata))

    @property
    def header(self) -> list[str]:
        if len(self._header) == 0:
            self._header.append(LOCATION_HEADING)
            self._header.append(PVALUE_HEADING)
            self._header.append(SCORE_HEADING)
            self._header.append(WIDTH_HEADING)

            # self._header.append(headings.REFSEQ_ID)
            # self._header.append(headings.ENTREZ_ID)
            # self._header.append(headings.GENE_SYMBOL)
            # self._header.append(f'{self._type} Relative To Gene {self._promoter_type}')
            # self._header.append(f'{self._type} {headings.TSS_DISTANCE}')

            #
            # Closest
            #
            # self._header.append(f'Closest {headings.REFSEQ_ID}')
            # self._header.append(f'Closest {headings.ENTREZ_ID}')
            # self._header.append(headings.CLOSEST_GENE_SYMBOL)
            # self._header.append(
            #     f'{self._type} Relative To Closest Gene {self._promoter_type}')
            # self._header.append(f'{self._type} TSS Closest Distance')

            #
            # Closest End
            #

            # self._header.append(f'Closest End {headings.REFSEQ_ID}')
            # self._header.append(f'Closest End {headings.ENTREZ_ID}')
            # self._header.append(f'Closest End {headings.GENE_SYMBOL}')
            # self._header.append(
            #     f'{self._type} Relative To Closest Gene End {self._promoter_type}')
            # self._header.append(f'{self._type} End Closest Distance')

            #
            # Other
            #

            #self._header.append(self.type + ' Centromere Telomere')
            #self._header.append('miR Symbol')
            #self._header.append(f'{self._type} Relative To miR')
            #self._header.append(f'{self._type} miR Start Closest Distance')
            #self._header.append(f'{self._type} miR Start Distance')

            #
            # Modules
            #

            for module in self._annotation_modules:
                self._header.extend(module.get_names())

            # for module in self._table_annotation_modules:
            #    self._header.extend(module.get_names())

        #
        # 5 closest genes
        #

        # for i in range(1, 6):
        #     ret.append(f'#{i} Closest {headings.REFSEQ_ID}')
        #     ret.append(f'#{i} Closest {headings.ENTREZ_ID}')
        #     ret.append(f'#{i} Closest {headings.GENE_SYMBOL}')
        #     ret.append(f'#{i} TSS Closest Distance')

        # End the header
        # sys.stdout.write('\n')

        return self._header

    def parse(self, file: str) -> pd.DataFrame:
        f = open(file, 'r')

        #
        # lets make a header
        #

        if 'TF_targets' in file:
            file_header = [LOCATION_HEADING, WIDTH_HEADING,
                           PVALUE_HEADING, SCORE_HEADING]
        elif 'seacr' in file:
            file_header = [LOCATION_HEADING,
                           WIDTH_HEADING, "Total Score", "Max Score"]
        elif 'overlap' in file:
            file_header = [LOCATION_HEADING,
                           WIDTH_HEADING, "Total Score", "Max Score", "Number Of Overlapping Peaks", "P1", "P2"]
            # skip header
            f.readline()
        elif 'bed' in file:
            file_header = [LOCATION_HEADING, WIDTH_HEADING]
        else:
            file_header = f.readline().strip().split('\t')

        if self._update_header:
            self._header.extend(file_header)

        # self._header.append(headings.REFSEQ_ID)
        # self._header.append(headings.ENTREZ_ID)
        # self._header.append(headings.GENE_SYMBOL)
        # self._header.append(f'Relative To Gene {self._promoter_type}')
        # self._header.append(f'{headings.TSS_DISTANCE}')

        #
        # Closest
        #
        # self._header.append(f'Closest {headings.REFSEQ_ID}')
        # self._header.append(f'Closest {headings.ENTREZ_ID}')
        # self._header.append(headings.CLOSEST_GENE_SYMBOL)
        # self._header.append(
        #     f'{self._type} Relative To Closest Gene {self._promoter_type}')
        # self._header.append(f'{self._type} TSS Closest Distance')

        #
        # Closest End
        #

        # self._header.append(f'Closest End {headings.REFSEQ_ID}')
        # self._header.append(f'Closest End {headings.ENTREZ_ID}')
        # self._header.append(f'Closest End {headings.GENE_SYMBOL}')
        # self._header.append(
        #     f'{self._type} Relative To Closest Gene End {self._promoter_type}')
        # self._header.append(f'{self._type} End Closest Distance')

        #
        # Other
        #

        #self._header.append(self.type + ' Centromere Telomere')
        # self._header.append('miR Symbol')
        # self._header.append(f'{self._type} Relative To miR')
        # self._header.append(f'{self._type} miR Start Closest Distance')
        # self._header.append(f'{self._type} miR Start Distance')

        #
        # Modules
        #

        if self._update_header:
            for module in self._annotation_modules:
                self._header.extend(module.get_names())

        #
        # Process some coordinates
        #

        locations = []
        annotations = []

        for line in f:
            tokens = line.strip().split('\t')

            first_col_is_loc = genomic.is_location(tokens[0])

            if first_col_is_loc:
                location = genomic.parse_location(tokens[0])
            else:
                location = genomic.Location(
                    tokens[0], int(tokens[1]), int(tokens[2]))

            #chr = tokens[0]

            # Skip chrM
            if 'chrM' in location.chr:
                continue

            locations.append(location)

            # supply already calculated items to each module
            row_map = {LOCATION_HEADING: location, 'type': self._type}

            #start = int(tokens[1])
            #end = int(tokens[2])
            width = location.end - location.start + 1

            if 'TF_targets' in file:
                # deal with infs
                if re.match(r'.*[Ii]nf.*', tokens[3]):
                    p = -1000
                else:
                    p = float(tokens[3])

                score = float(tokens[4])

                annotation = [str(location), width, p, score]

                row_map[PVALUE_HEADING] = p
                row_map[SCORE_HEADING] = score
                row_map[WIDTH_HEADING] = width
            elif 'seacr' in file:
                score = float(tokens[3])
                max_score = float(tokens[4])
                annotation = [str(location), width, score, max_score]
                row_map[WIDTH_HEADING] = width
                row_map["Total Score"] = score
                row_map["Max Score"] = max_score
            elif 'overlap' in file:
                score = float(tokens[2])
                max_score = float(tokens[3])
                overlaps = float(tokens[4])
                p1 = tokens[5]
                p2 = tokens[6]
                annotation = [str(location), width, score,
                              max_score, overlaps, p1, p2]
                row_map[WIDTH_HEADING] = width
                row_map["Total Score"] = score
                row_map["Max Score"] = max_score
                row_map["Overlaps"] = overlaps
                row_map["P1"] = p1
                row_map["P2"] = p2
            elif 'bed' in file:
                annotation = [str(location), width]
                row_map[WIDTH_HEADING] = width
            else:
                annotation = tokens[0:len(file_header), width]
                # annotation.append(location)

                # Use the columns as they are in the file
                for i in range(0, len(file_header)):
                    row_map[file_header[i]] = tokens[i]

            #max_height = int(tokens[6])

            # location = lib.genomic.Location(chr, start, end) #chr + ":" + str(start) + "-" + str(end)

            # Write the initial portion

            # Add the gene annotation
            annotation.extend(self.annotate(location, row_map=row_map))

            annotations.append(annotation)

        f.close()

        #header = self._header

        # turn the annotations into a table
        df = pd.DataFrame(annotations, columns=self._header)

        # apply any table level annotations, i.e. they need the
        # complete table to add extra info
        for module in self._table_annotation_modules:
            module.annotate(self._header, df)

        # prevent subsequent files for changing headers
        self._update_header = False

        self.write_metadata()

        #df.to_csv('test.tsv', sep='\t', header=True, index=False)
        return df

    def annotate(self, location, row_map: Mapping[str, Any]):
        row: list[Union[str, int, float]] = []

        # genes = self._refseq_annotation.annotate_location(location)
        # row_map['genes'] = genes

        # refseqs = []
        # symbols = []
        # entrezes = []
        # chrs = []
        # strands = []
        # starts = []
        # ends = []
        # types = []
        # tss_list = []

        # for gene in genes:
        #     gene_annotation = self._refseq_genes.get_gene(gene.id)

        #     chrs.append(gene_annotation.chr)
        #     starts.append(gene_annotation.start)
        #     strands.append(gene_annotation.strand)
        #     ends.append(gene_annotation.end)
        #     entrezes.append(gene_annotation.entrez)
        #     refseqs.append(gene_annotation.refseq)
        #     symbols.append(gene_annotation.symbol)

        #     # include tss distance for all peaks regardless
        #     tss_list.append(str(gene.d))

        #     types.append(','.join(reversed(sorted(gene.types))))

        # if len(entrezes) > 0:
        #     row.append(';'.join(refseqs))
        #     row.append(';'.join(entrezes))
        #     row.append(';'.join(symbols))
        #     row.append(';'.join(types))
        #     row.append(';'.join(tss_list))

        #     row_map[headings.GENE_SYMBOL] = ';'.join(symbols)
        # else:
        #     row.extend([text.NA] * 5)

        #     row_map[headings.GENE_SYMBOL] = text.NA

        #
        # Look at which gene TSS is closest
        #

        #sys.stderr.write(location.to_string() + '\n')

        #print(location.chr, location.start, 'loc', file=sys.stderr)

        #tss_gene = self._refseq_start.get_closest_gene(location)
        #gene_annotation = self._refseq_genes.get_gene(tss_gene.id)

        # Closest TSS gene symbol
        # row.append(gene_annotation.refseq)
        # row.append(gene_annotation.entrez)
        # row.append(gene_annotation.symbol)
        # row.append(','.join(reversed(sorted(tss_gene.types))))
        # row.append(str(tss_gene.d))
        #row_map[headings.CLOSEST_GENE_SYMBOL] = gene_annotation.symbol

        #
        # Look at which gene end is closest
        #

        # tss_gene = self._refseq_end.get_closest_gene(location)
        # gene_annotation = self._refseq_genes.get_gene(tss_gene.id)

        # # Closest TSS gene symbol
        # row.append(gene_annotation.refseq)

        # # Closest TSS gene entrez
        # row.append(gene_annotation.entrez)

        # # Closest TSS gene symbol
        # row.append(gene_annotation.symbol)

        # # Closest TSS gene types
        # row.append(','.join(reversed(sorted(tss_gene.types))))

        # # Closest TSS gene distances, since there can only be one closest distance
        # # all entries will have the same d so take the first. This is for
        # # odd entries where there might be two or more genes with the same starts
        # row.append(str(tss_gene.d))

        # #
        # # miRs
        # #

        # mir_annotations = collections.defaultdict(str)
        # mss = collections.defaultdict(str)
        # self.mir_annotation.annotate_location(
        #     location, self.mir_max_gap, mir_annotations, mss)

        # mirs = []
        # mir_types = []
        # mss_list = []

        # #
        # # If the peak is a promoter, check to see if there is a mir in the gene
        # #

        # for gene in genes:
        #     gene_annotation = self._refseq_genes.get_gene(gene.id)

        #     promoter = False

        #     for type in gene.types:
        #         if re.match(r'.*promoter.*', type):
        #             promoter = True
        #             break

        #     if promoter == False:
        #         continue

        #     transcript_mirs = self.mir_transcript_annotation.find_mirs(
        #         gene_annotation.entrez)

        #     if len(transcript_mirs) == 0:
        #         continue

        #     symbol = gene_annotation.symbol

        #     if symbol == text.NA:
        #         symbol = 'gene'

        #     for mir in sorted(transcript_mirs):
        #         mirs.append(mir)
        #         mir_types.append(symbol + '_promoter')
        #         mss_list.append(text.NA)

        # # Standard mir annotation

        # for mir in sorted(mir_annotations):
        #     mirs.append(mir)
        #     mir_types.append(mir_annotations[mir])
        #     mss_list.append(mss[mir])

        # if len(mirs) > 0:
        #     closest_mss = genomic.get_closest_tss(mss_list)

        #     row.append(';'.join(mirs))
        #     row.append(';'.join(mir_types))
        #     row.append(closest_mss)
        #     row.append(';'.join(mss_list))
        # else:
        #     row.append(text.NA)
        #     row.append(text.NA)
        #     row.append(text.NA)
        #     row.append(text.NA)

        #
        # Modules
        #

        for module in self._annotation_modules:
            #print("Running module", module.get_names(), file=sys.stderr)
            annotations = module.update_row(location, row_map)

            names = module.get_names()

            # update row map
            for i in range(0, len(names)):
                row_map[names[i]] = annotations[i]

            row.extend(annotations)

        #
        # Closest genes
        #

        # tss_genes = self._tss_refseq.get_n_closest_genes(location, n=50)

        # used = set()
        # n = 0

        # # find 5 unique entrez ids since there can be multiple refseqs with
        # # same gene id

        # for tss_gene in tss_genes:
        #     gene_annotation = self._refseq_genes.get_gene(tss_gene.id)

        #     if gene_annotation.entrez in used:
        #         continue

        #     ret.append(gene_annotation.refseq)
        #     ret.append(gene_annotation.entrez)
        #     ret.append(gene_annotation.symbol)
        #     ret.append(str(tss_gene.d))

        #     used.add(gene_annotation.entrez)

        #     n += 1

        #     if n == 5:
        #         break

        # ret.extend([text.NA] * 4 * (5 - n))

        # Finish the annotation
        return row
