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

import collections
import sys
from typing import Any, Mapping

from . import genomic
from . import text
from . import  headings

TADS_HEADING = 'TAD domains'
TAD_HEADING = 'Genes in same TAD domain'
IS_TAD_HEADING = 'Gene in TAD'
IS_CLOSEST_TAD_HEADING = 'Closest gene in TAD'

class TADAnnotation(genomic.Annotation):
    """
    Annotate peaks for all possible refseq genes they might overlap.
    """

    def __init__(self, file:str, bin_size: int = 100):
        """Create new tad annotation object.

        Args:
            file (str): file of TAD annotations.
            bin_size (int, optional): Search bin size. Defaults to 100.
        """        
        super().__init__('TAD')

        self._bin_size = bin_size

        self._tads = collections.defaultdict(list)
        self._tad_bins = collections.defaultdict(
            lambda: collections.defaultdict(int))

        print(f'Loading TADs from {file}...', file=sys.stderr)

        f = open(file, 'r')

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            tokens = line.split("\t")

            chr = tokens[0]
            start = int(tokens[1]) + 1
            end = int(tokens[2])
            location = genomic.Location(chr, start, end)

            ids = tokens[3]
            gene_symbols = tokens[4]

            start_bin = int(start / bin_size)
            end_bin = int(end / bin_size)

            for bin in range(start_bin, end_bin + 1):
                # apparently refseq genes from the ucsc always report
                # coordinates on the forward strand regardless of orientation
                # Each bin stores the index to a tad containing the genes
                # spanning this region
                self._tad_bins[chr][bin] = len(self._tads[chr])

                # if chr == 'chr1':
                #print(chr, location, bin, len(self._tads), file=sys.stderr)

            self._tads[chr].append([location, ids, gene_symbols])

            # if chr == 'chr1':
            #    print(self._tads[chr], file=sys.stderr)

        f.close()

        #print(self._tad_bins, file=sys.stderr)
        #print(self._tad_bins['chr10'], file=sys.stderr)

        # exit(0)

        print('Finished loading TADs.', file=sys.stderr)

    def get_names(self):
        return [TADS_HEADING, TAD_HEADING]

    def update_row(self, location:genomic.Location, row_map:Mapping[str, Any]):
        chr = location.chr

        if chr not in self._tad_bins:
            return text.NA

        start_bin = int(location.start / self._bin_size)
        end_bin = int(location.end / self._bin_size)

        # first find all the variants we might belong to
        genes = set()
        domains = []

        for bin in range(start_bin, end_bin + 1):
            if bin not in self._tad_bins[chr]:
                continue

            tad = self._tads[chr][self._tad_bins[chr][bin]]

            # if (chr == 'chr10'):
            #    print('lll', bin, location, tad[0], self._tad_bins[chr][bin], genomic.is_overlapping(location, tad[0]), file=sys.stderr)

            if genomic.is_overlapping(location, tad[0]):
                loc = str(tad[0])

                if loc not in domains:
                    domains.append(loc)

                # track all the unique genes we find
                genes.update(tad[2].split(';'))

        domains = ';'.join(domains)

        if domains == '':
            domains = text.NA

        genes = ';'.join(sorted(genes))

        if genes == '':
            genes = text.NA

        row_map[TADS_HEADING] = domains
        row_map[TAD_HEADING] = genes

        return [domains, genes]


class TADAnnotationFactory:
    _classifiers = collections.defaultdict(lambda: {})

    @staticmethod
    def getInstance(file: str, bin_size: int = 100):
        if bin_size not in TADAnnotationFactory._classifiers[file]:
            TADAnnotationFactory._classifiers[file][bin_size] = TADAnnotation(file, bin_size)

        return TADAnnotationFactory._classifiers[file][bin_size]

    def __init__(self):
        raise Exception("Call TADAnnotationFactory.getInstance()")


class IsTADAnnotation(genomic.Annotation):
    """
    Annotate whether the gene of a peak is listed in TAD domains the peak belongs to. 
    Mostly for Laura's benefit.
    """

    def __init__(self):
        super().__init__('is-tad')

    def get_names(self):
        return [IS_TAD_HEADING]

    def update_row(self, location:genomic.Location, row_map:Mapping[str, Any]):
        
        #print(row_map, file=sys.stderr)
        #print(row_map[headings.GENE_SYMBOL], file=sys.stderr)

        if TAD_HEADING not in row_map:
            return '0'
        
        if headings.GENE_SYMBOL not in row_map:
            return '0'

        gene = row_map[headings.GENE_SYMBOL]
        ret = 1 if gene in row_map[TAD_HEADING] else 0

        return [ret]


class IsClosestTADAnnotation(genomic.Annotation):
    """
    Annotate if the closest gene to a peak is in the TAD domain of this peak.
    """

    def __init__(self):
        super().__init__('is-closest-tad')

    def get_names(self):
        return [IS_CLOSEST_TAD_HEADING]

    def update_row(self, location:genomic.Location, row_map:Mapping[str, Any]):
        if TAD_HEADING not in row_map:
            return [0]
        
        if headings.CLOSEST_GENE_SYMBOL not in row_map:
            return [0]

        gene = row_map[headings.CLOSEST_GENE_SYMBOL]
        ret = 1 if gene in row_map[TAD_HEADING] else 0

        return [ret]
