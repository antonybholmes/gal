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
import re
from typing import Optional

from . import genomic
from . import headings
from . import text
from . import expression

GENE_ID = 'gene_id'
GENE_SYMBOL = 'gene_symbol'
GENE_REFSEQ = "refseq"
GENE_ENTREZ = 'entrez'
GENE_ENSEMBL = 'ensembl'


def parse_rdf_gene_id(text):
    """
    An RDF id contains a core gene id plus a decimal to indicate the variant
    """

    return re.match(r'(RDF\d+).*', text).group(1)


def parse_rdf_gene_variant_id(text):
    return int(re.match(r'.*(\d+)$', text).group(1))


def create_variant_id(id: str, chr: str, start: int, end: int):
    return f'{id}#{genomic.location_string(chr, start, end)}'


def parse_location_from_variant(location: str) -> genomic.Location:
    matcher = re.match(r'.+(chr.+):(\d+)-(\d+)', location)

    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))

    location = genomic.Location(chr, start, end)

    return location


def parse_id_from_variant(variant_id: str) -> str:
    matcher = re.match(r'^([^#]+).*', variant_id)

    return matcher.group(1)


def find_best_p_value(header, tokens):
    # Required for overlap files where both p-values from each
    # replicate are maintained. We need to select one.

    min_p = 1

    indices = text.find_indices(header, headings.P_VALUE)

    for i in indices:
        if tokens[i] != text.NA:
            p = float(tokens[i])

            if p < min_p:
                min_p = p

    return min_p


def find_best_score(header, tokens):
    # Required for overlap files where both scores from each
    # replicate are maintained. We need to select one.

    max_score = 0

    indices = text.find_indices(header, headings.SCORE)

    for i in indices:
        if tokens[i] != text.NA:
            s = float(tokens[i])

            if s > max_score:
                max_score = s

    return max_score


class Gene(genomic.Feature):
    """
    Represents a named genomic location with multiple ids.
    """

    def __init__(self, id: str, symbol: str, strand: str, chr: str, start: int, end: int):
        """
        Creates a gene genomic feature

        Args:
            id (str): gene id, such as an refseq id
            symbol (str): gene name, such as "BCL6"
            strand (str): '+' or '-'
            chr (str): chromosome
            start (int): 1 based start location
            end (int): 1 based end location
        """
        super().__init__(chr, start, end, strand)
        self.add_id(GENE_ID, id)
        self.add_id(GENE_SYMBOL, symbol)

    @property
    def id(self) -> str:
        return self.get_id(GENE_ID)

    @property
    def symbol(self) -> str:
        return self.get_id(GENE_SYMBOL)

    @property
    def name(self) -> str:
        return self.symbol


class RefSeqGene(Gene):
    def __init__(self, refseq: str, entrez: str, symbol: str, strand: str, chr: str, start: int, end: int):
        super().__init__(refseq, symbol, strand, chr, start, end)
        self.add_id(GENE_REFSEQ, refseq)
        self.add_id(GENE_ENTREZ, entrez)

    @property
    def refseq(self) -> str:
        return self.get_id(GENE_REFSEQ)

    @property
    def entrez(self) -> str:
        return self.get_id(GENE_ENTREZ)


class RefSeqGenes:
    """
    Gene lookup by variant id, other dbs store variant ids
    # instead of objects. This is to prevent redundancy
    """

    def __init__(self, file):
        self._genes = {}  # collections.defaultdict(Gene)
        self._file = file

        print(f'Loading genes from {file}...', file=sys.stderr)

        # annotation.REFSEQ_FILE
        f = open(file, 'r')

        # skip header
        f.readline()

        # To account for multiple versions of a gene, allocate each entrez
        # id a unique index

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')

            refseq = tokens[0]
            entrez = tokens[1]
            symbol = tokens[2]
            chr = tokens[3]
            strand = tokens[4]
            # ucsc convention
            start = int(tokens[5]) + 1
            end = int(tokens[6])

            if re.match(r'.*n/a.*', entrez):
                continue

            if re.match(r'.*n/a.*', symbol):
                continue

            if re.match(r'.*MIR.*', symbol):
                continue

            # index gene on an id and coordinates to keep it unique

            id = create_variant_id(refseq, chr, start, end)

            self._genes[id] = RefSeqGene(
                refseq, entrez, symbol, strand, chr, start, end)

        f.close()

    @property
    def genes(self) -> dict[str, RefSeqGene]:
        return self._genes

    @property
    def file(self) -> str:
        return self._file

    def get_gene(self, variant_id: str) -> Optional[RefSeqGene]:
        return self._genes.get(variant_id, None)


class GeneOrientatedPeaks:
    """
    Core annotation for gene oriented peaks
    """

    def __init__(self, type: str = "Peak", all_gene_mode: bool = False, show_scores: bool = True, collapse_id: str = 'gene_symbol'):
        """Create a new object.

        Args:
            type (str, optional): Indicate what is being annotated. Defaults to "Peak".
            all_gene_mode (bool, optional): If true, will use closest gene if peak is
            not overlapping a gene. Defaults to False.
        """

        self.type = type
        self._all_gene_mode = all_gene_mode
        self._show_scores = show_scores
        self._collapse_id = collapse_id
        self.expression_list: list[expression.GeneExpression] = []
        self.expression_list_headers: list[str] = []

        # self.affy_gene_cb_vs_m_expression = expression.AffyGeneCBvsMExpression()
        # self.affy_gene_cb_vs_n_expression = expression.AffyGeneCBvsNExpression()
        # self.rna_gene_cb_vs_m_expression = expression.RnaSeqGeneCBvsMExpression()
        # self.rna_gene_cb_vs_n_expression = expression.RnaSeqGeneCBvsNExpression()

        # self.expression_list_headers.append('GEP Affy GCvsN')
        # self.expression_list_headers.append('GEP Affy GCvsM')
        # self.expression_list.append(
        #    expression.AffyGeneCBvsNExpression())
        # self.expression_list.append(
        #    expression.AffyGeneCBvsMExpression())

        # self.expression_list_headers.append('RNA-seq GRCh38 GCvsN')
        # self.expression_list_headers.append('RNA-seq GRCh38 GCvsM')
        # self.expression_list.append(
        #    expression.RnaSeqGeneCBvsNExpression())
        # self.expression_list.append(
        #    expression.RnaSeqGeneCBvsMExpression())

        # RNA-seq from deseq2
        # self.expression_list_headers.append('RNA-seq GCvsN DeSeq2')
        # self.expression_list_headers.append('RNA-seq GCvsM DeSeq2')
        # self.expression_list.append(expression.RnaSeqGeneCBvsNDeSeq2Expression())
        # self.expression_list.append(expression.RnaSeqGeneCBvsMDeSeq2Expression())

        #
        # SUD10 stuff
        #
        # self.expression_list_headers.append('RNA-seq SUD10 WTvsD83V')
        # self.expression_list_headers.append('RNA-seq SUD10 WTvsSTOP')
        # self.expression_list_headers.append('RNA-seq SUD10 D83VvsSTOP')
        # self.expression_list.append(
        #    expression.RnaSeqGeneSUD10WTvsD83vExpression())
        # self.expression_list.append(
        #    expression.RnaSeqGeneSUD10WTvsStopExpression())
        # self.expression_list.append(
        #    expression.RnaSeqGeneSUD10D83vsStopExpression())

        #
        # Mouse MEF2B vs WT
        #

        # self.expression_list_headers.append('GEP Affy MEF2B Mouse WTvsKO')
        # self.expression_list.append(
        #    expression.GEPMEF2BMouseWTvsKOExpression())

        #
        # miR annotations
        #

        # self.solid_mir_expression = expression.SolidMirExpression()
        # self.solid_mir_cb_vs_n_expression = expression.SolidMirCBvsNExpression()
        # self.agilent_mir_cb_vs_m_expression = expression.AgilentMirCBvsMExpression()
        # self.agilent_mir_cb_vs_n_expression = expression.AgilentMirCBvsNExpression()

        #
        # Sets for annotations
        #

        self._collapsed_refseqs = collections.defaultdict(set)
        self._collapsed_entrezes = collections.defaultdict(set)
        self._collapsed_gene_symbols = collections.defaultdict(set)
        # self._collapsed_p = collections.defaultdict(list)
        # self._collapsed_scores = collections.defaultdict(list)
        self._collapsed_locations = collections.defaultdict(list)
        self._collapsed_tss = collections.defaultdict(list)
        self._collapsed_types = collections.defaultdict(list)
        self._collapsed_centromeres = collections.defaultdict(list)

        self._collapsed_scores = collections.defaultdict(
            lambda: collections.defaultdict(list))
        self._max_score_n = 0
        self._total_score_n = 0

        self.mirs = set()
        self.refseqs = set()

        # self._collapsed_mir_types = collections.defaultdict(list)
        # self._collapsed_mirs_p = collections.defaultdict(list)
        # self._collapsed_mirs_scores = collections.defaultdict(list)
        # self._collapsed_mirs_locations = collections.defaultdict(list)
        # self._collapsed_mss = collections.defaultdict(list)

    def add_expression(self, name: str, expression: expression.GeneExpression):
        """Adds a new expression object to the pipeline which will add
        an annotation column to the output table.

        Args:
            name (str): Name of annotation.
            expression (expression.GeneExpression): GeneExpression object that
            will be added into the workflow.
        """
        self.expression_list_headers.append(name)
        self.expression_list.append(expression)

    def load_peaks(self, file: str):
        f = open(file, 'r')

        # skip header
        header = f.readline().strip().split('\t')

        location_column = text.find_index(
            header, headings.LOCATION)
        entrez_column = text.find_index(
            header, headings.ENTREZ_ID)
        refseq_column = text.find_index(
            header, headings.REFSEQ_ID)
        symbol_column = text.find_index(
            header, headings.GENE_SYMBOL)
        type_column = text.find_index(header, 'Relative To Gene')

        # p_column = text.find_index(header, headings.P_VALUE)
        # score_column = text.find_index(header, headings.SCORE)

        tss_column = text.find_index(
            header, headings.TSS_DISTANCE)
        centromere_column = text.find_index(
            header, headings.CENTROMERE)

        closest_entrez_column = text.find_index(
            header, headings.CLOSEST_ENTREZ_ID)
        closest_refseq_column = text.find_index(
            header, headings.CLOSEST_REFSEQ_ID)
        closest_symbol_column = text.find_index(
            header, headings.CLOSEST_GENE_SYMBOL)
        closest_type_column = text.find_index(
            header, headings.CLOSEST_RELATIVE)
        closest_tss_column = text.find_index(
            header, headings.CLOSEST_TSS_DISTANCE)

        total_score_column = text.find_index(
            header, "Total Score")

        max_score_column = text.find_index(
            header, "Max Score")

        # mir_column = text.find_index(header, headings.MIR_SYMBOL)
        # mir_type_column = text.find_index(header, 'Relative To miR')
        # mss_column = text.find_index(header, 'mIR Start Distance')

        lc = 0

        for line in f:
            ls = line.rstrip()

            lc += 1

            if len(ls) == 0:
                continue

            tokens = ls.split('\t')

            # total_score = genes.find_best_p_value(
            #    header, tokens)  # float(tokens[p_column])
            # score = genes.find_best_score(
            #    header, tokens)  # = float(tokens[score_column])
            location = tokens[location_column]

            # Preference is given to peaks within a gene
            refseqs = tokens[refseq_column].split(';')
            entrezes = tokens[entrez_column].split(';')
            symbols = tokens[symbol_column].split(';')

            tsses = tokens[tss_column].split(';')
            types = tokens[type_column].split(';')
            centromere = tokens[centromere_column]

            total_score = float(
                tokens[total_score_column]) if total_score_column != -1 else -1
            max_score = float(tokens[max_score_column]) if max_score_column != -1 else -1

            for i in range(0, len(refseqs)):
                # The core id identifies the unique genes on a peak, rather than
                #
                id_map = {'gene_symbol':symbols[i], 'refseq':refseqs[i], 'entrez':entrezes[i]}

                id = id_map[self._collapse_id]

                tss = tsses[i]
                type = types[i]

                if 'SDCBP2' in id:
                    print('in', symbol_column, id, tss, type, file=sys.stderr)

                if id != text.NA:
                    self._collapsed_refseqs[id].add(id_map['refseq'])
                    self._collapsed_entrezes[id].add(id_map['entrez'])
                    self._collapsed_gene_symbols[id].add(id_map['gene_symbol'])
                    # self._collapsed_p[refseq].append(p)
                    # self._collapsed_scores[refseq].append(score)
                    self._collapsed_locations[id].append(location)
                    self._collapsed_tss[id].append(tss)
                    self._collapsed_types[id].append(type)
                    self._collapsed_centromeres[id].append(centromere)

                    # for all peaks we find keep track of every score
                    if max_score != "":
                        self._collapsed_scores[id]['max_score'].append(
                            max_score)

                    if total_score != "":
                        self._collapsed_scores[id]['total_score'].append(
                            total_score)

            if self._all_gene_mode: #and tokens[refseq_column] == text.NA:
                entrezes = tokens[closest_entrez_column].split(';')
                symbols = tokens[closest_symbol_column].split(';')
                refseqs = tokens[closest_refseq_column].split(';')
                tsses = tokens[closest_tss_column].split(';')
                types = tokens[closest_type_column].split(';')

                for i in range(0, len(refseqs)):
                    # The core id identifies the unique genes on a peak, rather than
                    #
                    id_map = {'gene_symbol':symbols[i], 'refseq':refseqs[i], 'entrez':entrezes[i]}

                    id = id_map[self._collapse_id]

                    tss = tsses[i]
                    type = types[i]

                    if id == 'ANO10':
                        print(id, tsses, types, file=sys.stderr)

                    if id != text.NA:
                        self._collapsed_refseqs[id].add(id_map['refseq'])
                        self._collapsed_entrezes[id].add(id_map['entrez'])
                        self._collapsed_gene_symbols[id].add(id_map['gene_symbol'])
                        # self._collapsed_p[refseq].append(p)
                        # self._collapsed_scores[refseq].append(score)
                        self._collapsed_locations[id].append(location)
                        self._collapsed_tss[id].append(tss)
                        self._collapsed_types[id].append(type)
                        self._collapsed_centromeres[id].append(centromere)

                        if max_score != "":
                            self._collapsed_scores[id]['max_score'].append(
                                max_score)

                        if total_score != "":
                            self._collapsed_scores[id]['total_score'].append(
                                total_score)

            # mir = tokens[mir_column]

            # if mir != text.NA:
            #     self.mirs.add(mir)
            #     self._collapsed_types[mir].append(tokens[mir_type_column])
            #     #self._collapsed_p[mir].append(p)
            #     #self._collapsed_scores[mir].append(score)
            #     self._collapsed_locations[mir].append(location)
            #     self._collapsed_tss[mir].append(tokens[mss_column])
            #     self._collapsed_centromeres[mir].append(centromere)

        f.close()

        print(lc, file=sys.stderr)

    def print_header(self):
        """Prints the table header to std out.
        """
        print('\t'.join(self.get_header()))

    def get_header(self) -> list[str]:
        """Returns the table header.

        Returns:
            list[str]: The table header as a list of string.

        """
        ret = []
        ret.append(headings.REFSEQ_ID)
        ret.append(headings.ENTREZ_ID)
        ret.append(headings.GENE_SYMBOL)
        # ret.append(f'\tGEP Affy GCvsN")
        # ret.append(f'\tGEP Affy GCvsM")
        # ret.append(f'\tRNAseq GCvsN")
        # ret.append(f'\tRNAseq GCvsM")

        for header in self.expression_list_headers:
            ret.append(header)

        ret.append(f'{self.type} Relative To Gene')
        ret.append(f'{self.type} {headings.TSS_DISTANCE}')
        ret.append(f'{self.type} TSS Closest Distance')
        ret.append(f'{self.type} {headings.CENTROMERE}')
        # ret.append('miR Symbol')
        # ret.append('miREP Agilent GCvsN')
        # ret.append('miREP Agilent GCvsM')
        # ret.append('sRE GCvsNM')
        # ret.append('sRE GCvsN')
        # ret.append(f'{self.type} Relative To miR')
        # ret.append(f'{self.type} miR Start Closest Distance')
        # ret.append(f'{self.type} miR Start Distance')
        # ret.append('Best P-value (ChIPseeqer)')
        # ret.append('Best Score (ChIPseeqer)')
        ret.append(f'{self.type} Count')
        ret.append(f'{self.type} Genomic Locations')

        if self._show_scores:
            ret.append(f'Highest Total x Max Score')
            ret.append(f'Highest Total Score')
            ret.append(f'Highest Max Score')

            self._total_score_n = 0

            for refseq in self._collapsed_scores:
                self._total_score_n = max(self._total_score_n, len(
                    self._collapsed_scores[refseq]['total_score']))

            if self._total_score_n > 0:
                ret.extend(
                    [f'Total Score {i + 1}' for i in range(0, self._total_score_n)])

            self._max_score_n = 0

            for refseq in self._collapsed_scores:
                self._max_score_n = max(self._max_score_n, len(
                    self._collapsed_scores[refseq]['max_score']))

            if self._max_score_n > 0:
                ret.extend(
                    [f'Max Score {i + 1}' for i in range(0, self._max_score_n)])

        return ret

    def get_ids(self) -> list[str]:
        """Returns the ids of the genes/gene ids to write

        Returns:
            list[str]: _description_
        """
        return list(sorted(self._collapsed_refseqs))

    def get_mirs(self):
        return sorted(self.mirs)

    def gene_orient_peak(self, id):
        refseqs = list(sorted(self._collapsed_refseqs[id]))
        entrezes = list(sorted(self._collapsed_entrezes[id]))
        symbols = list(sorted(self._collapsed_gene_symbols[id]))

        ret = []

        ret.append(';'.join(refseqs))
        ret.append(';'.join(entrezes))
        ret.append(';'.join(symbols))

        # ret.append(f'\t{self.affy_gene_cb_vs_n_expression.get_expression(entrez))
        # ret.append(f'\t{self.affy_gene_cb_vs_m_expression.get_expression(entrez))
        # ret.append(f'\t{self.rna_gene_cb_vs_n_expression.get_expression(entrez))
        # ret.append(f'\t{self.rna_gene_cb_vs_m_expression.get_expression(entrez))

        for expression in self.expression_list:
            ret.append(expression.get_expression(refseqs + symbols))

        ret.append(';'.join(self._collapsed_types[id]))

        ret.append(';'.join(self._collapsed_tss[id]))

        # if there are some nearest tss, print the closest
        ret.append(
            genomic.get_closest_tss(self._collapsed_tss[id]))

        # Centromeres
        ret.append(';'.join(self._collapsed_centromeres[id]))

        # no mir symbol
        # ret.append(text.NA)

        # no agilent expression
        # ret.append(text.NA)
        # ret.append(text.NA)

        # no small rna
        # ret.append(text.NA)
        # ret.append(text.NA)

        # no peak relative to mir
        # ret.append(text.NA)
        # ret.append(text.NA)
        # ret.append(text.NA)

        # pick the smallest p
        # p = sorted(self._collapsed_p[id])
        # ret.append(p[0])

        # pick the largest score
        # scores = sorted(self._collapsed_scores[id], reverse=True)
        # ret.append(scores[0])

        # peak count
        ret.append(len(self._collapsed_locations[id]))

        ret.append(';'.join(self._collapsed_locations[id]))

        highest_total_max = -1
        highest_total = -1
        highest_max = -1

        if len(self._collapsed_scores[id]['total_score']) > 0:
            for i in range(len(self._collapsed_scores[id]['total_score'])):
                s = self._collapsed_scores[id]['total_score'][i] * \
                    self._collapsed_scores[id]['max_score'][i]

                if s > highest_total_max:
                    highest_total_max = s
                    highest_total = self._collapsed_scores[id]['total_score'][i]
                    highest_max = self._collapsed_scores[id]['max_score'][i]

        if self._show_scores:
            ret.append(str(highest_total_max))
            ret.append(str(highest_total))
            ret.append(str(highest_max))

            # only add scores to table if there are
            if len(self._collapsed_scores[id]['total_score']) > 0:
                d = [text.NA] * self._total_score_n
                d[0:len(self._collapsed_scores[id]['total_score'])] = [str(x)
                                                                      for x in self._collapsed_scores[id]['total_score']]
                ret.extend(d)

            if len(self._collapsed_scores[id]['max_score']) > 0:
                d = [text.NA] * self._max_score_n
                d[0:len(self._collapsed_scores[id]['max_score'])] = [str(x)
                                                                    for x in self._collapsed_scores[id]['max_score']]
                ret.extend(d)

        return ret
    
    def print_gene_orient_peak(self, id):
        print('\t'.join([str(x) for x in self.gene_orient_peak(id)]))

    def print_gene_orient_peaks(self):
        self.print_header()

        for id in self.get_ids():
            self.print_gene_orient_peak(id)


    def mir_orient_peak(self, mir):
        ret = []

        ret.append(text.NA)

        ret.extend([text.NA] * len(self.expression_list_headers))

        # Fill in the gap
        ret.extend([text.NA] * 5)

        ret.append(';'.join(self._collapsed_centromeres[mir]))
        ret.append(mir)
        ret.append(
            self.agilent_mir_cb_vs_n_expression.get_expression(mir))
        ret.append(
            self.agilent_mir_cb_vs_m_expression.get_expression(mir))
        ret.append(self.solid_mir_expression.get_expression(mir))
        ret.append(self.solid_mir_cb_vs_n_expression.get_expression(mir))
        ret.append(';'.join(self._collapsed_types[mir]))
        ret.append(
            genomic.get_closest_tss(self._collapsed_tss[mir]))
        ret.append(';'.join(self._collapsed_tss[mir]))

        # pick the smallest p
        p = sorted(self._collapsed_p[mir])
        ret.append(p[0])

        # pick the largest score
        scores = sorted(self._collapsed_scores[mir], reverse=True)
        ret.append(scores[0])

        ret.append(len(self._collapsed_locations[mir]))
        ret.append(';'.join(self._collapsed_locations[mir]))

        print('\t'.join([str(x) for x in ret]))

        return ret

    def print_log(self):
        """
        Print a log to indicate what was used for annotations.
        """

        f = open('genes.log', 'w')

        f.write('Expression Type\tSource\n')

        for i in range(0, len(self.expression_list_headers)):
            f.write(
                f'{self.expression_list_headers[i]}\t{self.expression_list[i].get_file()}\n')

        f.close()
