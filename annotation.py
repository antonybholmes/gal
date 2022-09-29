# -*- coding: utf-8 -*-
"""
Annotation functions

Created on Thu Jun 26 09:43:29 2014

@author: Antony Holmes
"""
import sys
import collections
import re

from . import text

#import radix_tree


#REFSEQ_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_refseq_exons_entrez_hg19_20150217.txt"
#REFSEQ_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_refseq_exons_entrez_hg19.txt"
#REFSEQ_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_refseq_exons_entrez_canonical_only_hg19.txt"

#ENSEMBL_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_ensembl_exons_gene_symbol_hg19.txt"

#RDF_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/rdf_ucsc_refseq_ensembl_genes_hg19.txt"

# refseqs only
#RDF_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/rdf_ucsc_refseq_genes_hg19.txt"


def get_mir_id(id):
    mir = id.lower()

    #  if not re.match('(hsa)-(.+?)-(.+?).*', s):
    #    return s
    #
    #  matcher = re.match('(hsa)-([^\-]+)-([^\-]+)', s)
    #
    #  mir = matcher.group(1) + "-" + matcher.group(2) + "-" + matcher.group(3)

    return mir


def get_mir_precursor_id(id):
    mir = id.lower()

    if not re.match('(hsa)-(.+?)-(.+?).*', mir):
        return mir

    matcher = re.match('(hsa)-([^\-]+)-([^\-]+)', mir)
    mir = matcher.group(1) + "-" + matcher.group(2) + "-" + matcher.group(3)

    return mir


def gene_orient_peak_locations(file,
                               entrez_gene_map,
                               collapsed_genes,
                               collapsed_locations,
                               collapsed_p):
    f = open(file, 'r')

    # skip header
    header = f.readline().split('\t')

    for line in f:
        ls = line.strip()

        if len(ls) == 0:
            continue

        tokens = line.split('\t')

        id = tokens[0]

        location = tokens[1]
        p = float(tokens[3])

        entrezes = tokens[text.find_index(header, 'entrez')].split(';')

        for entrez in entrezes:
            if entrez == text.NA:
                continue

            if entrez in entrez_gene_map:
                collapsed_genes[entrez] = entrez_gene_map[entrez]
            else:
                collapsed_genes[entrez] = text.NA

            collapsed_locations[entrez].append(location)

            if entrez in collapsed_p:
                if p < collapsed_p[entrez]:
                    collapsed_p[entrez] = p
            else:
                collapsed_p[entrez] = p

    f.close()


def unique_symbols(line):
    tokens = line.split(',')

    s = set()

    for item in tokens:
        s.add(item)

    return ','.join(sorted(s))


class EntrezGeneLookup(object):
    """
    Loads a mapping between entrez ids and gene symbols.
    """

    def __init__(self, file):
        self.entrez_gene_map = collections.defaultdict(str)

        print(f'Loading genes from {file}...', file=sys.stderr)

        f = open(file, 'r')

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')

            id = tokens[0]
            symbol = unique_symbols(tokens[1])

            if symbol == text.NA:
                continue

            self.entrez_gene_map[id] = symbol

        f.close()

        print('Finished.', file=sys.stderr)

    def get_symbol(self, entrez):
        if entrez not in self.entrez_gene_map:
            return text.NA

        return self.entrez_gene_map[entrez]


class MirPreviousIds(object):
    def __init__(self):
        self.previous_ids = collections.defaultdict(
            lambda: collections.defaultdict(bool))

        f = open(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirna_mature.txt", "r")

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')

            mir = tokens[1]

            if not re.match(r'hsa.*', mir):
                continue

            previous = tokens[2]

            if len(previous) < 1:
                continue

            old_ids = previous.split(";")

            for id in old_ids:
                #sys.stderr.write(get_mir_id(mir) + " " + get_mir_id(id) + "\n")
                self.previous_ids[get_mir_id(mir)][get_mir_id(id)] = True

        f.close()

    def get_previous_ids(self, id):
        id = get_mir_id(id)

        ret = []

        if id in self.previous_ids:
            for old_id in self.previous_ids[id]:
                ret.append(old_id)

        return ret



