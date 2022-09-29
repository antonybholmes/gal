from abc import ABC, abstractmethod
import sys
import collections
import re
import os.path

from . import annotation
from . import text
from . import headings


def gene_exp(genes_file, up_gene_file, down_gene_file):
    """
    Takes three gene lists for mapping genes to whether they are up
    or down regulated in a particular gene expression type e.g.
    """

    gene_exp = collections.defaultdict(str)

    load_gene_exp(genes_file, gene_exp, "not_moving")
    load_gene_exp(up_gene_file, gene_exp, "up")
    load_gene_exp(down_gene_file, gene_exp, "down")


def load_gene_exp(file, genes_map, type):
    if not os.path.exists(file):
        return

    f = open(file, "r")

    for line in f:
        ls = line.strip()

        if len(ls) == 0:
            continue

        tokens = ls.split("\t")

        if len(tokens[0]) == 0 or re.match(r'n\/a', tokens[0]):
            continue

        genes_map[tokens[0].lower()] = type

    f.close()


def load_classification(file, genes_map):
    """
    Load a classification file. The first row is a header, the
    first column is an entrez id, the second a gene symbol and
    the third whether the gene moves up or down etc.
    """
    if not os.path.exists(file):
        return

    sys.stderr.write("Opening classification file " + file + "...\n")

    f = open(file, "r")

    # skip header
    f.readline()

    for line in f:
        tokens = line.strip().split("\t")

        #
        # Entrez
        #
        ids = tokens[0].lower().split(";")

        for id in ids:
            if "n/a" not in id:
                # Map an entrez id to a gene symbol
                genes_map[id] = tokens[2]

        #
        # Symbol
        #
        ids = tokens[1].lower().split(";")

        for id in ids:
            if "n/a" not in id:
                # Map an entrez id to a gene symbol
                genes_map[id] = tokens[2]

    f.close()


def load_mir_exp(file, genes_map, type):

    f = open(file, "r")

    for line in f:
        ls = line.strip()

        if len(ls) == 0:
            continue

        tokens = ls.split("\t")

        if len(tokens[0]) == 0 or re.match(r'n\/a', tokens[0]):
            continue

        genes_map[annotation.get_mir_id(tokens[0])] = type

    f.close()


class PeakExpression(ABC):
    """ 
    Adds expression values to peak tables by looking up gene symbols
    and matching in RNA-seq data etc.
    """
    
    @abstractmethod
    def annotate(self, file:str):
        """
        Adds expression values to peak tables.

        Args:
            file (str): Peak table file to add expression columns to.
        """        
        ...


class Expression:
    """
    Maps symbols/ids to their expression type e.g. up or down regulated
    """

    def __init__(self, dir):
        self.dir = dir

        self.expression_map = collections.defaultdict(str)

        load_gene_exp(dir + "not_moving_genes.txt",
                      self.expression_map,
                      "not_moving")

        # now mark as down-regulated
        load_gene_exp(dir + "down_genes.txt",
                      self.expression_map,
                      "down")

        # now mark as up-regulated
        load_gene_exp(dir + "up_genes.txt",
                      self.expression_map,
                      "up")

        # Anything ambigous is loaded last (e.g. up/down) to override
        # any other choices
        load_gene_exp(dir + "ambiguous_genes.txt",
                      self.expression_map,
                      "ambiguous")

    def get_expression(self, id):
        s = id.lower()

        #sys.stderr.write("get exp for " + s + "\n")

        if s not in self.expression_map:
            return "n/a"

        return self.expression_map[s]

    def get_file(self):
        return self.dir


class RnaSeqExpression(Expression):
    def __init__(self, dir):
        super().__init__(dir)

        load_gene_exp(dir + "not_expressed_genes.txt",
                      self.expression_map,
                      "no_exp")


class GeneExpression:
    """
    Maps symbols/ids to their expression type e.g. up or down regulated
    using a table definition file. The file must contain a header and
    three columns corresponding to an entrez id, a gene symbol and a
    classification
    """

    def __init__(self, dir):
        self.dir = dir
        self.file = dir + "classification.txt"
        self.expression_map = collections.defaultdict(str)
        load_classification(self.file, self.expression_map)

    def get_expression(self, id):
        s = id.lower()

        if s in self.expression_map:
            return self.expression_map[s]
        else:
            return "n/a"

    def get_file(self):
        return self.dir


class MirExpression:
    def __init__(self):
        self.expression_map = collections.defaultdict(str)
        self.previous_ids = annotation.MirPreviousIds()

    def get_expression(self, id):
        s = annotation.get_mir_id(id)

        if s in self.expression_map:
            return self.expression_map[s]

        previous_ids = self.previous_ids.get_previous_ids(id)

        for pid in previous_ids:
            if pid in self.expression_map:
                return self.expression_map[pid]

        return "n/a"


#
# Microarray
#

class AffyGeneCTRLvsSIExpression(Expression):
    """
    Annotate whether a gene is up or down regulated in affymetrix data
    """

    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/bcl6_silencing/")


class AffyGeneCBvsMExpression(GeneExpression):
    """
    Annotate whether a gene is up or down regulated in affymetrix data
    """

    def __init__(self):
        #super(AffyGeneCBvsMExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references_cb_vs_m/")
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/cb_vs_m_20160608/")


class AffyGeneCBvsNExpression(GeneExpression):
    """
    Annotate whether a gene is up or down regulated in affymetrix data
    """

    def __init__(self):
        #super(AffyGeneCBvsNExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references_cb_vs_n/")
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/cb_vs_n_20160608/")


class AffyGeneDZvsLZExpression(GeneExpression):
    """
    Annotate whether a gene is up or down regulated in affymetrix data
    """

    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/dz_vs_lz_20160720/")

#
# RNA-seq classes
#


class RnaSeqGeneSUD10D83vsStopExpression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_sud10_d83v_vs_stop_20160610/")


class RnaSeqGeneSUD10WTvsD83vExpression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_sud10_wt_vs_d83v_20160610/")


class RnaSeqGeneSUD10WTvsStopExpression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_sud10_wt_vs_stop_20160610/")


class RnaSeqGeneSUD10D83vsStopDeSeq2Expression(GeneExpression):
    """
    Maps entrez ids to gene expression up or down. This version allows
    ambiguous classifications such as up;down or up;not_moving
    """

    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_sud10_d83v_vs_stop_deseq2_20160411/")


class RnaSeqGeneSUD10WTvsD83vDeSeq2Expression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_sud10_wt_vs_d83v_deseq2_20160411/")


class RnaSeqGeneSUD10WTvsStopDeSeq2Expression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_sud10_wt_vs_stop_deseq2_20160411/")


class RnaSeqGeneCBvsMExpression(GeneExpression):
    def __init__(self):
        super().__init__(
            "/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_cb_vs_m_20160610/")


class RnaSeqGeneCBvsNExpression(GeneExpression):
    def __init__(self):
        super().__init__(
            "/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_cb_vs_n_20160610/")


class RnaSeqGeneCBvsMDeSeq2Expression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_cb_vs_m_deseq2_20160411/")


class RnaSeqGeneCBvsNDeSeq2Expression(GeneExpression):
    def __init__(self):
        super().__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/rna_seq_cb_vs_n_deseq2_20160411/")

# 20171001

# class RnaSeqTPMCBvsNExpression(GeneExpression):
    # def __init__(self):
        #super(RnaSeqTPMCBvsNExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/tpm_cb_vs_n_20171001/")

# class RnaSeqTPMCBvsMExpression(GeneExpression):
    # def __init__(self):
        #super(RnaSeqTPMCBvsMExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/tpm_cb_vs_m_20171001/")


# 20171121

class RnaSeqTPMCBvsNExpression(GeneExpression):
    def __init__(self):
        super().__init__(
            "/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/tpm_cb_vs_n/")


class RnaSeqTPMCBvsMExpression(GeneExpression):
    def __init__(self):
        super().__init__(
            "/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/tpm_cb_vs_m/")


#
# GEP
#

# Old possibly wrong version
# class GEPMEF2BMouseWTvsKOExpression(Expression):
#  def __init__(self):
#    super(GEPMEF2BMouseWTvsKOExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references_mouse_mef2b_wt_vs_ko_20160411/")


# class GEPMEF2BMouseWTvsKOExpression(Expression):
    # """
    # Version 2 using standard max standard deviation collapsed table
    # """
    # def __init__(self):
        #super(GEPMEF2BMouseWTvsKOExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references_mouse_mef2b_wt_vs_ko_20160603/")

class GEPMEF2BMouseWTvsKOExpression(GeneExpression):
    """
    Maps entrez ids to gene expression up or down. This version allows
    ambiguous classifications such as up;down or up;not_moving
    """

    def __init__(self):
        super().__init__(
            "/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references_mouse_mef2b_wt_vs_ko_20160606/")


class DzLzGeneExpression(Expression):
    """
    Annotate whether a gene is up or down regulated in affymetrix data
    """

    def __init__(self):
        super().__init__(
            "/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/dark_zone_vs_light_zone/")


class SolidMirExpression(MirExpression):
    """
    Annotate whether a miR is up or down regulated in our SOLiD miR data
    """

    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/not_moving_mirs.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/down_mirs.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/up_mirs.txt",
                     self.expression_map,
                     "up")


class SolidMirLog2Expression(MirExpression):
    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/not_moving_mirs_log2.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/down_mirs_log2.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/up_mirs_log2.txt",
                     self.expression_map,
                     "up")


class SolidMirCBvsNExpression(MirExpression):
    """
    Annotate whether a miR is up or down regulated in our SOLiD miR data
    """

    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/not_moving_mirs_cb_vs_n.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/down_mirs_cb_vs_n.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/up_mirs_cb_vs_n.txt",
                     self.expression_map,
                     "up")


class SolidMirCBvsNLog2Expression(MirExpression):
    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/not_moving_mirs_cb_vs_n_log2.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/down_mirs_cb_vs_n_log2.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/solid_mir/up_mirs_cb_vs_n_log2.txt",
                     self.expression_map,
                     "up")


class AgilentMirCBvsMExpression(MirExpression):
    """
    Annotate whether a miR is up or down regulated in our SOLiD miR data
    """

    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/not_moving_mirs_cb_vs_m.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/down_mirs_cb_vs_m.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/up_mirs_cb_vs_m.txt",
                     self.expression_map,
                     "up")


class AgilentMirCBvsNExpression(MirExpression):
    """
    Annotate whether a miR is up or down regulated in our SOLiD miR data
    """

    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/not_moving_mirs_cb_vs_n.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/down_mirs_cb_vs_n.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/up_mirs_cb_vs_n.txt",
                     self.expression_map,
                     "up")


class AgilentMirCBvsMLog2Expression(MirExpression):
    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/not_moving_mirs_cb_vs_m_log2.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/down_mirs_cb_vs_m_log2.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/up_mirs_cb_vs_m_log2.txt",
                     self.expression_map,
                     "up")


class AgilentMirCBvsNLog2Expression(MirExpression):
    def __init__(self):
        super().__init__()

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/not_moving_mirs_cb_vs_n_log2.txt",
                     self.expression_map,
                     "not_moving")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/down_mirs_cb_vs_n_log2.txt",
                     self.expression_map,
                     "down")

        load_mir_exp("/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/references/agilent_mir/up_mirs_cb_vs_n_log2.txt",
                     self.expression_map,
                     "up")


class CustomExpression:
    """
    Add custom expression to a gene file.
    """

    def __init__(self, type, expression):
        self.type = type
        self.expression = expression

    def annotate(self, file):
        f = open(file, 'r')

        # skip header
        header = f.readline().strip().split("\t")
        header.append(self.type)

        print("\t".join(header))

        entrez_column = text.find_index(
            header, headings.ENTREZ_ID)

        for line in f:
            ls = line.strip()

            if len(ls) == 0:
                continue

            tokens = ls.split("\t")

            entrez = tokens[entrez_column]

            tokens.append(self.expression.get_expression(entrez))

            print("\t".join(tokens))

        f.close()