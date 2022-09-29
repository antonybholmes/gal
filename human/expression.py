import sys

from .. import text
from .. import headings
from .. import expression


class PeakExpression(expression.PeakExpression):
    """
    Add expression to peaks with gene symbols
    """

    def __init__(self):
        self.expression_list = []
        self.expression_list_headers = []

        #self.affy_gene_cb_vs_m_expression = expression.AffyGeneCBvsMExpression()
        #self.affy_gene_cb_vs_n_expression = expression.AffyGeneCBvsNExpression()
        #self.rna_gene_cb_vs_m_expression = expression.RnaSeqGeneCBvsMExpression()
        #self.rna_gene_cb_vs_n_expression = expression.RnaSeqGeneCBvsNExpression()

        self.expression_list_headers.append("GEP Affy GCvsN")
        self.expression_list_headers.append("GEP Affy GCvsM")
        self.expression_list.append(expression.AffyGeneCBvsNExpression())
        self.expression_list.append(expression.AffyGeneCBvsMExpression())

        self.expression_list_headers.append("RNA-seq GCvsN")
        self.expression_list_headers.append("RNA-seq GCvsM")
        self.expression_list.append(expression.RnaSeqGeneCBvsNExpression())
        self.expression_list.append(expression.RnaSeqGeneCBvsMExpression())

        # RNA-seq from deseq2
        self.expression_list_headers.append("RNA-seq GCvsN DeSeq2")
        self.expression_list_headers.append("RNA-seq GCvsM DeSeq2")
        self.expression_list.append(
            expression.RnaSeqGeneCBvsNDeSeq2Expression())
        self.expression_list.append(
            expression.RnaSeqGeneCBvsMDeSeq2Expression())

        #
        # SUD10 stuff
        #
        self.expression_list_headers.append("RNA-seq SUD10 WTvsD83V")
        self.expression_list_headers.append("RNA-seq SUD10 WTvsSTOP")
        self.expression_list_headers.append("RNA-seq SUD10 D83VvsSTOP")
        self.expression_list.append(
            expression.RnaSeqGeneSUD10WTvsD83vExpression())
        self.expression_list.append(
            expression.RnaSeqGeneSUD10WTvsStopExpression())
        self.expression_list.append(
            expression.RnaSeqGeneSUD10D83vsStopExpression())

        self.expression_list_headers.append("RNA-seq SUD10 WTvsD83V DeSeq2")
        self.expression_list_headers.append("RNA-seq SUD10 WTvsSTOP DeSeq2")
        self.expression_list_headers.append("RNA-seq SUD10 D83VvsSTOP DeSeq2")
        self.expression_list.append(
            expression.RnaSeqGeneSUD10WTvsD83vDeSeq2Expression())
        self.expression_list.append(
            expression.RnaSeqGeneSUD10WTvsStopDeSeq2Expression())
        self.expression_list.append(
            expression.RnaSeqGeneSUD10D83vsStopDeSeq2Expression())

        #
        # Mouse MEF2B vs WT
        #

        self.expression_list_headers.append("GEP Affy MEF2B Mouse WTvsKO")
        self.expression_list.append(expression.GEPMEF2BMouseWTvsKOExpression())


    def annotate(self, file):
        f = open(file, 'r')

        # skip header
        header = f.readline().strip().split("\t")

        location_column = text.find_index(header, headings.LOCATION)
        entrez_column = text.find_index(header, headings.ENTREZ_ID)
        symbol_column = text.find_index(header, headings.GENE_SYMBOL)

        # Print header
        sys.stdout.write("\t".join(header))

        for header in self.expression_list_headers:
            sys.stdout.write("\t" + header)

        sys.stdout.write("\n")

        for line in f:
            line = line.strip()

            tokens = line.split("\t")

            if entrez_column != -1:
                ids = sorted(tokens[entrez_column].split(";"))
            else:
                ids = sorted(tokens[symbol_column].split(";"))

            sys.stdout.write(line)

            for expression in self.expression_list:
                exp = []

                for id in ids:
                    exp.append(expression.get_expression(id))

                sys.stdout.write("\t" + ";".join(exp))

            sys.stdout.write("\n")
