from typing import Any, Mapping

from . import genomic

class SearchMirs(genomic.BlockSearch):
    def __init__(self, file:str):
        super().__init__()

        f = open(file, 'r')

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')

            location = genomic.Feature(tokens[0], int(tokens[1]), int(tokens[2]))
            location.add_id('name', tokens[3])

            # we store locations for determining if they overlap or not
            self.add_feature(location, location)

        f.close()

    def get_overlapping_mirs(self, location: genomic.Location):
        features = self.get_features(location)

        ret = []

        for fs in features:
            for feature in fs:
                overlap = genomic.overlap_locations(location, feature)
                if overlap is not None:
                    ret.append(feature)

        return ret


class MirAnnotation(genomic.Annotation):
    def __init__(self, mirs:SearchMirs):
        super().__init__('mir')
        self._mirs = mirs

    def get_names(self):
        return ['miR Symbol']

    def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
        mirs = self._mirs.get_overlapping_mirs(location)

        ret = [';'.join(sorted([mir.get_id('name') for mir in mirs]))]

        return ret


# class MirAnnotation(genomic.Annotation):
#     """
#     Loads a mir database.
#     """

#     def __init__(self, refseq_genes:genes.RefSeqGenes, mir_max_gap:int, bin_size:int):
#         super().__init__('mir')

#         self._refseq_genes = refseq_genes
#         self._mir_max_gap = mir_max_gap
#         self._bin_size = bin_size
#         self._mir_types = collections.defaultdict(
#             lambda: collections.defaultdict(lambda: collections.defaultdict(str)))
#         self._mir_starts = collections.defaultdict(
#             lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
#         self._mir_ends = collections.defaultdict(
#             lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
#         self._mir_strands = collections.defaultdict(
#             lambda: collections.defaultdict(lambda: collections.defaultdict(str)))

#         self._precursors = collections.defaultdict(
#             lambda: collections.defaultdict(bool))
#         self._mature = collections.defaultdict(
#             lambda: collections.defaultdict(bool))

#         self._header = ['miR Symbol', 'Relative To miR',
#                         'miR Start Closest Distance', 'miR Start Distance']

#         print('Loading miRs...', file=sys.stderr)

#         f = open(
#             '/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirbase_19_hg19_genomic_locations.txt', 'r')

#         # skip header
#         f.readline()

#         for line in f:
#             line = line.strip()

#             if len(line) == 0:
#                 continue

#             tokens = line.split('\t')

#             name = tokens[0]
#             group = tokens[1]
#             type = tokens[2]
#             chr = tokens[3]
#             start = int(tokens[4])
#             end = int(tokens[5])
#             strand = tokens[6]

#             # if type != "miRNA_primary_transcript":
#             #  continue

#             # only annotate the mature
#             if type != 'miRNA':
#                 continue

#             start_bin = int(start / bin_size)
#             end_bin = int(end / bin_size)

#             for bin in range(start_bin, end_bin + 1):
#                 self._mir_types[chr][bin][name] = type
#                 self._mir_starts[chr][bin][name] = start
#                 self._mir_ends[chr][bin][name] = end
#                 self._mir_strands[chr][bin][name] = strand

#         f.close()

#         print('Loading Transcript miRs...', file=sys.stderr)

#         f = open('/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirbase_19_transcript_locations_hg19_entrez.txt', 'r')

#         # skip header
#         f.readline()

#         for line in f:
#             line = line.strip()

#             if len(line) == 0:
#                 continue

#             tokens = line.split('\t')

#             mir = tokens[0]
#             entrez = tokens[6]

#             if entrez == text.NA:
#                 continue

#             self._precursors[entrez][mir] = True

#         f.close()

#         #
#         # open a mapping between precursor and mature
#         #

#         f = open(HSA_FILE, 'r')

#         id_precursor_map = collections.defaultdict(str)

#         for line in f:
#             if re.match(r'^#.*', line):
#                 continue

#             line = line.strip()

#             if len(line) == 0:
#                 continue

#             tokens = line.split('\t')

#             type = tokens[2]

#             matcher = re.match(r'^ID=([^\;]+).*Name=([^\;]+).*', tokens[8])

#             id = matcher.group(1)
#             mir = matcher.group(2)

#             if type == 'miRNA_primary_transcript':
#                 id_precursor_map[id] = mir
#             else:
#                 matcher = re.match(r'.*derives_from=([^\;]+).*', tokens[8])

#                 precursor_id = matcher.group(1)

#                 precursor = id_precursor_map[precursor_id]

#                 self._mature[precursor][mir] = True

#         f.close()

#     def get_names(self):
#         return self._header

#     def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
#         genes = row_map['genes']

#         mir_annotations = collections.defaultdict(str)
#         mss = collections.defaultdict(str)
#         self._annotate_location(
#             location, self._mir_max_gap, mir_annotations, mss)

#         mirs = []
#         mir_types = []
#         mss_list = []

#         #
#         # If the peak is a promoter, check to see if there is a mir in the gene
#         #

#         for gene in genes:
#             gene_annotation = self._refseq_genes.get_gene(gene.id)

#             promoter = False

#             for type in gene.types:
#                 if re.match(r'.*promoter.*', type):
#                     promoter = True
#                     break

#             if promoter == False:
#                 continue

#             transcript_mirs = self._find_mirs(gene_annotation.entrez)

#             if len(transcript_mirs) == 0:
#                 continue

#             symbol = gene_annotation.symbol

#             if symbol == text.NA:
#                 symbol = 'gene'

#             for mir in sorted(transcript_mirs):
#                 mirs.append(mir)
#                 mir_types.append(symbol + '_promoter')
#                 mss_list.append(text.NA)

#         # Standard mir annotation

#         for mir in sorted(mir_annotations):
#             mirs.append(mir)
#             mir_types.append(mir_annotations[mir])
#             mss_list.append(mss[mir])

#         ret = []

#         if len(mirs) > 0:
#             closest_mss = genomic.get_closest_tss(mss_list)

#             ret.append(';'.join(mirs))
#             ret.append(';'.join(mir_types))
#             ret.append(closest_mss)
#             ret.append(';'.join(mss_list))
#         else:
#             ret.append(text.NA)
#             ret.append(text.NA)
#             ret.append(text.NA)
#             ret.append(text.NA)

#         return ret

#     def _annotate_location(self, location, max_gap, annotation_types, mss):
#         if not chr in self._mir_types:
#             return

#         mid_point = int((location.start + location.end) / 2)

#         start_bin = int((location.start - max_gap) / self._bin_size)
#         end_bin = int((location.end + max_gap) / self._bin_size)

#         #sys.stderr.write('mir ' + chr + ' ' + str(mid_point) + '\n')

#         for bin in range(start_bin, end_bin + 1):
#             if not bin in self._mir_starts[chr]:
#                 continue

#             for name in self._mir_starts[chr][bin]:
#                 mir_start = self._mir_starts[chr][bin][name]
#                 mir_end = self._mir_ends[chr][bin][name]
#                 mir_type = self._mir_types[chr][bin][name]
#                 mir_strand = self._mir_strands[chr][bin][name]

#                 found = True

#                 if mir_strand == '+':
#                     if mir_start - mid_point <= max_gap and mir_start - mid_point > 0:
#                         annotation_types[name] = ','.join([mir_type, 'before'])
#                     elif mid_point - mir_end <= max_gap and mid_point - mir_end > 0:
#                         annotation_types[name] = ','.join([mir_type, 'after'])
#                     elif mid_point >= mir_start and mid_point <= mir_end:
#                         annotation_types[name] = ','.join([mir_type, 'over'])
#                     else:
#                         found = False

#                     if found == True:
#                         mss[name] = str(mid_point - mir_start)
#                     else:
#                         mss[name] = text.NA
#                 else:
#                     # neg strand
#                     if mir_start - mid_point <= max_gap and mir_start - mid_point > 0:
#                         annotation_types[name] = ",".join([mir_type, "after"])
#                     elif mid_point - mir_end <= max_gap and mid_point - mir_end > 0:
#                         annotation_types[name] = ",".join([mir_type, "before"])
#                     elif mid_point >= mir_start and mid_point <= mir_end:
#                         annotation_types[name] = ",".join([mir_type, "over"])
#                     else:
#                         found = False

#                     if found == True:
#                         mss[name] = str(mir_end - mid_point)
#                     else:
#                         mss[name] = text.NA

#     def _find_mirs(self, entrez):
#         ret = []

#         if entrez not in self._precursors:
#             return ret

#         for precursor in sorted(self._precursors[entrez]):
#             for mature in self._mature[precursor]:
#                 ret.append(mature)

#         return ret


# class MirTranscriptAnnotation:
#     """
#     Can be used to check if a gene contains a mir or not
#     """

#     def __init__(self):
#         self.precursors = collections.defaultdict(
#             lambda: collections.defaultdict(bool))
#         self.mature = collections.defaultdict(
#             lambda: collections.defaultdict(bool))

#         print('Loading Transcript miRs...', file=sys.stderr)

#         f = open('/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirbase_19_transcript_locations_hg19_entrez.txt', 'r')

#         # skip header
#         f.readline()

#         for line in f:
#             line = line.strip()

#             if len(line) == 0:
#                 continue

#             tokens = line.split('\t')

#             mir = tokens[0]
#             entrez = tokens[6]

#             if entrez == text.NA:
#                 continue

#             self.precursors[entrez][mir] = True

#         f.close()

#         #
#         # open a mapping between precursor and mature
#         #

#         f = open(
#             '/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/hsa.gff3', 'r')

#         id_precursor_map = collections.defaultdict(str)

#         for line in f:
#             if re.match(r'^#.*', line):
#                 continue

#             line = line.strip()

#             if len(line) == 0:
#                 continue

#             tokens = line.split('\t')

#             type = tokens[2]

#             matcher = re.match(r'^ID=([^\;]+).*Name=([^\;]+).*', tokens[8])

#             id = matcher.group(1)
#             mir = matcher.group(2)

#             if type == 'miRNA_primary_transcript':
#                 id_precursor_map[id] = mir
#             else:
#                 matcher = re.match(r'.*derives_from=([^\;]+).*', tokens[8])

#                 precursor_id = matcher.group(1)

#                 precursor = id_precursor_map[precursor_id]

#                 self.mature[precursor][mir] = True

#         f.close()

#     def find_mirs(self, entrez):
#         ret = []

#         if entrez not in self.precursors:
#             return ret

#         for precursor in sorted(self.precursors[entrez]):
#             for mature in self.mature[precursor]:
#                 ret.append(mature)

#         return ret
