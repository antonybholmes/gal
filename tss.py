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
from typing import Any, Mapping, Optional, Union
import numpy as np

from . import genomic
from . import genes
from . import headings
from . import text

NO_GENE_ROW = np.array([text.NA] * 5)


class TssGene:
	def __init__(self, id, types, d):
		self.id = id
		self.types = types
		self.d = d


class TssGenes:
	def __init__(self, d, genes):
		self.d = d
		self.genes = genes


class RefSeqTssClassification:
	"""
	Classify genes by type
	"""

	def __init__(self, file, prom_ext_5p, prom_ext_3p):
		self._variants = collections.defaultdict(set)
		self._gene_strands = collections.defaultdict()
		self._gene_start_map = collections.defaultdict(int)
		self._gene_end_map = collections.defaultdict(int)
		self._promoter_starts = collections.defaultdict(int)
		self._promoter_ends = collections.defaultdict(int)
		self._exon_counts = collections.defaultdict(int)
		self._exon_starts = collections.defaultdict(list)
		self._exon_ends = collections.defaultdict(list)

		# how to label promoters
		# "promoter_-" + str(prom_ext_5p / 1000) + "/+" + str(prom_ext_3p / 1000)
		self._promoter_label = 'promoter'

		print(
			f'RefSeqTssClassification: Loadings from {file}...', file=sys.stderr)

		f = open(file, 'r')

		# skip header
		f.readline()

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
			exon_count = int(tokens[7])

			if re.match(r'.*n/a.*', refseq):
				continue

			if re.match(r'.*n/a.*', entrez):
				continue

			if re.match(r'.*n/a.*', symbol):
				continue

			if re.match(r'.*MIR.*', symbol):
				continue

			# UCSC convention
			ex_starts = [int(p) + 1 for p in tokens[8].split(',')]
			ex_ends = [int(p) for p in tokens[9].split(',')]

			if strand == '+':
				promoter_start = start - prom_ext_5p
				promoter_end = start + prom_ext_3p
			else:
				promoter_start = end - prom_ext_3p
				promoter_end = end + prom_ext_5p

			variant_id = genes.create_variant_id(
				refseq, chr, start, end)

			self._variants[refseq].add(variant_id)
			self._gene_strands[variant_id] = strand
			self._gene_start_map[variant_id] = start
			self._gene_end_map[variant_id] = end
			self._promoter_starts[variant_id] = promoter_start
			self._promoter_ends[variant_id] = promoter_end
			self._exon_counts[variant_id] = exon_count

			for p in ex_starts:
				self._exon_starts[variant_id].append(p)

			for p in ex_ends:
				self._exon_ends[variant_id].append(p)

		f.close()

		print('Finished.', file=sys.stderr)

	def classify(self, variant_id, mid_point):
		refseq = genes.parse_id_from_variant(variant_id)

		types = set()

		# Test all variants with the same refseq and start as our variant
		# Although we picked a closest variant, we should test all variants
		# with the same start/end to determine the type of the gene

		min_d = -1
		min_abs_d = -1

		for var_id in self._variants[refseq]:
			strand = self._gene_strands[var_id]

			# skip variants that do not begin or end where the reference
			# variant does.
			if strand == '+':
				if self._gene_start_map[var_id] != self._gene_start_map[variant_id]:
					continue
			else:
				if self._gene_end_map[var_id] != self._gene_end_map[variant_id]:
					continue

			gene_start = self._gene_start_map[var_id]
			gene_end = self._gene_end_map[var_id]

			promoter_start = self._promoter_starts[var_id]
			promoter_end = self._promoter_ends[var_id]

			if strand == '+':
				d = mid_point - gene_start
			else:
				d = gene_end - mid_point

			abs_d = abs(d)

			if abs_d < min_abs_d or min_abs_d == -1:
				min_abs_d = abs_d
				min_d = d

			# Peak must fall within gene boundaries
			within_bounds = False

			if strand == '+':
				if mid_point >= promoter_start and mid_point <= gene_end:
					within_bounds = True
			else:
				# negative strand
				if mid_point >= gene_start and mid_point <= promoter_end:
					within_bounds = True

			if within_bounds:
				# Peak must be within gene boundary to be called promoter, exon or intron
				if mid_point >= promoter_start and mid_point <= promoter_end:
					types.add(self._promoter_label)

				# exonic
				in_exon = False

				for i in range(0, len(self._exon_starts[var_id])):
					if mid_point >= self._exon_starts[var_id][i] and mid_point <= self._exon_ends[var_id][i]:
						types.add('exonic')
						in_exon = True
						break

				# if you're not exonic and you are within the gene, you must
				# be intronic (in gene)
				if not in_exon and mid_point >= gene_start and mid_point <= gene_end:
					types.add('intronic')

		if 'exonic' in types and 'intronic' in types:
			# We favor classifying as intronic where possible
			types.remove('exonic')

		# classificatio of promoter really means promoter intergenic (outside gene tx)

		if len(types) == 0:
			# No types implies integenic
			types.add('intergenic')

		gene = TssGene(variant_id, types, min_d)

		# Return in sorted so that promoter comes first
		# return reversed(sorted(types))

		return gene


class RefSeqTssClassificationFactory:
	_classifiers = collections.defaultdict(
		lambda: collections.defaultdict(lambda: {}))

	@staticmethod
	def getInstance(file: str, prom_ext_5p: int, prom_ext_3p: int):
		if prom_ext_3p not in RefSeqTssClassificationFactory._classifiers[file][prom_ext_5p]:
			RefSeqTssClassificationFactory._classifiers[file][prom_ext_5p][prom_ext_3p] = RefSeqTssClassification(
				file, prom_ext_5p, prom_ext_3p)

		return RefSeqTssClassificationFactory._classifiers[file][prom_ext_5p][prom_ext_3p]

	def __init__(self):  # , prom_ext_5p, prom_ext_3p):
		# super().__init__(annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)
		raise Exception("Call RefSeqTssClassificationFactory.getInstance()")


class RefSeqAnnotation:
	"""
	Annotate peaks for all possible refseq genes they might overlap.
	"""

	def __init__(self, file: str, prom_ext_5p: int = 2000, prom_ext_3p: int = 1000, bin_size: int = 1000):
		self._metadata = {}

		"""
		Create a new annotation object.

		Args:
			file (str): _description_
			prom_ext_5p (int, optional): How far upstream should be considered a promoter. Defaults to 2000.
			prom_ext_3p (int, optional): How far downstream should be considered a promoter. Defaults to 1000.
			bin_size (int, optional): The size of the bins to group genes into. Defaults to 1000.
		"""

		self._bin_size = bin_size

		self._transcript_to_gene_map = collections.defaultdict(str)
		self._gene_strands = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(str)))
		self._gene_start_map = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
		self._gene_end_map = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
		self._promoter_starts = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
		self._promoter_ends = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
		self._exon_counts = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
		self.exon_starts = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(list[int])))
		self._exon_ends = collections.defaultdict(
			lambda: collections.defaultdict(lambda: collections.defaultdict(list[int])))
		self._refseq_class = RefSeqTssClassificationFactory.getInstance(
			file, prom_ext_5p, prom_ext_3p)

		print(f'RefSeqAnnotation: Loading from {file}...', file=sys.stderr)

		self._metadata['db'] = file
		self._metadata['prom_ext_5p'] = prom_ext_5p
		self._metadata['prom_ext_3p'] = prom_ext_3p
		self._metadata['bin_size'] = bin_size
		self._metadata['version'] = re.search(
			r'v([A-Za-z0-9]+)', file).group(1)

		f = open(file, 'r')

		# skip header
		f.readline()

		for line in f:
			line = line.strip()

			if len(line) == 0:
				continue

			tokens = line.split("\t")
			refseq = tokens[0]
			entrez = tokens[1]
			symbol = tokens[2]
			chr = tokens[3]
			strand = tokens[4]
			# ucsc convention of having first coordinate be 0 based
			start = int(tokens[5]) + 1
			end = int(tokens[6])
			exon_count = int(tokens[7])

			if re.match(r'.*n/a.*', entrez):
				continue

			if re.match(r'.*n/a.*', symbol):
				continue

			if re.match(r'.*MIR.*', symbol):
				continue

			self._transcript_to_gene_map[refseq] = entrez

			ex_starts = [int(p) + 1 for p in tokens[8].split(",")]
			ex_ends = [int(p) for p in tokens[9].split(",")]

			max_offset = max(prom_ext_5p, prom_ext_3p)

			if strand == "+":
				promoter_start = start - prom_ext_5p
				promoter_end = start + prom_ext_3p
			else:
				promoter_start = end - prom_ext_3p
				promoter_end = end + prom_ext_5p

			start_bin = int((start - max_offset) / bin_size)
			end_bin = int((end + max_offset) / bin_size)

			# essentially transcript id
			transcript_id = genes.create_variant_id(
				refseq, chr, start, end)

			for bin in range(start_bin, end_bin + 1):
				# apparently refseq genes from the ucsc always report
				# coordinates on the forward strand regardless of orientation
				self._gene_strands[chr][bin][transcript_id] = strand
				self._gene_start_map[chr][bin][transcript_id] = start
				self._gene_end_map[chr][bin][transcript_id] = end
				self._promoter_starts[chr][bin][transcript_id] = promoter_start
				self._promoter_ends[chr][bin][transcript_id] = promoter_end
				self._exon_counts[chr][bin][transcript_id] = exon_count

				for p in ex_starts:
					self.exon_starts[chr][bin][transcript_id].append(p)

				for p in ex_ends:
					self._exon_ends[chr][bin][transcript_id].append(p)

		f.close()

		print('Finished.', file=sys.stderr)

	@property
	def metadata(self):
		return self._metadata

	def annotate_location(self, location):
		if not location.chr in self._gene_start_map:
			return []

		mid_point = int((location.start + location.end) / 2)

		start_bin = int(location.start / self._bin_size)
		end_bin = int(location.end / self._bin_size)

		# first find all the variants we might belong to
		closest_transcripts = collections.defaultdict()
		closest_d = collections.defaultdict(int)
		closest_abs_d = collections.defaultdict(int)

		for bin in range(start_bin, end_bin + 1):
			if not bin in self._gene_start_map[location.chr]:
				continue

			for transcript_id in self._gene_start_map[location.chr][bin]:
				refseq = genes.parse_id_from_variant(transcript_id)

				gene_id = self._transcript_to_gene_map[refseq]

				strand = self._gene_strands[location.chr][bin][transcript_id]

				gene_start = self._gene_start_map[location.chr][bin][transcript_id]
				gene_end = self._gene_end_map[location.chr][bin][transcript_id]

				#
				# Deal with a peak being in a promoter
				#

				promoter_start = self._promoter_starts[location.chr][bin][transcript_id]
				promoter_end = self._promoter_ends[location.chr][bin][transcript_id]

				absd = -1

				# we are considered as being within a gene if we fall inside
				# the gene coordinates with the promoter extension, i.e. even
				# if outside the gene a bit, we are still considered overlapping
				# the gene if in the promoter region
				if strand == "+":
					# for short genes we need to check the max distance of
					# gene end or promoter end

					if mid_point >= promoter_start and mid_point <= gene_end:
						d = mid_point - gene_start
						absd = abs(d)
				else:
					# negative strand
					if mid_point >= gene_start and mid_point <= promoter_end:
						d = gene_end - mid_point
						absd = abs(d)

				if absd != -1:
					# we are in a variant so record the closest in the entrez group
					if gene_id in closest_abs_d:
						if absd < closest_abs_d[gene_id]:  # closest_d[entrez]
							closest_transcripts[gene_id] = transcript_id
							closest_d[gene_id] = d
							closest_abs_d[gene_id] = absd
					else:
						closest_transcripts[gene_id] = transcript_id
						closest_d[gene_id] = d
						closest_abs_d[gene_id] = absd

		# Now classify each separate closest variant/transcript (which must all have a different entrez id)

		ret = []

		for gene_id in sorted(closest_transcripts):
			transcript_id = closest_transcripts[gene_id]
			gene = self._refseq_class.classify(transcript_id, mid_point)
			ret.append(gene)

		return ret


class RefSeqAnnotationFactory:
	_classifiers = collections.defaultdict(
		lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: {})))

	@staticmethod
	def getInstance(file: str, prom_ext_5p: int = 2000, prom_ext_3p: int = 1000, bin_size: int = 1000):
		if bin_size not in RefSeqAnnotationFactory._classifiers[file][prom_ext_5p][prom_ext_3p]:
			RefSeqAnnotationFactory._classifiers[file][prom_ext_5p][prom_ext_3p][bin_size] = RefSeqAnnotation(
				file, prom_ext_5p, prom_ext_3p, bin_size)

		return RefSeqAnnotationFactory._classifiers[file][prom_ext_5p][prom_ext_3p][bin_size]

	def __init__(self):
		raise Exception("Call RefSeqAnnotationFactory.getInstance()")


class RefSeqTss(genomic.Annotation):
	"""
	Determines the closest gene(s) to a given position.
	"""

	def __init__(self, name: str, file: str, refseq_genes: genes.RefSeqGenes, prom_ext_5p: int = 2000, prom_ext_3p: int = 1000, bin_size: Optional[int] = 1000):
		super().__init__(name)
		self._refseq_genes = refseq_genes
		self._gene_start_map = collections.defaultdict(
			lambda: collections.defaultdict(set))
		self._gene_index_map = collections.defaultdict(
			lambda: collections.defaultdict(str))
		self._starts = collections.defaultdict(list)
		self._bin_size = bin_size
		self._bins = collections.defaultdict(
			lambda: collections.defaultdict(list))
		self._refseq = RefSeqTssClassificationFactory.getInstance(
			file, prom_ext_5p, prom_ext_3p)
		# f'(prom=-{prom_ext_5p / 1000}/+{prom_ext_3p / 1000}kb)'
		self._promoter_type = genomic.promoter_type(prom_ext_5p, prom_ext_3p)
		self._header = []

	def get_names(self):
		return self._header

	def update_row(self, location: genomic.Location, row_map: Mapping[str, Any]):
		#print(location, file=sys.stderr)
		
		tss_gene = self.get_closest_gene(location)

		if tss_gene is None:
			return NO_GENE_ROW

		gene_annotation = self._refseq_genes.get_gene(tss_gene.id)

		ret = [gene_annotation.refseq, gene_annotation.entrez, gene_annotation.symbol, ','.join(
			reversed(sorted(tss_gene.types))), str(tss_gene.d)]

		return ret

	def get_closest_gene(self, location: genomic.Location) -> TssGene:
		"""
		Returns the genes that are closest to a given location.

		Args:
				chr:     The chromosome.
				start:   The start position.
				end:     The end position
		Returns:
				A list of genes with the distance of their TSS from the peak.
		"""
		mid_point = genomic.mid_point(location)

		# use a binary style search to find closest peak,
		# we keep track of closest end and start we find
		# hence we do not update start and end by adding
		# or subtracting 1, but instead use the test point
		# as the new range end or start.
		# Once we have the closest start before the test
		# point and the closest end, we pick the closest
		# of the two.

		n = len(self._starts[location.chr])

		if n == 0:
			return None

		ps = 0
		pe = n - 1

		#print(location, file=sys.stderr)

		if mid_point < self._starts[location.chr][ps]:
			return self.get_annotation(location.chr, mid_point, self._starts[location.chr][ps])

		if mid_point > self._starts[location.chr][pe]:
			return self.get_annotation(location.chr, mid_point, self._starts[location.chr][pe])

		while pe - ps > 1:  # while pe >= ps:
			pm = int((ps + pe) / 2)

			test_start = self._starts[location.chr][pm]

			#print('test', ps, pe, pm, test_start, self._gene_start_map[location.chr][test_start], file=sys.stderr)

			if test_start > mid_point:
				pe = pm  # - 1
			elif test_start < mid_point:
				ps = pm  # + 1
			else:
				# perfect match
				return self.get_annotation(location.chr, mid_point, test_start)

		# The peak lies between two starts, pick the closest
		abss = abs(mid_point - self._starts[location.chr][ps])
		abse = abs(mid_point - self._starts[location.chr][pe])

		#print('hmm', ps, pe, abss, abse, file=sys.stderr)

		if abss <= abse:
			start = self._starts[location.chr][ps]
		else:
			start = self._starts[location.chr][pe]

		return self.get_annotation(location.chr, mid_point, start)

	def get_annotation(self, chr, mid_point, start) -> TssGene:
		# Pick the first one we find
		variant_id = sorted(self._gene_start_map[chr][start])[0]
		gene = self._refseq.classify(variant_id, mid_point)
		return gene

	def get_n_closest_genes(self, location: genomic.Location, n: int = 5) -> list[TssGene]:
		"""
		Find the n closest genes to a location.

		Args:
				location:   a location object
				n:          number of closest genes to return
		Returns:
				A list of closest gene objects.
		"""

		closest_gene = self.get_closest_gene(location)

		if closest_gene is None:
			return np.array([])

		# find index of gene in list of starts
		idx = self._gene_index_map[location.chr][closest_gene.id]

		#print(closest_gene.id, idx, file=sys.stderr)

		# pick 25 start positions before and 25 after
		idx1 = max(0, idx - 25)
		idx2 = min(self._starts[location.chr].size, idx + 25)

		before_variants = []

		for start in self._starts[location.chr][idx1:idx]:
			before_variants.extend(self._gene_start_map[location.chr][start])

		before_variants = np.array(before_variants)

		after_variants = []

		for start in self._starts[location.chr][idx:idx2]:
			after_variants.extend(self._gene_start_map[location.chr][start])

		after_variants = np.array(after_variants)

		variants = np.concatenate([before_variants, after_variants])

		mid_point = genomic.mid_point(location)

		genes = [self._refseq.classify(v, mid_point) for v in variants]

		gene_dist_map = collections.defaultdict(list)

		for gene in genes:
			gene_dist_map[abs(gene.d)].append(gene)

		genes_by_dist = []

		for d in sorted(gene_dist_map):
			genes_by_dist.extend(gene_dist_map[d])

		ret = [closest_gene]
		ret.extend(genes_by_dist[0:min(len(genes_by_dist), n-1)])

		return np.array(ret)


class RefSeqStart(RefSeqTss):
	"""
	Determines the closest gene(s) to a given position by TSS.
	"""

	def __init__(self, file, refseq_genes: genes.RefSeqGenes, prom_ext_5p: int = 2000, prom_ext_3p: int = 1000):
		super().__init__('refseq-tss-start', file, refseq_genes, prom_ext_5p, prom_ext_3p)

		print(f'RefSeqStart: Loading from {file}...', file=sys.stderr)

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

			if re.match(r'.*n/a.*', refseq):
				continue

			if re.match(r'.*n/a.*', entrez):
				continue

			if re.match(r'.*n/a.*', symbol):
				continue

			if re.match(r'.*MIR.*', symbol):
				continue

			# If the gene is on the negative strand, the promoter is around the
			# end coordinate on the forward strand

			variant_id = genes.create_variant_id(
				refseq, chr, start, end)

			if strand == "+":
				self._gene_start_map[chr][start].add(variant_id)
			else:
				self._gene_start_map[chr][end].add(variant_id)

		f.close()

		# Create sorted lists of start locations for each gene
		for chr in self._gene_start_map:
			self._starts[chr] = np.array(sorted(self._gene_start_map[chr]))

		# create index so we can determine where in the list each variant is

		for chr in self._starts:
			for i in range(0, len(self._starts[chr])):
				start = self._starts[chr][i]
				for variant_id in self._gene_start_map[chr][start]:
					self._gene_index_map[chr][variant_id] = i

				#bin = int(start / self._bin_size)
				# self._bins[chr][bin].append(start)

		print('Finished.', file=sys.stderr)

		self._header.append(f'Closest {headings.REFSEQ_ID}')
		self._header.append(f'Closest {headings.ENTREZ_ID}')
		self._header.append(headings.CLOSEST_GENE_SYMBOL)
		self._header.append(
			f'Relative To Closest Gene ({self._promoter_type})')
		self._header.append(f'TSS Closest Distance')


class RefSeqEnd(RefSeqTss):
	"""
	Find the closest gene to a given location using the genes end. Secondary to
	RefSeqStart.
	"""

	def __init__(self, file, refseq_genes: genes.RefSeqGenes, prom_ext_5p: int = 2000, prom_ext_3p: int = 1000):
		super().__init__('refseq-tss-end', file, refseq_genes, prom_ext_5p, prom_ext_3p)

		print(
			f'RefSeqEnd: Loading from {file} {prom_ext_5p} {prom_ext_3p}...', file=sys.stderr)

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

			# If the gene is on the negative strand, the promoter is around the
			# end coordinate on the forward strand

			variant_id = genes.create_variant_id(
				refseq, chr, start, end)

			#sys.stderr.write("var id " + variant_id + "\n")

			if strand == "+":
				self._gene_start_map[chr][end].add(variant_id)
			else:
				self._gene_start_map[chr][start].add(variant_id)

		f.close()

		# Create sorted lists of start locations for each gene
		for chr in self._gene_start_map:
			self._starts[chr] = sorted(self._gene_start_map[chr])

		# create index so we can determine where in the list each variant is

		for chr in self._starts:
			for i in range(0, len(self._starts[chr])):
				start = self._starts[chr][i]
				for variant_id in self._gene_start_map[chr][start]:
					self._gene_index_map[chr][variant_id] = i

				#bin = int(start / self._bin_size)
				# self._bins[chr][bin].append(start)

		print('Finished.', file=sys.stderr)

		self._header.append(f'Closest End {headings.REFSEQ_ID}')
		self._header.append(f'Closest End {headings.ENTREZ_ID}')
		self._header.append(f'Closest End {headings.GENE_SYMBOL}')
		self._header.append(
			f'Relative To Closest Gene End ({self._promoter_type})')
		self._header.append(f'End Closest Distance')


class OverlapTss(genomic.Annotation):
	"""
	Determines if a region overlaps a TSS
	"""

	def __init__(self, file: str, block_size: int = 100):
		super().__init__('overlap-tss')

		self._search = genomic.BlockSearch(block_size)

		print(f'OverlapTss: Loading from {file}...', file=sys.stderr)

		f = open(file, 'r')

		# skip header
		f.readline()

		# To account for multiple versions of a gene, allocate each entrez
		# id a unique index

		for line in f:
			line = line.strip()

			if len(line) == 0:
				continue

			tokens = line.split("\t")

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

			# If the gene is on the negative strand, the promoter is around the
			# end coordinate on the forward strand

			if strand == '+':
				s = start
			else:
				s = end

			# only interested in tss
			gene = genes.RefSeqGene(
				refseq, entrez, symbol, strand, chr, s, s)

			self._search.add_feature(gene, gene)

		f.close()

		print('Finished.', file=sys.stderr)

	def get_names(self):
		return ["Overlaps TSS"]

	def update_row(self, location: genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
		featureSets = self._search.get_features(location)

		ret = set()

		for featureSet in featureSets:
			for gene in featureSet:
				if genomic.is_overlapping(location, gene):
					ret.add(gene.symbol)

		if len(ret) == 0:
			ret.add(text.NA)

		return [';'.join(sorted(ret))]
