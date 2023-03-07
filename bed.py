# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 16:39:28 2015

@author: antony
"""

import sys
import collections
import re

from . import sample
from .genomic import Feature, Location, GappedSearch, BlockSearch, GenomicFeaturesOverlap, overlap_locations

class SearchBed(BlockSearch):  # (lib.search.GappedSearch):
	"""
	Instance of gapped search for bed files allowing for named regions
	"""

	def __init__(self, file: str):
		super().__init__()  # GappedSearch

		f = open(file, 'r')

		for line in f:
			line = line.strip()

			if len(line) == 0:
					continue

			if 'track' in line:
					continue

			tokens = line.split('\t')

			location = Feature(tokens[0], int(tokens[1]), int(tokens[2]))

			if len(tokens) > 3:
				name = tokens[3]
			else:
				name = ''

			location.add_id('name', name)

			# we store locations for determining if they overlap or not
			self.add_feature(location, location)

		f.close()


# (lib.search.GappedSearch):
class SearchGenomicBedFeatures(BlockSearch):
	"""
	Instance of gapped search for bed files
	"""

	def __init__(self, bed_file: str):
		super().__init__()  # GappedSearch

		f = open(bed_file, 'r')

		for line in f:
			line = line.strip()

			if len(line) == 0:
				continue

			if 'track' in line:
				continue

			tokens = line.split('\t')

			location = Location(tokens[0], int(tokens[1]), int(tokens[2]))

			# we store locations for determining if they overlap or not
			self.add_feature(location, location)

		f.close()


class GenomicBedOverlap(GenomicFeaturesOverlap):
	"""
	Uses a gapped search to determine by how much a location overlaps a feature
	"""

	def __init__(self, gapped_search: GappedSearch):
		super().__init__(gapped_search)

	def get_regions(self, location: Location):
		features = self.get_features(location)

		ret = []

		for feature in features:
			for l in feature.values:
				ret.append(l)

		return ret

	def get_overlapping_regions(self, location: Location):
		regions = self.get_regions(location)

		ret = []

		for r in regions:
			overlap = overlap_locations(location, r[0])
			if overlap is not None:
				ret.append(r)

		return ret
	

def create_bed_header(name, color="255,0,0"):
	# Need to truncate the name to 16 chars
	
	if len(name) == 16:
		short_name = name
	else:
		short_name = name[0:13] + "..." #pychipseq.sample.get_sample_id(name)
	
	return "track type=bed name=\"" + short_name  + "\" description=\"" + name + " peaks\" color=\"" + color + "\""


def create_bedgraph_header(name, color="255,0,0"):
	# Need to truncate the name to 16 chars
	short_name = sample.get_sample_id(name)
	
	return "track type=bedGraph name=\"" + short_name  + "\" description=\"" + name + " peaks\" visibility=full autoScale=on alwaysZero=on color=\"" + color + "\""


def write_bed_line(chr, start, end, value):
	# UCSC convention
	sys.stdout.write("\t".join([chr, str(start - 1), str(end), value]) + "\n")

def write_bed_location(location):
	write_bed_line(location.chr, location.start, location.end, location.to_string())


def write_bedgraph_line(chr, start, end, value):
	# UCSC convention
	sys.stdout.write("\t".join([chr, str(start - 1), str(end), "{:.2f}".format(value)]) + "\n")


def load_bed_peaks(file):
	"""
	Loads just the peaks from a bed file
	"""
	
	sys.stderr.write("Reading from bed " + file + "\n")
	
	peaks = collections.defaultdict(float)
	
	f = open(file, "r")

	for line in f:
		line = line.strip()
		
		if len(line) == 0:
			continue

		if re.match(r'^track.*', line):
			continue
		
		tokens = line.split("\t")
		
		chr = tokens[0]
		start = int(tokens[1]) + 1
		end = int(tokens[2])
		
		location = chr + ":" + str(start) + "-" + str(end)
		
		rpm = float(tokens[3])
		
		peaks[location] = rpm
	
	f.close()
	
	return peaks


def load_bed(file:str):
	"""
	Loads BED file as locations
	"""
	
	print(f'Loading peaks from {file}...', file=sys.stdout)

	peaks = []
	
	f = open(file, "r")

	for line in f:
		line = line.strip()
		
		if len(line) == 0:
			continue

		if re.match(r'^track.*', line):
			continue
		
		tokens = line.split("\t")
		
		chr = tokens[0]
		start = int(tokens[1]) + 1
		end = int(tokens[2])
		
		location = Location(chr, start, end)
		
		peaks.append(location)
	
	f.close()
	
	return peaks
