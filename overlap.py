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
from operator import itemgetter
import sys

from . import genomic
from . import peaks

BIN_SIZE = 1000


class Overlap:
		"""
		Uses a gapped search to determine by how much a location overlaps a feature
		"""

		def __init__(self, gapped_search):
				self.gapped_search = gapped_search

		def get_max_overlap(self, location):
				features = self.gapped_search.get_features(location)

				max_overlap_width = -1
				max_overlap = None

				for feature in features:
						for l in feature.values:
								overlap = genomic.overlap_locations(location, l)

								if overlap is not None:
										if overlap.width > max_overlap_width:
												max_overlap_width = overlap.width
												max_overlap = overlap

				return max_overlap


class Overlapping:
		"""
		Finds the closest peak to another set of peaks.
		"""

		def __init__(self, file):
				self.bins = collections.defaultdict(
						lambda: collections.defaultdict(set))

				self.load(file)

		def load(self, file):
				sys.stderr.write("Loading peaks from " + file + "...\n")

				f = open(file, 'r')

				# skip header
				f.readline()

				for line in f:
						line = line.strip()

						if len(line) == 0:
								continue

						tokens = line.split("\t")

						location = genomic.parse_location(tokens[0])

						sbin = int(location.start / BIN_SIZE)
						ebin = int(location.end / BIN_SIZE)

						for i in range(sbin, ebin + 1):
								self.bins[location.chr][i].add(location)

				f.close()

		def overlaps(self, location):
				sbin = int(location.start / BIN_SIZE)
				ebin = int(location.end / BIN_SIZE)

				overlaps = []

				for i in range(sbin, ebin + 1):
						for l in self.bins[location.chr][i]:
								if genomic.is_overlapping(location, l):
										overlaps.append(l)

				return overlaps


def overlap(ref_file, query_file):
		peaks = peaks.parse_peaks(ref_file)

		f = open(query_file, 'r')

		# skip header
		# f.readline()

		lines = []

		for line in f:
				tokens = line.split('\t')

				chr = tokens[0]
				start = int(tokens[1])
				end = int(tokens[2])

				features = peak_overlap(peaks, chr, start, end)

				lines.append('\t'.join(
						[f'{chr}:{start}-{end}', str(len(features)), ';'.join(features)]))

		f.close()

		return lines


def peak_overlap(peaks, chr, start, end):
		features = []

		for id in sorted(peaks[chr]):
				p_start = int(peaks[chr][id][0])
				p_end = int(peaks[chr][id][1])

				if start >= p_start and end <= p_end:
						features.append(';'.join([id, 'within']))
				elif start < p_start and end > p_end:
						features.append(';'.join([id, 'over']))
				elif start < p_start and end > p_start:
						features.append(';'.join([id, 'upstream']))
				elif start < p_end and end > p_end:
						features.append(';'.join([id, 'downstream']))

		return features


def _overlapping(location_id_map: dict[str, str], location_map: dict[str, genomic.Location], location_bins: dict[int, list[str]], locations: list[str]) -> dict[str, dict[str, str]]:
	"""
	Calculates the maximum overlaps between a set of sample locations

	Args:
			location_id_map (dict[str, str]): _description_
			location_map (dict[str, genomic.Location]): _description_
			location_bins (dict[int, list[str]]): _description_
			locations (list[str]): _description_

	Returns:
			dict[str, dict[str, str]]: _description_
	"""
	
	# lets see what overlaps 
	
	location_core_map = collections.defaultdict(lambda: collections.defaultdict(str))

	# debug for testing to end remove as it truncates list
	# locations = locations[1:(len(locations) / 4)]

	total = len(locations)

	total *= total

	print(f'Processing {len(locations)}...', file=sys.stderr)

	p = 1

	# keep track of all locations that have been allocated at least once
	allocated = set()

	for location1 in locations:
			# of the form id=chrN:start-end, an lid

			# get its location
			loc1 = location_map[location1]
			# group1 = location_group_map[location1]

			# if group1 != 'none':
			#	continue

			# find possible overlapping locations

			test_locations = set()

			bin_start = int(loc1.start / BIN_SIZE)
			bin_end = int(loc1.end / BIN_SIZE)

			for bin in range(bin_start, bin_end + 1):
					for location in location_bins[bin]:
							test_locations.add(location)

			used = set()

			# Form the largest group of overlapping peaks
			exhausted = False

			while not exhausted:
					grouped_locations = [location1]

					# reset for each group search
					loc1 = location_map[location1]

					for location2 in test_locations:
							if location1 == location2:
									continue

							if location2 in used:
									continue

							# sys.stderr.write(location1 + ' ' + location2 + '\n')

							loc2 = location_map[location2]
							# group2 = location_group_map[location2]

							if p % 10000000 == 0:
									print(f'p ({p})', file=sys.stderr)

							p += 1

							# if group2 != 'none':
							#	continue

							if loc1.chr != loc2.chr:
									continue

							# Given the two locations we are testing, sort them so
							# the starts and ends are in order.
							if (loc1.start <= loc2.start):
									min_loc = loc1
									max_loc = loc2
							else:
									min_loc = loc2
									max_loc = loc1

							overlap = -1
							overlap_start = -1

							if max_loc.start == min_loc.start and max_loc.end > min_loc.end:
									# Peak 1 and 2 start at the same point but peak 2 is wider
									# so the overlap region is peak 1
									overlap = min_loc.end - min_loc.start + 1
									overlap_start = min_loc.start
							elif max_loc.start >= min_loc.start and max_loc.end <= min_loc.end:
									# Peak 1 is wider than peak 2 and contains it
									# so the overlap region is peak 2
									overlap = max_loc.end - max_loc.start + 1
									overlap_start = max_loc.start
							elif min_loc.start < max_loc.start and min_loc.end > max_loc.start:
									# Peak 1 starts before peak 2 but ends within peak 2
									# so the overlap is the start of peak 2 to the end of
									# peak 1
									overlap = min_loc.end - max_loc.start + 1
									overlap_start = max_loc.start
							else:
									pass

							# We have not found an overlap yet so continue
							if overlap == -1:
									continue

							# change the start1 and end1 coordinates to reflect the overlap
							# region so that each subsequent match must be within this region
							# this prevents long peaks that overlap two smaller peaks who
							# themselves do not overlap each other
							loc1 = genomic.Location(
									loc1.chr, overlap_start, overlap_start + overlap - 1)

							grouped_locations.append(location2)

					# now we have a list of all locations that overlap each other

					# if we have a group of entries, merge them, otherwise if the
					# location is by itself, only add it if it has not been allocated
					# to another group. This prevents duplicate entries of the whole
					# region by itself plus any overlapping regions
					if len(grouped_locations) > 1 or location1 not in allocated:
							overlap_location = str(loc1)	# f'{chr1}:{start1}-{end1}'

							for location in grouped_locations:
									# id is a sample id
									sid = location_id_map[location]

									# sys.stderr.write('overlap ' + overlap_location + ' ' + id + ' ' + location + '\n')

									# .add(location)
									location_core_map[overlap_location][sid] = location

									used.add(location)
									allocated.add(location)

					if len(grouped_locations) == 1:
							# no more to add so quit looping
							exhausted = True

	# after iterating over everything, group locations by group

	return location_core_map


def overlapping_peaks(sids: list[tuple[str, str]]) -> dict[str, dict[str, str]]:
		"""
		Takes a list of file and finds the common (if any) overlapping regions

		Args:
			sids (list[tuple[str, str]]):	 list of sample id and file to match

		Returns:
			dict[str, dict[str, str]]: _description_
		"""

		location_id_map = collections.defaultdict(str)
		location_map = collections.defaultdict(int)
		location_bins = collections.defaultdict(set)
		locations = []

		for item in sids:
			sid, file = item

			f = open(file, 'r')

			# Skip header
			if 'Peaks' in file:
				f.readline()

			for line in f:
				tokens = line.strip().split('\t')

				if genomic.is_location(tokens[0]):
					location = genomic.parse_location(tokens[0])
				else:
						if genomic.is_chr(tokens[0]):
							location = genomic.Location(
									tokens[0], int(tokens[1]), int(tokens[2]))
						else:
							print(f'Invalid line: {line}', file=sys.stderr)
							continue

				lid = f'{sid}={str(location)}'

				# sys.stderr.write('lid ' + lid + '\n')

				locations.append(lid)

				# mapping from location id to sample id
				location_id_map[lid] = sid

				location_map[lid] = location

				bin_start = int(location.start / BIN_SIZE)
				bin_end = int(location.end / BIN_SIZE)

				for bin in range(bin_start, bin_end + 1):
					location_bins[bin].add(lid)

		f.close()

		return _overlapping(location_id_map, location_map, location_bins, locations)


def overlapping_peak_tables(files, ids):
		"""
		Takes a list of file and finds the common (if any) overlapping regions
		"""

		location_id_map = collections.defaultdict(str)
		location_map = collections.defaultdict(int)
		location_bins = collections.defaultdict(set)
		locations = []

		for i in range(0, len(files)):
			file = files[i]
			id = ids[i]

			f = open(file, 'r')

			f.readline()

			for line in f:
				ls = line.strip()

				if len(ls) == 0:
						continue

				tokens = line.split('\t')

				location = genomic.parse_location(tokens[0])

				lid = f'{id}={location.chr}:{location.start}-{location.end}'

				locations.append(lid)

				# sys.stderr.write(id + '\n')
				location_id_map[lid] = id

				location_map[lid] = location

				bin_start = int(location.start / BIN_SIZE)
				bin_end = int(location.end / BIN_SIZE)

				for bin in range(bin_start, bin_end + 1):
					location_bins[bin].add(lid)

		f.close()

		return _overlapping(location_id_map, location_map, location_bins, locations)
