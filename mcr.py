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
# from operator import itemgetter
import os
import sys
import pandas as pd
import numpy as np
import re

from . import genomic
from . import text

BIN_SIZE = 1000


def get_uid(sid: str, loc: genomic.Location) -> str:
    return f'{sid}={loc}'


def parse_uid(uid: str) -> list[str, str]:
    """Convert uid to sample id and location

    Args:
        uid (str): A uid to split

    Returns:
        list[str, str]: the sample id and location as strings
    """
    return uid.split('=')

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
        print(f"Loading peaks from {file}...", file=sys.stderr)

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


def get_test_uids(uid1: str, loc1: genomic.Location, bin_to_uids_map: dict[str, dict[int, list[str]]]) -> set[str]:
    sid1, _ = parse_uid(uid1)

    test_uids = set()

    bin_start = int(loc1.start / BIN_SIZE)
    bin_end = int(loc1.end / BIN_SIZE)

    for bin in range(bin_start, bin_end + 1):
        for uid in bin_to_uids_map[loc1.chr][bin]:
            sid2, _ = parse_uid(uid)
            if sid2 != sid1:
                # only test samples from other files, note that
                # we allow peaks from the allocated list since
                # the peak we are testing might overlap a peak
                # that was already allocated. The allocated peak
                # won't be check again, but we will check the
                # contrary
                test_uids.add(uid)

    return test_uids


def _min_common_regions(uids: list[str],
                        uid_to_loc_map: dict[str, genomic.Location],
                        bin_to_uids_map: dict[str, dict[int, list[str]]]) -> dict[str, dict[str, set[str]]]:
    """
    Calculates all overlaps between a set of sample locations

    Args:
        uids (list[str]): all uids across all samples
        uid_to_loc_map (dict[str, genomic.Location]): uid to genomic location
        bin_to_uids_map (dict[int, list[str]]): location bins to uids

    Returns:
        dict[str, dict[str, str]]: _description_
    """

    # lets see what overlaps

    location_core_map = collections.defaultdict(
        lambda: collections.defaultdict(set))

    # debug for testing to end remove as it truncates list
    # locations = locations[1:(len(locations) / 4)]

    total = len(uids)

    total *= total

    print(f'Processing {len(uids)}...', file=sys.stderr)

    p = 1

    # keep track of all locations that have been allocated at least once
    allocated = set()

    # test all locations in first sample
    for uid1 in uids:
        # of the form id=chrN:start-end, an lid

        # if the peaks has already been assigned to
        # an overlap, don't bother checking it as it
        # will duplicate regions
        if uid1 in allocated:
            continue

        # get its location
        loc1 = uid_to_loc_map[uid1]

        # find possible overlapping locations

        # test_uids = set()

        # bin_start = int(loc1.start / BIN_SIZE)
        # bin_end = int(loc1.end / BIN_SIZE)

        # for bin in range(bin_start, bin_end + 1):
        #     for uid in bin_to_uids_map[bin]:
        #         sid2, _ = parse_uid(uid1)
        #         if sid2 != sid1:
        #             # only test samples from other files, note that
        #             # we allow peaks from the allocated list since
        #             # the peak we are testing might overlap a peak
        #             # that was already allocated. The allocated peak
        #             # won't be check again, but we will check the
        #             # contrary
        #             test_uids.add(uid)

        test_uids = get_test_uids(uid1, loc1, bin_to_uids_map)

        used = set()

        # Form the largest group of overlapping peaks
        exhausted = False

        if 'chr10:110153776-110153972' in loc1.__str__():
            print('ww')

        while not exhausted:
            grouped_uids = {uid1}

            # reset for each group search
            loc1 = uid_to_loc_map[uid1]

            for uid2 in test_uids:
                # if uid1 == uid2:
                # continue

                if uid2 in used:
                    continue

                # sys.stderr.write(uid1 + ' ' + uid2 + '\n')

                loc2 = uid_to_loc_map[uid2]

                # Given the two locations we are testing, sort them so
                # the starts and ends are in order.
                # if (loc1.start <= loc2.start):
                #     min_loc = loc1
                #     max_loc = loc2
                # else:
                #     min_loc = loc2
                #     max_loc = loc1

                # overlap = -1
                # overlap_start = -1

                # if max_loc.start <= min_loc.end:
                #     overlap_start = max_loc.start
                #     overlap = min(min_loc.end, max_loc.end) - overlap_start + 1

                overlap = genomic.overlap_locations(loc1, loc2)
                
                if 'chr10:110153776-110153972' in loc1.__str__() or 'chr10:110153776-110153972' in loc2.__str__():
                    print('overlap', loc1, loc2, overlap)

                if overlap is not None:
                    # change the start1 and end1 coordinates to reflect the overlap
                    # regions so that each subsequent match must be within this region
                    # this prevents long peaks that overlap two smaller peaks who
                    # themselves do not overlap each other
                    # genomic.Location(loc1.chr, overlap_start, overlap_start + overlap - 1)
                    loc1 = overlap

                    # we found someone we are overlapping
                    grouped_uids.add(uid2)

                # if max_loc.start == min_loc.start and max_loc.end > min_loc.end:
                # 	# Peak 1 and 2 start at the same point but peak 2 is wider
                # 	# so the overlap region is peak 1
                # 	overlap = min_loc.end - min_loc.start + 1
                # 	overlap_start = min_loc.start
                # elif max_loc.start >= min_loc.start and max_loc.end <= min_loc.end:
                # 	# Peak 1 is wider than peak 2 and contains it
                # 	# so the overlap region is peak 2
                # 	overlap = max_loc.end - max_loc.start + 1
                # 	overlap_start = max_loc.start
                # elif min_loc.start < max_loc.start and min_loc.end > max_loc.start:
                # 	# Peak 1 starts before peak 2 but ends within peak 2
                # 	# so the overlap is the start of peak 2 to the end of
                # 	# peak 1
                # 	overlap = min_loc.end - max_loc.start + 1
                # 	overlap_start = max_loc.start
                # else:
                # 	pass

                # We have not found an overlap yet so continue
                # if overlap == -1:
                #    continue

                # change the start1 and end1 coordinates to reflect the overlap
                # regions so that each subsequent match must be within this region
                # this prevents long peaks that overlap two smaller peaks who
                # themselves do not overlap each other
                # loc1 = overlap #genomic.Location(loc1.chr, overlap_start, overlap_start + overlap - 1)

                # we found someone we are overlapping
                # grouped_uids.add(uid2)

            # now we have a list of all locations that overlap each other

            # if we have a group of entries, merge them, otherwise if the
            # location is by itself, only add it if it has not been allocated
            # to another group. This prevents duplicate entries of the whole
            # regions by itself plus any overlapping regions
            # if len(grouped_locations) > 1 or uid1 not in allocated:

            overlap_location = str(loc1)  # f'{chr1}:{start1}-{end1}'

            # if (location_map[uid1].start == 100008783):
            #	print(loc1, loc2, max_loc.start <= min_loc.end, overlap, overlap_start, grouped_locations)
            # exit(0)

            # if there are multiple locations, group them and mark as
            # allocated. If there is only one location, it either means
            # we found nothing or the loops are exhausted. We check that
            # the grouped location (same as uid1) has not been allocated
            # in a previous loop and if not, it means this is a single
            # location and does not overlap anything so we can mark it
            # allocated and move on

            if len(grouped_uids) > 1 or next(iter(grouped_uids)) not in allocated:
                
            
                for uid in grouped_uids:
                    # sid is a sample id
                    sid, _ = parse_uid(uid)  # sample_id_map[uid]

                    # if the uid has not been allocated yet
                    location_core_map[overlap_location][sid].add(uid)

                    if 'chr10:110153776-110153972' in uid:
                        print('allocated 1', uid1, uid)

                    if 'chr10:110153588-110153799' in uid:
                        print('allocated 2', uid1, uid)

                    used.add(uid)
                    allocated.add(uid)

            if len(grouped_uids) == 1:
                # we can stop looking
                exhausted = True

    # after iterating over everything, group locations by group

    return location_core_map


def min_common_regions(fids: list[tuple[str, str]], core_regions=_min_common_regions) -> tuple[dict[str, dict[str, str]], dict[str, genomic.Location]]:
    """
    Takes a list of file and finds the common (if any) overlapping regions

    Args:
            fids (list[tuple[str, str]]):	 list of sample id and file to match

    Returns:
            dict[str, dict[str, str]]: _description_
    """
    
    loc_sample_map = collections.defaultdict(list[str])
    locations = []

    for item in fids:
        sid, file = item

        with open(file, 'r') as f:

            # Skip header
            if 'Peaks' in file or 'tsv' in file:
                f.readline()

            for line in f:
                tokens = line.strip().split('\t')

                if genomic.is_location(tokens[0]):
                    location = genomic.parse_location(tokens[0])
                else:
                    if genomic.is_chr(tokens[0]):
                        location = genomic.Location(
                            tokens[0], int(tokens[1]), int(tokens[2]))
                        locations.append(location)
                        loc_sample_map[str(location)].append(sid)
                    else:
                        print(f'Invalid line: {line}', file=sys.stderr)
    
    # sort all locations
    locations = genomic.sort_locations(locations)

    # map sorted locations back to sample and run in genomic location
    # order

    uids = []
    location_map = collections.defaultdict(genomic.Location)
    bin_to_uids_map = collections.defaultdict(lambda: collections.defaultdict(list[str]))

    for location in locations:
        sids = loc_sample_map[str(location)]

        for sid in sids:
            uid = get_uid(sid, location)
            
            uids.append(uid)

            location_map[uid] = location

            bin_start = int(location.start / BIN_SIZE)
            bin_end = int(location.end / BIN_SIZE)

            for bin in range(bin_start, bin_end + 1):
                bin_to_uids_map[location.chr][bin].append(uid)

    location_core_map = core_regions(uids, location_map, bin_to_uids_map)

    return location_core_map, location_map


def create_overlap_table(files: list[str], core_regions=min_common_regions):
    sids = []
    fids = []

    # total_score_map = collections.defaultdict(float)
    # max_score_map = collections.defaultdict(float)

    ext_cols = []
    ext_data = collections.defaultdict(lambda: collections.defaultdict(float))

    for file in files:
        print(f'file: {file}', file=sys.stderr)

        if ',' in file:
            sid, file = file.split(',')
        else:
            # sid = os.path.basename(os.path.dirname(
            # os.path.dirname(os.path.dirname(file))))

            sid = re.sub(r'\.[^\.]+$', '', os.path.basename(file))

        sids.append(sid)
        fids.append([sid, file])

        #print(f'ids {sid}', file=sys.stderr)

        # now make a list of locations and best p-values

        f = open(file, 'r')

        # total_score_col = -1  # 3
        # max_score_col = -1  # 4

        ext_col_indexes = {}

        # Adjust colums to look it for peak files
        if "seacr" in file and 'tsv' in file:
            # print(f'what gives', file=sys.stderr)
            tokens = f.readline().strip().split("\t")

            # total_score_col = gal.text.find_index(tokens, "Total Score")
            # max_score_col = gal.text.find_index(tokens, "Max Score")

            ext_cols = ["Total Score", "Max Score"]
            ext_col_indexes['Total Score'] = text.find_index(
                tokens, "Total Score")
            ext_col_indexes['Max Score'] = text.find_index(
                tokens, "Max Score")
        elif ("narrowPeak" in file or "broadPeak" in file) and not file.endswith('bed'):
            ext_cols = ["fold_change", "-log10pvalue", "-log10qvalue"]
            ext_col_indexes['fold_change'] = 6
            ext_col_indexes["-log10pvalue"] = 7
            ext_col_indexes["-log10qvalue"] = 8
        elif file.endswith('bed'):
            pass
        else:
            # assume table so skip first line
            f.readline()

        for line in f:
            tokens = line.strip().split("\t")

            if genomic.is_location(tokens[0]):
                location = genomic.parse_location(tokens[0])
            else:
                if genomic.is_chr(tokens[0]):
                    location = genomic.Location(
                        tokens[0], int(tokens[1]), int(tokens[2]))
                else:
                    print(f'Invalid line: {line}', file=sys.stderr)

                    continue

            # total_score = 0
            # max_score = 0

            # if total_score_col != -1:
            #	total_score = float(tokens[total_score_col])

            # if max_score_col != -1:
            #	max_score = float(tokens[max_score_col])

            uid = get_uid(sid, location)

            

            for col in ext_cols:
                ext_data[uid][col] = float(tokens[ext_col_indexes[col]])

            # total_score_map[lid] = total_score
            # max_score_map[lid] = max_score

            # sys.stderr.write(str(location) + " " + lid + "\n")
            # exit(0)

        f.close()

    location_core_map, location_map = core_regions(
        fids)  # min_common_regions(fids)

    # keep the ids sorted
    # ids = sorted(ids)

    header = ["Genomic Location", "Width", "Region", "Region Width"]

    for sid in sids:
        header.extend([f'{sid} {c}' for c in ext_cols])

    header.append("# Overlapping Samples")
    header.append("# Overlapping Peaks")

    header.extend([f'Sample {s}' for s in sids])

    #header.extend([f'Peak {s}' for s in sids])

    header.extend([f'Overlap % {s}' for s in sids])

    # print("\t".join(header))

    data = []

    for core_location in sorted(location_core_map):
        overlap_location = genomic.parse_location(core_location)

        overlap_width = overlap_location.end - overlap_location.start + 1

        sample_count = len(location_core_map[core_location])

        peak_count = 0
        for sid in location_core_map[core_location]:
            peak_count += len(location_core_map[core_location][sid])

        locs: list[genomic.Location] = []

        for sid in sids:
            if sid in location_core_map[core_location]:
                uid = location_core_map[core_location][sid]
                locs.extend(genomic.sort_locations(
                    [location_map[uid] for uid in location_core_map[core_location][sid]]))

        start = min([loc.start for loc in locs])
        end = max([loc.end for loc in locs])
        region = genomic.Location(locs[0].chr, start, end)
        region_width = region.end - region.start + 1

        row = [core_location, str(overlap_width), region, str(region_width)]

        # for id in location_core_map[core_location]:
        # for location in location_core_map[core_location][id]:
        #	c += 1

        # total_score = -1

        # max_score = -1

        for sid in sids:
            if sid in location_core_map[core_location]:
                # for location in location_core_map[core_location][id]:
                uid = location_core_map[core_location][sid]
                row.extend([str(ext_data[uid][col]) for col in ext_cols])
            else:
                row.extend([text.NA] * len(ext_cols))

            # row.extend([str(ext_data[lid][col]) for col in ext_cols])

        # add count/number of cols we appear in
        row.append(str(sample_count))
        row.append(str(peak_count))

        # place ids in table skip and just use peaks locations as peak ids
        # for sid in sids:
        # 	if sid in location_core_map[core_location]:
        # 		uid = location_core_map[core_location][sid]

        # 		row.append(uid)
        # 	else:
        # 		row.append(text.NA)

        for sid in sids:
            if sid in location_core_map[core_location]:
                # ";".join(sorted(location_core_map[core_location][id]))
                #uid = location_core_map[core_location][sid]

                locs = genomic.sort_locations(
                    [location_map[uid] for uid in sorted(location_core_map[core_location][sid])])

                row.append(';'.join([str(loc) for loc in locs]))
            else:
                row.append(text.NA)

            # % overlap

        for sid in sids:
            if sid in location_core_map[core_location]:
                # ";".join(sorted(location_core_map[core_location][id]))
                #uid = location_core_map[core_location][sid]

                #loc1 = location_map[uid]

                locs = genomic.sort_locations(
                    [location_map[uid] for uid in sorted(location_core_map[core_location][sid])])

                fracs = [genomic.overlap_fraction(
                    overlap_location, loc) for loc in locs]

                # overlap_fraction = genomic.overlap_fraction(
                #    overlap_location, loc1)

                row.append(
                    ';'.join([str(min(100, max(0, np.round(f * 100, 2)))) for f in fracs]))
            else:
                row.append('0')

        # print("\t".join(row))

        # find the spanning region of the two peaks
        # i.e. min start, max end

        data.append(row)

    return pd.DataFrame(data, columns=header)
