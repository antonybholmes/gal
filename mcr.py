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
from functools import reduce, total_ordering

# from operator import itemgetter
import os
from pathlib import Path
import re
import sys
from typing import Callable, Iterable, Union

import numpy as np
import pandas as pd

from . import genomic, text

BIN_SIZE = 1000


# def create_uid(sid: str, loc: genomic.Location) -> str:
#     """
#     Create a unique identifier for a location in a given sample
#     """
#     return f"{sid}={loc}"


# def parse_uid(uid: str) -> list[str, str]:
#     """Convert uid to sample id and location

#     Args:
#         uid (str): A uid to split

#     Returns:
#         list[str, str]: the sample id and location as strings
#     """
#     return uid.split("=")


@total_ordering
class ULoc(genomic.Location):
    def __init__(self, id: str, chr: str, start: int, end: int, strand: str = "+"):
        super().__init__(chr, start, end, strand)
        self._id = id

    @property
    def id(self) -> str:
        """
        The sample id
        """
        return self._id

    def __str__(self) -> str:
        return f"{self._id}={str(super().__str__())}"

    def __hash__(self) -> int:
        return hash(self.__str__())

    def __eq__(self, other) -> bool:
        if not isinstance(other, ULoc):
            return False

        return self._id == other.id and super().__eq__(other)

    def __lt__(self, other) -> bool:
        if not isinstance(other, ULoc):
            return True

        if self._id != other.id:
            return self._id < other.id

        return super().__lt__(other)

    def to_loc(self) -> genomic.Location:
        """
        Turns the UID back to a simple genomic location. Useful in
        instances where you want to compare the location portion of
        the uid without considering the sid.
        """
        return genomic.Location(self.chr, self.start, self.end, self.strand)

    @classmethod
    def create(cls, loc: genomic.Location, id: str = "uid") -> "ULoc":
        return cls(id, loc.chr, loc.start, loc.end, loc.strand)

    @classmethod
    def from_list(
        cls, locs: Iterable[genomic.Location], id: str = "uid"
    ) -> list["ULoc"]:
        return [ULoc.create(loc) for loc in locs]


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

    def __init__(self, file: str, bin_size: int = BIN_SIZE):
        self._bins = collections.defaultdict(lambda: collections.defaultdict(set))
        self._bin_size = bin_size
        self.load(file)

    def load(self, file):
        print(f"Loading peaks from {file}...", file=sys.stderr)

        f = open(file, "r")

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split("\t")

            location = genomic.parse_location(tokens[0])

            sbin = int(location.start / self._bin_size)
            ebin = int(location.end / self._bin_size)

            for i in range(sbin, ebin + 1):
                self._bins[location.chr][i].add(location)

        f.close()

    def overlaps(self, location):
        sbin = int(location.start / self._bin_size)
        ebin = int(location.end / self._bin_size)

        overlaps = []

        for i in range(sbin, ebin + 1):
            for l in self._bins[location.chr][i]:
                if genomic.is_overlapping(location, l):
                    overlaps.append(l)

        return overlaps


def overlap(ref_file, query_file):
    peaks = peaks.parse_peaks(ref_file)

    f = open(query_file, "r")

    # skip header
    # f.readline()

    lines = []

    for line in f:
        tokens = line.split("\t")

        chr = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])

        features = peak_overlap(peaks, chr, start, end)

        lines.append(
            "\t".join([f"{chr}:{start}-{end}", str(len(features)), ";".join(features)])
        )

    f.close()

    return lines


def peak_overlap(peaks, chr, start, end):
    features = []

    for id in sorted(peaks[chr]):
        p_start = int(peaks[chr][id][0])
        p_end = int(peaks[chr][id][1])

        if start >= p_start and end <= p_end:
            features.append(";".join([id, "within"]))
        elif start < p_start and end > p_end:
            features.append(";".join([id, "over"]))
        elif start < p_start and end > p_start:
            features.append(";".join([id, "upstream"]))
        elif start < p_end and end > p_end:
            features.append(";".join([id, "downstream"]))

    return features


def get_test_ulocs(
    uloc1: ULoc,
    loc: genomic.Location,
    bin_to_uloc_map: dict[str, dict[int, list[str]]],
    bin_size: int = BIN_SIZE,
) -> set[str]:
    """
    Since we don't know what overlaps, we bin all locations into fixed bp gaps and then
    each time we need to test a location, we get all other sequences in the same bin as
    the test location and test those instead of testing all sequences since only the
    sequences in the same bins have any chance of overlapping our location. Note that
    locations returned here may not overlap the location, they are just close to it
    so further overlap tests are required. This is just an initial filter to reduce
    the search space. This set will not contain the input location, hence it is all
    locations close to uid1, but not including it (since there is no reason to test
    uid1 against itself).

    args:
        loc -   location to test, this can be independent of the uid, for example
                in the min/max where we are scanning a larger region than the uid
    """

    # sid1, _ = parse_uid(uid1)

    test_uids: list[ULoc] = []
    used: set[ULoc] = set()

    bin_start = int(loc.start / bin_size)
    bin_end = int(loc.end / bin_size)

    for bin in range(bin_start, bin_end + 1):
        for uid2 in bin_to_uloc_map[loc.chr][bin]:
            # sid2, _ = parse_uid(uid)
            # if sid2 != sid1:
            if uid2 != uloc1 and uid2 not in used:
                # only test samples from other files, note that
                # we allow peaks from the allocated list since
                # the peak we are testing might overlap a peak
                # that was already allocated. The allocated peak
                # won't be check again, but we will check the
                # contrary
                test_uids.append(uid2)
                used.add(uid2)

    return test_uids


def _min_common_regions(
    ulocs: list[ULoc],
    bin_to_uloc_map: dict[str, dict[int, list[ULoc]]],
    bin_size: int = BIN_SIZE,
) -> dict[genomic.Location, dict[str, list[ULoc]]]:
    """
    Calculates minimal common regions between a set of sample locations

    Args:
        uids (list[str]): all uids across all samples
        bin_to_uids_map (dict[int, list[str]]): location bins to uids
        bin_size - size of bins to partition genome into

    Returns:
        dict[genomic.Location, dict[str, list[UID]]]: a map of locations to sets of uids grouped by sample.
    """

    # lets see what overlaps

    location_core_map: dict[genomic.Location, dict[str, list[ULoc]]] = (
        collections.defaultdict(lambda: collections.defaultdict(list))
    )

    # debug for testing to end remove as it truncates list
    # locations = locations[1:(len(locations) / 4)]

    total = len(ulocs)

    total *= total

    print(f"Processing {len(ulocs)} regions...", file=sys.stderr)

    p = 1

    # keep track of all locations that have been allocated at least once
    allocated = set()

    # test all locations in first sample
    for uloc1 in ulocs:
        # of the form id=chrN:start-end, an lid

        # a new overlap cannot start at a location that has
        # already been assigned to cluster, but clusters
        # can be formed from locations already assigned since
        # that means a location can be split into multiple pieces.
        if uloc1 in allocated:
            continue

        # print(uid1, loc1)

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

        test_uids = get_test_ulocs(
            uloc1=uloc1, loc=uloc1, bin_to_uloc_map=bin_to_uloc_map, bin_size=bin_size
        )

        # print(loc1, test_uids)

        #
        # simplified
        #

        grouped_uids = [uloc1]

        # reset for each group search
        loc1 = uloc1.to_loc()

        for uid2 in test_uids:
            # test ids are everything but uid1

            overlap = genomic.overlap_locations(loc1, uid2)

            if overlap is not None:
                # reduce size of overlap each time we overlap something
                loc1 = overlap
                # we found someone we are overlapping
                grouped_uids.append(uid2)

        # if there are multiple locations, group them and mark as
        # allocated. If there is only one location, it either means
        # we found nothing or the loops are exhausted.
        for uid in grouped_uids:
            # sid is a sample id
            # sid, _ = parse_uid(uid)
            # if the uid has not been allocated yet
            location_core_map[loc1][uid.id].append(uid)
            allocated.add(uid)

        #
        # end simplified
        #

        #
        # long winded can we get rid of it?
        #

        # used = {uid1}

        # # Form the largest group of overlapping peaks
        # exhausted = False

        # while not exhausted:
        #     # seed the group with first location
        #     grouped_uids = [uid1]

        #     # reset for each group search
        #     loc1 = uid1.loc

        #     for uid2 in test_uids:
        #         if uid2 in used:
        #             continue

        #         overlap = genomic.overlap_locations(loc1, uid2.loc)

        #         # print(loc1, loc2, overlap, file=sys.stderr)

        #         if overlap is not None:
        #             # change the start1 and end1 coordinates to reflect the overlap
        #             # regions so that each subsequent match must be within this region
        #             # this prevents long peaks that overlap two smaller peaks who
        #             # themselves do not overlap each other
        #             # genomic.Location(loc1.chr, overlap_start, overlap_start + overlap - 1)
        #             loc1 = overlap

        #             # we found someone we are overlapping
        #             grouped_uids.append(uid2)
        #             # used.add(uid2)

        #     # if there are multiple locations, group them and mark as
        #     # allocated. If there is only one location, it either means
        #     # we found nothing or the loops are exhausted.

        #     if len(grouped_uids) > 1:
        #         for uid in grouped_uids:
        #             # sid is a sample id
        #             # sid, _ = parse_uid(uid)  # sample_id_map[uid]

        #             # if the uid has not been allocated yet
        #             location_core_map[loc1][uid.id].add(uid)

        #             used.add(uid)
        #             allocated.add(uid)

        #             # print(uid1)

        #             if "136874906" in str(uid1):
        #                 print(uid, file=sys.stderr)
        #     else:
        #         # group contains 1 item, the seed uid
        #         # We check that the grouped location (same as uid1) has not been
        #         # allocated in a previous loop and if not, it means this is a single
        #         # location and does not overlap anything so we can mark it
        #         # allocated and move on
        #         if uid1 not in allocated:
        #             location_core_map[loc1][uid1.sid].add(uid1)
        #             allocated.add(uid1)

        #         # we can stop looking
        #         exhausted = True
        #         break

        # if "136874906" in str(uid1):
        #     print(location_core_map[str(loc1)][uid1.sid], file=sys.stderr)

        #
        # end long winded
        #

    # after iterating over everything, group locations by group

     
    # exit(0)

    return location_core_map


def make_bin_to_uloc_map(
    ulocs: list[ULoc], bin_size: int = BIN_SIZE
) -> dict[str, dict[int, list[ULoc]]]:
    """Converts a list of uids into a searchable map

    Args:
        uids (list[UID]): A sorted (by location) list of UIDs
        bin_size (int, optional): _description_. Defaults to BIN_SIZE.

    Returns:
        dict[str, dict[int, list[UID]]]: a chr->bin->uids map
    """

    # uids = list(sorted(uids))

    bin_to_ulocs_map = collections.defaultdict(
        lambda: collections.defaultdict(list[ULoc])
    )

    for uloc in ulocs:
        bin_start = int(uloc.start / bin_size)
        bin_end = int(uloc.end / bin_size)

        for bin in range(bin_start, bin_end + 1):
            bin_to_ulocs_map[uloc.chr][bin].append(uloc)

    return bin_to_ulocs_map


def group_ulocs_by_core_location(
    ulocs: list[ULoc],
    func_core_regions: Callable[
        [list[ULoc], dict[str, dict[int, list[ULoc]]], int],
        dict[genomic.Location, dict[str, list[ULoc]]],
    ] = _min_common_regions,
    bin_size: int = BIN_SIZE,
) -> dict[genomic.Location, dict[str, list[ULoc]]]:
    """_summary_

    Args:
        uids (list[UID]): list of UIDS to group
        func_core_regions (Callable[ [list[UID], dict[str, dict[int, list[UID]]], int], dict[genomic.Location, dict[str, list[UID]]], ], optional): A function to do the actual grouping. Defaults to _min_common_regions.
        bin_size (int, optional): Size of genomic bins for searching. Smaller reduces search space and speeds up comparisons. Defaults to BIN_SIZE.

    Returns:
        dict[genomic.Location, dict[str, list[UID]]]: _description_
    """

    if not all(isinstance(uloc, ULoc) for uloc in ulocs):
        raise TypeError("ulocs must be a list of type ULoc.")

    # ensure all locations are sorted
    ulocs = list(sorted(ulocs))  # genomic.sort_locations(locations)

    # map sorted locations back to sample and run in genomic location
    # order

    # uid_map = collections.defaultdict(list)
    # uid_to_loc_map = collections.defaultdict(genomic.Location)
    bin_to_uloc_map = make_bin_to_uloc_map(ulocs, bin_size=bin_size)

    # lets keep uids sorted by sample and coordinate
    # use reduce to flatten list of lists
    # uids: list[str] = list(
    #     reduce(
    #         lambda x, y: x + y,
    #         [
    #             [create_uid(sid, loc) for loc in sorted(uid_map[sid])]
    #             for sid in sorted(uid_map)
    #         ],
    #     )
    # )

    # print(uids)

    location_core_map = func_core_regions(ulocs, bin_to_uloc_map, bin_size)

    return location_core_map


def files_to_ulocs(
    fids: list[tuple[str, str]],
) -> list[ULoc]:
    """Takes a list of file and finds the common (if any) overlapping regions

    Args:
        fids (list[tuple[str, str]]): list of sample id and file to match

    Returns:
        list[ULoc]: A list of uids generated from the files
    """

    # loc_sample_map = collections.defaultdict(list[str])
    ulocs: list[ULoc] = []

    for item in fids:
        sid, file = item

        with open(file, "r") as f:
            # Skip header
            if "Peaks" in file or "tsv" in file or "txt" in file:
                print(
                    f"{file} appears to have a header line, skipping...",
                    file=sys.stderr,
                )
                f.readline()

            for line in f:
                tokens = line.strip().split("\t")

                location = None

                if genomic.is_location(tokens[0]):
                    location = genomic.parse_location(tokens[0])
                else:
                    if genomic.is_chr(tokens[0]):
                        location = genomic.Location(
                            tokens[0], int(tokens[1]), int(tokens[2])
                        )
                    else:
                        print(f"{line} is an invalid location.", file=sys.stderr)

                if location is not None:
                    ulocs.append(ULoc.create(location, sid))

    return list(sorted(ulocs))


def common_regions(
    fids: list[tuple[str, str]],
    func_core_regions: Callable[
        [list[ULoc], dict[str, dict[int, list[ULoc]]], int],
        dict[genomic.Location, dict[str, list[ULoc]]],
    ],
    bin_size: int = BIN_SIZE,
) -> dict[str, dict[str, str]]:
    """
    Takes a list of file and finds the common (if any) overlapping regions

    Args:
            fids (list[tuple[str, str]]):	 list of sample id and file to match

    Returns:
            dict[str, dict[str, str]]: _description_
    """

    # loc_sample_map = collections.defaultdict(list[str])
    ulocs: list[ULoc] = files_to_ulocs(fids)

    group_ulocs_by_core_location(ulocs, func_core_regions, bin_size)


def min_common_regions(
    uids: list[ULoc], bin_size: int = BIN_SIZE
) -> dict[genomic.Location, dict[str, list[ULoc]]]:
    print("MCR mode.", file=sys.stderr)
    return group_ulocs_by_core_location(
        uids, func_core_regions=_min_common_regions, bin_size=bin_size
    )


def create_mcr_table(
    files: list[str],
    func_core_regions: Callable[
        [list[ULoc], int],
        dict[genomic.Location, dict[str, list[ULoc]]],
    ] = min_common_regions,
    bin_size: int = BIN_SIZE,
):
    return create_overlap_table(
        files, func_core_regions=func_core_regions, bin_size=bin_size
    )


def create_overlap_table(
    files: list[str],
    func_core_regions: Callable[
        [list[ULoc], int],
        dict[genomic.Location, dict[str, list[ULoc]]],
    ] = min_common_regions,
    bin_size: int = BIN_SIZE,
) -> pd.DataFrame:
    """_summary_

    Args:
        files (list[str]): A list of files containing coordinates to overlap.
        func_core_regions (Callable[ [list[UID], int], dict[genomic.Location, dict[str, list[UID]]], ], optional): A function to do the overlap either MCR or maximal. Defaults to min_common_regions.
        bin_size (int, optional): _description_. Defaults to BIN_SIZE.

    Returns:
        pd.DataFrame: A dataframe highlighting all the overlaps
    """

    sids: list[str] = []
    fids: list[tuple[str, str]] = []

    # total_score_map = collections.defaultdict(float)
    # max_score_map = collections.defaultdict(float)

    ext_cols = []
    ext_data: map[ULoc, map[str, Union[str, int, float]]] = collections.defaultdict(
        lambda: collections.defaultdict(float)
    )

    all_ulocs: list[ULoc] = []

    for file in files:
        print(f"file: {file}", file=sys.stderr)

        if "," in file:
            sid, file = file.split(",")[0:2]
        else:
            # sid = os.path.basename(os.path.dirname(
            # os.path.dirname(os.path.dirname(file))))

            # remove file extension to get sample/file id
            sid = Path(file).stem
            # , _ = os.path.splitext(os.path.basename(file))
            # re.sub(r"\.[^\.]+$", "", os.path.basename(file))

        sids.append(sid)
        fids.append((sid, file))

        # print(f'ids {sid}', file=sys.stderr)

        # now make a list of locations and best p-values

        with open(file, "r") as f:
            # total_score_col = -1  # 3
            # max_score_col = -1  # 4

            ext_col_indexes = {}

            # Adjust colums to look it for peak files
            if "seacr" in file and "tsv" in file:
                print(
                    f"{file} appears to come from SEACR, adding extra columns to output table...",
                    file=sys.stderr,
                )

                tokens = f.readline().strip().split("\t")

                # total_score_col = gal.text.find_index(tokens, "Total Score")
                # max_score_col = gal.text.find_index(tokens, "Max Score")

                ext_cols = ["Total Score", "Max Score"]
                ext_col_indexes["Total Score"] = text.find_index(tokens, "Total Score")
                ext_col_indexes["Max Score"] = text.find_index(tokens, "Max Score")

            elif ("narrowPeak" in file or "broadPeak" in file) and not file.endswith(
                "bed"
            ):
                print(
                    f"{file} appears to come from MACS, adding extra columns to output table...",
                    file=sys.stderr,
                )

                ext_cols = ["fold_change", "-log10pvalue", "-log10qvalue"]
                ext_col_indexes["fold_change"] = 6
                ext_col_indexes["-log10pvalue"] = 7
                ext_col_indexes["-log10qvalue"] = 8
            elif file.endswith("bed"):
                pass
            else:
                print(
                    f"{file} seems to be a text file, skipping header...",
                    file=sys.stderr,
                )
                # assume table so skip first line
                f.readline()

            for line in f:
                tokens = line.strip().split("\t")

                if genomic.is_location(tokens[0]):
                    location = genomic.parse_location(tokens[0])
                else:
                    # bed format which assumes start is zero based so
                    # need to convert it to one based
                    if genomic.is_chr(tokens[0]):
                        start = int(tokens[1])

                        if file.endswith("bed"):
                            start += 1

                        location = genomic.Location(
                            tokens[0], start, int(tokens[2])
                        )
                    else:
                        print(f"Invalid line: {line}", file=sys.stderr)

                        continue

                # total_score = 0
                # max_score = 0

                # if total_score_col != -1:
                # 	total_score = float(tokens[total_score_col])

                # if max_score_col != -1:
                # 	max_score = float(tokens[max_score_col])

                uloc = ULoc.create(location, sid)
                all_ulocs.append(uloc)

                for col in ext_cols:
                    ext_data[uloc][col] = float(tokens[ext_col_indexes[col]])

                # total_score_map[lid] = total_score
                # max_score_map[lid] = max_score

                # sys.stderr.write(str(location) + " " + lid + "\n")
                # exit(0)

    all_ulocs = list(sorted(all_ulocs))

    core_location_map = func_core_regions(
        all_ulocs, bin_size
    )  # min_common_regions(fids)

    # keep the ids sorted
    # ids = sorted(ids)

    header = ["Genomic Location", "Width", "Region", "Region Width"]

    for sid in sids:
        header.extend([f"{sid} {c}" for c in ext_cols])

    header.append("# Overlapping Samples")
    header.append("# Overlapping Peaks")

    header.extend([f"Sample {s}" for s in sids])

    # header.extend([f'Peak {s}' for s in sids])

    header.extend([f"Overlap % {s}" for s in sids])

    # print("\t".join(header))

    data = []

    for core_location in sorted(core_location_map):
        # overlap_location = genomic.parse_location(core_location)

        overlap_width = len(core_location)

        sample_count = len(core_location_map[core_location])

        peak_count = 0
        for sid in core_location_map[core_location]:
            peak_count += len(core_location_map[core_location][sid])

        core_ulocs: list[ULoc] = []

        # a uid is effectively a peak

        for sid in sids:
            # we check for key existence since we do not want
            # defaultdict to create an empty entry if sid does
            # not exist
            if sid in core_location_map[core_location]:
                core_ulocs.extend(core_location_map[core_location][sid])

        # sort locations
        # core_uids = list(sorted(core_uids))

        # the max region spanned by the overlapping locs
        #print(sids)

        start = min([uid.start for uid in core_ulocs])
        end = max([uid.end for uid in core_ulocs])
        region = genomic.Location(core_ulocs[0].chr, start, end)
        region_width = region.end - region.start + 1

        row = [str(core_location), str(overlap_width), region, str(region_width)]

        # for id in location_core_map[core_location]:
        # for location in location_core_map[core_location][id]:
        # 	c += 1

        # total_score = -1

        # max_score = -1

        for sid in sids:
            if sid in core_location_map[core_location]:
                # for location in location_core_map[core_location][id]:

                # core location will return a set of peaks for a sample,
                # in the mcr case this will be 1 sample, but just in case
                # sort by name and pick first
                #print(sid, core_location, core_location_map[core_location])

                uloc = core_location_map[core_location][sid][0]

                # print(ext_cols, uid)

                row.extend([str(ext_data[uloc][col]) for col in ext_cols])
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
            if sid in core_location_map[core_location]:
                # ";".join(sorted(location_core_map[core_location][id]))
                # uid = location_core_map[core_location][sid]

                core_ulocs: list[ULoc] = core_location_map[core_location][sid]

                row.append(";".join([str(uloc.to_loc()) for uloc in core_ulocs]))
            else:
                row.append(text.NA)

        # % overlap

        for sid in sids:
            if sid in core_location_map[core_location]:
                # ";".join(sorted(location_core_map[core_location][id]))
                # uid = location_core_map[core_location][sid]

                # loc1 = location_map[uid]

                core_ulocs: list[ULoc] = list(
                    sorted(core_location_map[core_location][sid])
                )

                fracs = [
                    genomic.overlap_fraction(core_location, uloc) for uloc in core_ulocs
                ]

                # overlap_fraction = genomic.overlap_fraction(
                #    overlap_location, loc1)

                row.append(
                    ";".join(
                        [str(min(100, max(0, np.round(f * 100, 2)))) for f in fracs]
                    )
                )
            else:
                row.append("0")

        # print("\t".join(row))

        # find the spanning region of the two peaks
        # i.e. min start, max end

        data.append(row)

    return pd.DataFrame(data, columns=header)
