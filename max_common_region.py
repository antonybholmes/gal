import collections
import os
import re
import sys
import pandas as pd
import numpy as np
from . import genomic
from . import text
from . import mcr


def _max_common_regions(uids: list[str],
                        uid_to_loc_map: dict[str, genomic.Location],
                        bin_to_uids_map: dict[int, list[str]]) -> dict[str, dict[str, set[str]]]:
    """
    Join all overlapping peaks into one, i.e. keep expanding the region as we
    find overlapping peaks. Unlike finding the minimum

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

        sid1, _ = mcr.parse_uid(uid1)

        # get its location
        loc1 = uid_to_loc_map[uid1]

        # find possible overlapping locations
        test_uids = mcr.get_test_uids(uid1, loc1, bin_to_uids_map)

        used = set()

        # Form the largest group of overlapping peaks
        exhausted = False

        while not exhausted:
            grouped_locations = {uid1}

            # reset for each group search
            loc1 = uid_to_loc_map[uid1]

            for uid2 in test_uids:
                if uid2 in used:
                    continue

                loc2 = uid_to_loc_map[uid2]

                overlap = genomic.overlap_locations(loc1, loc2)

                if overlap is not None:
                    # expand search to region covered by both peaks
                    grouped_locations.add(uid2)

                    loc1 = genomic.max_region(loc1, loc2)
                    test_uids = mcr.get_test_uids(uid1, loc1, bin_to_uids_map)

            # now we have a list of all locations that overlap each other

            # if we have a group of entries, merge them, otherwise if the
            # location is by itself, only add it if it has not been allocated
            # to another group. This prevents duplicate entries of the whole
            # regions by itself plus any overlapping regions
            # if len(grouped_locations) > 1 or uid1 not in allocated:

            overlap_location = str(loc1)  # f'{chr1}:{start1}-{end1}'
            
            if 'chr10:104809368' in str(loc1):
                print(loc1, len(grouped_locations), grouped_locations, len(loc1))

            for uid in grouped_locations:
                # sid is a sample id
                sid, _ = mcr.parse_uid(uid)  # sample_id_map[uid]

                # if the uid has not been allocated yet
                location_core_map[overlap_location][sid].add(uid)

                used.add(uid)
                allocated.add(uid)

            if len(grouped_locations) == 1:
                # we can stop looking
                exhausted = True

    # after iterating over everything, group locations by group

    return location_core_map


def max_common_regions(fids: list[tuple[str, str]]) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    return mcr.min_common_regions(fids, _max_common_regions)


def create_max_common_region_table(files: list[str]):
    return mcr.create_overlap_table(files, core_regions = max_common_regions)