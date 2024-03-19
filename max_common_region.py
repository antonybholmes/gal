import collections
import os
import re
import sys
import pandas as pd
import numpy as np
from . import genomic
from . import text
from . import mcr


def _max_common_regions(
    ulocs: list[mcr.ULoc],
    bin_to_uloc_map: dict[str, dict[int, list[mcr.ULoc]]],
    bin_size: int = mcr.BIN_SIZE,
) -> dict[genomic.Location, dict[str, list[mcr.ULoc]]]:
    """
    Join all overlapping peaks into one, i.e. keep expanding the region as we
    find overlapping peaks. Unlike finding the minimum

    Args:
        uids (list[str]): all uids across all samples
        uid_to_loc_map (dict[str, genomic.Location]): uid to genomic location
        bin_to_uloc_map (dict[int, list[str]]): location bins to uids

    Returns:
        dict[str, dict[str, str]]: _description_
    """

    # lets see what overlaps

    location_core_map: dict[
        genomic.Location, dict[str, list[mcr.ULoc]]
    ] = collections.defaultdict(lambda: collections.defaultdict(list))

    # debug for testing to end remove as it truncates list
    # locations = locations[1:(len(locations) / 4)]

    total = len(ulocs)

    total *= total

    print(f"Processing {len(ulocs)}...", file=sys.stderr)

    p = 1

    # keep track of all locations that have been allocated at least once
    allocated: set[mcr.ULoc] = set()

    # test all locations in first sample
    for uloc1 in ulocs:
        # of the form id=chrN:start-end, an lid

        # if the peaks has already been assigned to
        # an overlap, don't bother checking it as it
        # will duplicate regions
        if uloc1 in allocated:
            continue

        # sid1, _ = mcr.parse_uid(uloc1)

        # find possible overlapping locations
        # test_ulocs = mcr.get_test_ulocs(uloc1, bin_to_uloc_map, bin_size)

        # Form the largest group of overlapping peaks
        # exhausted = False

        # grouped_locations = [uloc1]

        # reset for each group search

        #
        # simplified
        #

        # test_uids = mcr.get_test_uids(uid1, loc1, bin_to_uids_map)

        # for uid2 in test_uids:
        #     # only test locations we haven't added yet
        #     if uid2 in used:
        #         continue

        #     loc2 = uid2.loc

        #     overlap = genomic.overlap_locations(loc1, loc2)

        #     if overlap is not None:
        #         # expand search to region covered by both peaks
        #         loc1 = genomic.max_region(loc1, loc2)

        #         # we found an overlap so update coords and keep searching
        #         grouped_locations.append(uid2)
        #         #used.add(uid2)

        #
        # end simplified
        #

        used = {uloc1}
        loc1 = uloc1.to_loc()

        run = True

        while run:
            test_ulocs = mcr.get_test_ulocs(uloc1, loc1, bin_to_uloc_map, bin_size)

            print("test", uloc1, test_ulocs)

            # assume we want to stop on this loop
            run = False

            for uloc2 in test_ulocs:
                # only test locations we haven't added yet
                if uloc2 in used:
                    continue

                overlap = genomic.is_overlapping(loc1, uloc2)

                #print(loc1, uloc2, overlap)

                if overlap:
                    # expand search to region covered by both peaks
                    loc1 = genomic.max_region(loc1, uloc2)

                    # we found an overlap so update coords and keep searching
                    # grouped_locations.append(uloc2)
                    used.add(uloc2)
                    # since we changed the overlap, keep running
                    run = True

                    # if there is an overlap, reset for changes
                    # break

        # now we have a list of all locations that overlap each other

        # if we have a group of entries, merge them, otherwise if the
        # location is by itself, only add it if it has not been allocated
        # to another group. This prevents duplicate entries of the whole
        # regions by itself plus any overlapping regions
        # if len(grouped_locations) > 1 or uid1 not in allocated:

        # overlap_location = str(loc1)  # f'{chr1}:{start1}-{end1}'

        for uloc in sorted(used):  # grouped_locations:
            # sid is a sample id
            # sid, _ = mcr.parse_uid(uid)  # sample_id_map[uid]

            # if the uid has not been allocated yet
            print(loc1, uloc.id, uloc)
            location_core_map[loc1][uloc.id].append(uloc)

            # used.add(uid)
            allocated.add(uloc)

    # after iterating over everything, group locations by group

    return location_core_map


def max_common_regions(
    ulocs: list[mcr.ULoc], bin_size: int = mcr.BIN_SIZE
) -> dict[genomic.Location, dict[str, list[mcr.ULoc]]]:
    print("Max region mode.", file=sys.stderr)
    return mcr.group_ulocs_by_core_location(
        ulocs, func_core_regions=_max_common_regions, bin_size=bin_size
    )


def create_max_common_region_table(
    files: list[str], bin_size: int = mcr.BIN_SIZE
) -> pd.DataFrame:
    return mcr.create_overlap_table(
        files, func_core_regions=max_common_regions, bin_size=bin_size
    )
