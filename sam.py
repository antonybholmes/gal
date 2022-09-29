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

import sys
import collections
import tools
import numpy


SAMTOOLS = "/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-0.1.19/samtools"


READ_COUNT = [SAMTOOLS, "view", "-F", "0x04", "-c"]

READS_IN_REGION = [SAMTOOLS, "view", "-F", "0x04"]

PSEUDO_READ_COUNT = 10000000


def get_unique_read_count(bam_file):
    cmd = list(READ_COUNT)
    cmd.append(bam_file)

    sys.stderr.write(" ".join(cmd) + "\n")

    out = tools.cmd(cmd)

    return int(out)


def get_starts(bam_file, location):
    cmd = list(READS_IN_REGION)
    cmd.append(bam_file)
    cmd.append(location.to_string())

    #sys.stderr.write(";".join(cmd) + "\n")

    lines = tools.cmd_lines(cmd)

    ret = []

    for line in lines:
        tokens = line.strip().split("\t")

        start = int(tokens[3])

        ret.append(start)

    return ret


def bin_starts(location, starts, window, read_length):
    start_bin = location.start // window
    end_bin = location.end // window
    l = end_bin - start_bin + 1

    count_map = collections.Counter()

    for start in starts:
        sbin = start // window - start_bin
        ebin = (start + read_length) // window - start_bin

        for b in range(sbin, ebin + 1):
            if b >= 0 and b < l:
                count_map[b] += 1

    ret = numpy.zeros(l)

    for b in count_map:
        ret[b] = count_map[b]

    return ret


def get_bins(padding5p, padding3p, window):
    c = -padding5p

    ret = [c]

    c += window

    while c <= padding3p:
        ret.append(c)
        c += window

    return ret


def print_header(padding5p, padding3p, window):
    bins = get_bins(padding5p, padding3p, window)

    sys.stdout.write("ID")

    for i in range(len(bins)):
        sys.stdout.write("\t" + str(bins[i]))

    sys.stdout.write("\n")


class Sam(object):
    def __init__(self, bam_file):
        self.bam_file = bam_file
        self.read_count = get_unique_read_count(bam_file)

    def get_read_count(self):
        return self.read_count

    def get_starts(self, location):
        return get_starts(self.bam_file, location)

    def get_start_bins(self, location, window=200, read_length=101):
        return bin_starts(location, self.get_starts(location), window, read_length)

    def get_normalized_start_bins(self, location, window=200, read_length=101):
        """
        Normalize the start bins to a pseudo read count so that samples
        are comparable

        Args:
          location        a lib.genomic.Location object.
          window          integer argument for the bin size.
          read_length     integer length of each read.
        """

        bins = self.get_start_bins(location, window, read_length)

        for i in range(len(bins)):
            bins[i] = bins[i] / self.read_count * PSEUDO_READ_COUNT

        return bins
