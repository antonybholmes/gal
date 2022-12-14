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
import pandas as pd
import re

from . import genomic
from . import peaks

BIN_SIZE = 1000


def merge_peaks(file:str, out:str):
  """_summary_

  Args:
      file (str): merge overlaps in overlap file
  """

  df = pd.read_csv(file, sep='\t', header=0, keep_default_na=False)

  location_map = collections.defaultdict(list)

  c = 0

  for i in range(0, df.shape[0]):
    loc = genomic.parse_location(df.iloc[i, 0])

    # keep track of original index
    location_map[loc.chr].append((loc, i))

    c += 1

  for chr in sorted(location_map):
    location_map[chr] = sorted(location_map[chr], key=lambda loc: loc[0].start)


  data = []

  for chr in sorted(location_map):
    
    i1 = 0
    while i1 < len(location_map[chr]):
      i2 = i1 + 1

      #print(chr, i1, i2, location_map[chr])

      index1 = location_map[chr][i1][1]
      

      loc1 = location_map[chr][i1][0]

      loc2 = location_map[chr][i2][0] if i2 < len(location_map[chr]) else None

      

      if loc2 is not None and loc1.end >= loc2.start:
        index2 = location_map[chr][i2][1]

        start = loc1.start
        end = max(loc1.end, loc2.end)
        new_loc = genomic.Location(chr, start, end)

        # merge
        #sample11 = df.iloc[index1, -2]
        #sample21 = df.iloc[index1, -1]
        #sample12 = df.iloc[index2, -2]
        #sample22 = df.iloc[index2, -1]

        #print(sample21)
        #s11id, s11loc = sample11.split('=')
        #s21id, s21loc = sample21.split('=')

        #s12id, s12loc = sample12.split('=')
        #s22id, s22loc = sample22.split('=')

        # update to new location
        d1 = df.iloc[index1, :].values.copy()
        d2 = df.iloc[index2, :].values.copy()
        
        d1[0] = str(new_loc) + '//merged'
        # update width
        d1[1] = new_loc.end - new_loc.start + 1

        l1 = None
        l2 = None
        if d1[-1] != 'n/a':
          l1 = genomic.parse_location(re.search(r'(chr.+)', d1[-1]).group(1))
        if d2[-1] != 'n/a':
          l2 = genomic.parse_location(re.search(r'(chr.+)', d2[-1]).group(1))

        # two peaks to merge
        if l1 is not None and l2 is not None:
          l3 = genomic.Location(l1.chr, min(l1.start, l2.start), max(l1.end, l2.end))
          d1[-1] = re.sub(r'chr.+', str(l3), d1[-1])
        else:
          if d2[-1] != 'n/a':
            d1[-1] = d2[-1]

        # merge second column

        l1 = None
        l2 = None
        if d1[-2] != 'n/a':
          l1 = genomic.parse_location(re.search(r'(chr.+)', d1[-2]).group(1))
        if d2[-2] != 'n/a':
          l2 = genomic.parse_location(re.search(r'(chr.+)', d2[-2]).group(1))

        # two peaks to merge
        if l1 is not None and l2 is not None:
          l3 = genomic.Location(l1.chr, min(l1.start, l2.start), max(l1.end, l2.end))
          d1[-2] = re.sub(r'chr.+', str(l3), d1[-2])
        else:
          if d2[-2] != 'n/a':
            d1[-2] = d2[-2]
        

        data.append(d1)

        i1 += 1
      else:
        # copy as is
        data.append(df.iloc[index1, :].values)

      i1 += 1

  df_out = pd.DataFrame(data, columns=df.columns)
  df_out.to_csv(out, sep='\t', header=True, index=False)

