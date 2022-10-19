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


def get_overlaps(files, ids):
  """
  Takes a list of file and finds the common (if any) overlapping regions
  """

  bin_size = 10000  
  
  location_id_map = collections.defaultdict(str)
  location_core_map = collections.defaultdict(lambda: collections.defaultdict(set))
  location_chrs = collections.defaultdict(str)
  location_starts = collections.defaultdict(int)
  location_ends = collections.defaultdict(int)

  location_bins = collections.defaultdict(set)  
  
  locations = []
  
  for i in range(0, len(files)):
    file = files[i]
    id = ids[i]

    f = open(file, 'r')
    
    for line in f:
      if line[0] != '#':
        break
    
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = line.split("\t")
    
      chr = tokens[0]
      start = int(tokens[1])
      end = int(tokens[2])
    
      location = id + "=" + chr + ":" + str(start) + "-" + str(end)
    
      locations.append(location)
      
      location_id_map[location] = id
      
      location_chrs[location] = chr
      location_starts[location] = start
      location_ends[location] = end
      
      bin_start = int(start / bin_size)
      bin_end = int(end / bin_size)
      
      for bin in range(bin_start, bin_end + 1):
        location_bins[bin].add(location)
      
    f.close()
  
  # lets see what overlaps

  # debug for testing to end remove as it truncates list
  #locations = locations[1:(len(locations) / 4)]

  total = len(locations)
  
  total *= total
  
  print("Processing", str(len(locations)), str(total), "locations", file=sys.stderr);

  p = 0
  
  # keep track of all locations that have been allocated at least once
  allocated = set()
  
  for i in range(0, len(locations)):
    location1 = locations[i]
    
    chr1 = location_chrs[location1]
    start1 = location_starts[location1]
    end1 = location_ends[location1]

    # find possible overlapping locations

    

    test_locations = set()
    
    bin_start = int(start1 / bin_size)
    bin_end = int(end1 / bin_size)
    
    for bin in range(bin_start, bin_end + 1):
      for location in location_bins[bin]:
        test_locations.add(location)
    
    exhausted = False

    used = set()
    
    while not exhausted:
      # Keep restarting from a given location until all possible
      # overlaps have been found
      grouped_locations = [location1]
      start1 = location_starts[location1]
      end1 = location_ends[location1]
    
      for location2 in test_locations:
        if location1 == location2:
          continue
        
        if location2 in used:
          continue

        #sys.stderr.write(location1 + " " + location2 + "\n")      
      
        chr2 = location_chrs[location2]
        start2 = location_starts[location2]
        end2 = location_ends[location2]
        #group2 = location_group_map[location2]

        if p % 10000000 == 0:
          sys.stderr.write("p " + str(p) + "\n")

        p += 1      
      
        #if group2 != "none":
        #  continue      
        
        if chr1 != chr2:
          continue
        
        if (start1 <= start2):
          min_start = start1
          min_end = end1
          max_start = start2
          max_end = end2
        else:
          min_start = start2
          min_end = end2
          max_start = start1
          max_end = end1
          
        overlap = -1
        overlap_start = -1
        
        if max_start >= min_start and max_end <= min_end:
          overlap = max_end - max_start + 1
          overlap_start = max_start
        elif min_start < max_start and min_end > max_start:
          overlap = min_end - max_start + 1
          overlap_start = max_start
        else:
          pass
        
          
        if overlap == -1:
          continue
        
      
        # change the start1 and end1 coordinates to reflect the overlap
        # region so that each subsequent match must be within this region
        # this prevents long peaks that overlap two smaller peaks who
        # themselves do not overlap each other
        start1 = overlap_start
        end1 = overlap_start + overlap - 1
      
        grouped_locations.append(location2)
    
      # now we have a list of all locations that overlap each other
      
      # if we have a group of entries, merge them. otherwise if the
      # location is by itself, only add it if it has not been allocated
      # to another group. This prevents duplicate entries of the whole
      # region by itself plus any overlapping regions
      if len(grouped_locations) > 1 or location1 not in allocated:
        overlap_location = chr1 + ":" + str(start1) + "-" + str(end1)

        for location in grouped_locations:
          id = location_id_map[location]
      
          location_core_map[overlap_location][id].add(location)
        
          used.add(location)
          allocated.add(location)
       
      if len(grouped_locations) == 1:
        # no more to add so quit looping
        exhausted = True
      
      #sys.stderr.write("overlap " + overlap_location + " " + str(overlap) + " " + location1 + " " + location2 + " " + group + " " + str(group_id) + "\n")
  
  # after iterating over everything, group locations by group
  
  return location_core_map
