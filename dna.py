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
import urllib
import json

#link = "http://156.145.14.248:8080/dna/api/v2/hg19/chr1/100000/101000/f/u/n"
BASE_LINK = "http://156.145.14.248:8080/dna/api/v2/hg19/"

def get_sequence_with_mask(chr, start, end):
  link = BASE_LINK + chr + "/" + str(start) + "/" + str(end) + "/f/u/n"  

  #sys.stderr.write(link + "\n")  
  
  f = urllib.urlopen(link)
  json_data = f.read()

  data = json.loads(json_data)
  
  sequence = data[0]["sequence"]
  
  return sequence


def score_mask_seq(seq):
  '''
  Score a sequence for percentage of Ns
  '''
  
  p = 0.0
  
  for c in seq:
    #sys.stderr.write(c + "\n")
    if c == 'N':
      p += 1.0
  
  p /= len(seq)
  
  #sys.stderr.write(str(p) + "\n")
  
  return p