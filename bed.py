# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 16:39:28 2015

@author: antony
"""

import sys
import collections
import re

import gal.sample


def create_bed_header(name, color="255,0,0"):
  # Need to truncate the name to 16 chars
  
  if len(name) == 16:
    short_name = name
  else:
    short_name = name[0:13] + "..." #pychipseq.sample.get_sample_id(name)
  
  return "track type=bed name=\"" + short_name  + "\" description=\"" + name + " peaks\" color=\"" + color + "\""


def create_bedgraph_header(name, color="255,0,0"):
  # Need to truncate the name to 16 chars
  short_name = gal.sample.get_sample_id(name)
  
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
