# -*- coding: utf-8 -*-
"""
Encode read counts per base in 2 bytes

@author: Antony Holmes
"""

import sys
import collections
import numpy
import math


def encode_sam_16bit(chr_size_file, file, chromosome, read_length, window):
  chr_sizes = collections.defaultdict(int)
  
  f = open(chr_size_file, 'r')

  f.readline()
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    size = int(tokens[1])
    
    chr_sizes[chr] = size
  
  f.close()
  
  # size of this chromosome in windows
  size = (chr_sizes[chromosome] // window) + 1
  
  sys.stderr.write("Parsing " + file + " " + str(read_length) + " " + str(size) + " " + str(window) + "...\n")

  counts = numpy.zeros(size, int)
  
  f = open(file, 'r')
  
  lc = 0
  
  for line in f:
    line = line.strip()

    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[2]
    
    if chr != chromosome:
      continue

    start = int(tokens[3]) - 1
    end = start + read_length - 1
    
    win_start = start // window
    win_end = end // window
    

    for i in range(win_start, win_end + 1):
      counts[i] += 1
      
    lc += 1
    #break
    
  f.close()
  
  sys.stderr.write("Encoding...\n")
  
  # we need 2 bytes per position
  bytes = bytearray(2 * size)

  for i in range(0, size):
    p = 2 * i
    
    encode = numpy.minimum(counts[i], 65535)
    
    # upper 8 bit mask (65280) which we shift into 1 byte
    bytes[p] = (encode & 0b1111111100000000) >> 8
    
    # lower 8 bit mask (255)
    bytes[p + 1] = encode & 0b11111111

  sys.stderr.write("Writing...\n")
  
  sys.stdout.write(bytes)
  

def get_chr_sizes(chr_size_file):
  chr_sizes = collections.defaultdict(int)
  
  f = open(chr_size_file, 'r')

  f.readline()
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    size = int(tokens[1])
    
    chr_sizes[chr] = size
  
  f.close()
  
  return chr_sizes
  

def max_bits(file, chr, read_length, window):
  max_count = 0
  
  bins = collections.defaultdict(int)
  
  f = open(file, 'r')
  
  lc = 0
  
  for line in f:
    line = line.strip()

    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    c = tokens[2]
    
    if c != chr:
      continue

    start = int(tokens[3]) - 1
    end = start + read_length - 1
    
    win_start = start // window
    win_end = end // window
    

    for i in range(win_start, win_end + 1):
      bins[i] += 1
      
      if bins[i] > max_count:
        max_index = i
        max_count = bins[i]
        
  f.close()
  
  bits = math.floor(math.log(max_count, 2))
  
  if bits < 4:
    return 4
  elif bits < 8:
    return 8
  elif bits < 12:
    return 12
  elif bits < 16:
    return 16
  elif bits < 20:
    return 20
  elif bits < 24:
    return 24
  else:
    return 32
    

def encode_counts(counts, chr, read_length, bit_size):
  
  size = len(counts)
  
  sys.stderr.write("Encoding...\n")
  
  byte_size = bit_size / 8.0
  
  sys.stderr.write(str(byte_size) + " " + str(size) + "\n")
  
  # add one extra byte for overhang
  bytes = bytearray(int(byte_size * size) + 1)

  sys.stderr.write("pwer " + str(2) + " " + str(bit_size) + "\n")
  
  max_val = int(math.pow(2, bit_size) - 1)
  
  if bit_size == 4:
    p = 0
    
    even = True
    
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      if even:
        bytes[p] = encode << 4
      else:
        bytes[p] |= encode
        p += 1
        
      even = not even
      
  elif bit_size == 8:
    # One number per byte, easy case
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      bytes[i] = encode
  elif bit_size == 12:
    p = 0
    
    even = True
    
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      
      if even:
        # encode upper 8 bits in this byte and lower 4 bits in next byte
        bytes[p] = (encode & 0b111111110000) >> 4
        bytes[p + 1] = (encode & 0b1111) << 4
      else:
        # encode upper 4 bits of number in last 4 bits of this byte and
        # then 8 bits in the next byte
        
        # Need to OR with existing upper 4 bits from previous position
        bytes[p] |= ((encode & 0b111100000000) >> 8)
        bytes[p + 1] = encode & 0b11111111
        
        # Every other byte increase position count since two numbers
        # requires 3 bytes
        p += 1
        
      p += 1
      
      even = not even
      
  elif bit_size == 16:
    p = 0
    
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      # upper 8 bit mask (65280) which we shift into 1 byte
      bytes[p] = (encode & 0b1111111100000000) >> 8
      
      # lower 8 bit mask (255)
      bytes[p + 1] = encode & 0b11111111
      
      p += 2
  elif bit_size == 20:
    # We need 2.5 bytes per number
    p = 0
    
    even = True
    
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      if even:
        # encode number in all of this byte plus 4 bits of second
        bytes[p] = (encode & 0b11111111000000000000) >> 12
        bytes[p + 1] = (encode & 0b111111110000) >> 4
        bytes[p + 2] = (encode & 0b1111) << 4
      else:
        # encode 4 bits of number in last 4 bits of this byte and
        # then 8 bits in the next byte
        
        bytes[p] |= (encode & 0b11110000000000000000) >> 16
        bytes[p + 1] = (encode & 0b1111111100000000) >> 8
        bytes[p + 2] = encode & 0b11111111
        
        p += 1
       
      p += 2
      
      even = not even
  elif bit_size == 24:
    # 3 bytes per number
    p = 0
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      bytes[p] = (encode & 0b111111110000000000000000) >> 16
      bytes[p + 1] = (encode & 0b1111111100000000) >> 8
      bytes[p + 2] = encode & 0b11111111
      
      p += 3
  elif bit_size == 28:
    # We need 3.5 bytes per number
    p = 0
    
    even = True
    
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      if even:
        bytes[p] = (encode & 0b1111111100000000000000000000) >> 20
        bytes[p + 1] = (encode & 0b11111111000000000000) >> 12
        bytes[p + 2] = (encode & 0b111111110000) >> 4
        bytes[p + 3] = (encode & 0b1111) << 4
      else:
        bytes[p] |= (encode & 0b111100000000000000000000) >> 24
        bytes[p + 1] = (encode & 0b111111110000000000000000) >> 16
        bytes[p + 2] = (encode & 0b1111111100000000) >> 8
        bytes[p + 3] = encode & 0b11111111
        
        p += 1
       
      p += 3
      
      even = not even
  else:
    # 4 bytes per number 32 bits
    p = 0
    
    for i in range(0, size):
      encode = numpy.minimum(counts[i], max_val)
      
      bytes[p] = (encode & 0b11111111000000000000000000000000) >> 24
      bytes[p + 1] = (encode & 0b111111110000000000000000) >> 16
      bytes[p + 2] = (encode & 0b1111111100000000) >> 8
      bytes[p + 3] = encode & 0b11111111
      
      p += 4

  return bytes
  
  
def encode_sam(chr_size_file, file, chr, read_length, window):
  chr_sizes = get_chr_sizes(chr_size_file)
  
  # Determine how many bits we need to encode the file
  bit_size = max_bits(file, chr, read_length, window)
  
  
  # size of this chromosome in windows
  size = (chr_sizes[chr] // window) + 1
  
  sys.stderr.write("Parsing " + file + " " + str(read_length) + " " + str(size) + " " + str(window) + " " + str(bit_size) + "...\n")

  counts = numpy.zeros(size, int)
  
  f = open(file, 'r')
  
  lc = 0
  
  for line in f:
    line = line.strip()

    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    c = tokens[2]
    
    if c != chr:
      continue

    start = int(tokens[3]) - 1
    end = start + read_length - 1
    
    win_start = start // window
    win_end = end // window
    

    for i in range(win_start, win_end + 1):
      counts[i] += 1
      
    lc += 1
    #break
    
  f.close()
  
  bytes = encode_counts(counts, chr, read_length, bit_size)
  
  sys.stderr.write("Writing...\n")
  
  f = open('{}.counts.win.{}.{}bit'.format(chr, window, bit_size), 'wb')
  f.write(bytes)
  f.close()


def decode_12bit(bytes, start, l):
  p = 0
  
  even = start % 2 == 0

  ret = []
  
  for i in range(0, l):
    if even:
      ret.append((bytes[p] << 4) | ((bytes[p + 1] & 0b11110000) >> 4))
    else:
      ret.append(((bytes[p] & 0b1111) << 8) | (bytes[p + 1] & 0b11111111))
				
      p += 1
			
    p += 1
			
    even = not even
    
  return ret
