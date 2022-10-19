# -*- coding: utf-8 -*-
"""
Functions related to overlaps between peaks and genes.

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

from .. import text
from .. import headings
from .. import regions
from .. import genes
  
class GeneOrientatedRegionIntersections(genes.GeneOrientatedPeaks):
  """
  Annotate genes that are part of intersections only.
  """
  
  def __init__(self):
    super().__init__("Region")
    
  
  def load_annotations(self, file):
    """
    Override to cope with checking for true intersections
    """
    
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = text.find_index(header, headings.LOCATION)
    entrez_column = text.find_index(header, headings.ENTREZ_ID)
    refseq_column = text.find_index(header, headings.REFSEQ_ID)
    symbol_column = text.find_index(header, headings.GENE_SYMBOL)
    type_column = text.find_index(header, "Relative To Gene")
    #p_column = text.find_index(header, headings.P_VALUE)
    #score_column = text.find_index(header, headings.SCORE)
    tss_column = text.find_index(header, headings.TSS_DISTANCE)
    sample_column = text.find_index(header, headings.SAMPLE)
    centromere_column = text.find_index(header, headings.CENTROMERE)
    
    sample_columns = regions.get_sample_column_count(header)

    mir_column = text.find_index(header, headings.MIR_SYMBOL)
    mir_type_column = text.find_index(header, "Relative To miR")
    mss_column = text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      type = tokens[type_column]
      p = genes.find_best_p_value(header, tokens) #float(tokens[p_column])
      score = genes.find_best_score(header, tokens) #float(tokens[score_column])
      location = tokens[location_column]
      centromere = tokens[centromere_column]
      
      intersection = True

      #sys.stderr.write("sample column " + str(sample_column) + " " + str(sample_columns) + "\n")      
      
      for i in range(0, sample_columns):
        if tokens[sample_column + i] == text.NA:
          intersection = False
          break
        
      if not intersection:
        continue
      
      entrezes = tokens[entrez_column].split(";")
      symbols = tokens[symbol_column].split(";")
      refseqs = tokens[refseq_column].split(";")
      tsses = tokens[tss_column].split(";")
      
      for i in range(0, len(refseqs)):
        refseq = refseqs[i]
        entrez = entrezes[i]
        symbol = symbols[i]
        tss = tsses[i]
        
        if refseq != "n/a":
          self.refseqs.add(refseq)
          self.collapsed_entrezes[refseq] = entrez
          self.collapsed_symbols[refseq] = symbol
          self.collapsed_p[refseq].append(p)
          self.collapsed_scores[refseq].append(score)
          self.collapsed_locations[refseq].append(location)
          self.collapsed_tss[refseq].append(tss)
          self.collapsed_types[refseq].append(type)
          self.collapsed_centromeres[refseq].append(centromere)
        
      mir = tokens[mir_column]
    
      if mir != "n/a":
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
    

class ClosestGeneOrientatedRegionIntersections(genes.GeneOrientatedPeaks):
  def __init__(self):
    super().__init__("Region")
    
  
  def load_annotations(self, file):
    """
    Override to cope with checking for true intersections
    """
    
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = text.find_index(header, headings.LOCATION)
    p_column = text.find_index(header, headings.P_VALUE)
    score_column = text.find_index(header, headings.SCORE)
    centromere_column = text.find_index(header, headings.CENTROMERE)

    entrez_column = text.find_index(header, headings.ENTREZ_ID)
    refseq_column = text.find_index(header, headings.REFSEQ_ID)
    symbol_column = text.find_index(header, headings.GENE_SYMBOL)
    type_column = text.find_index(header, headings.RELATIVE)
    tss_column = text.find_index(header, headings.TSS_DISTANCE)
    
    closest_entrez_column = text.find_index(header, headings.CLOSEST_ENTREZ_ID)
    closest_refseq_column = text.find_index(header, headings.CLOSEST_REFSEQ_ID)
    closest_symbol_column = text.find_index(header, headings.CLOSEST_GENE_SYMBOL)
    closest_type_column = text.find_index(header, headings.CLOSEST_RELATIVE)
    closest_tss_column = text.find_index(header, headings.CLOSEST_TSS_DISTANCE)
    
    
    mir_column = text.find_index(header, headings.MIR_SYMBOL)
    mir_type_column = text.find_index(header, "Relative To miR")
    mss_column = text.find_index(header, "mIR Start Distance")
    
    
    sample_column = text.find_index(header, headings.SAMPLE)
    
    sample_columns = regions.get_sample_column_count(header)

    mir_column = text.find_index(header, headings.MIR_SYMBOL)
    mir_type_column = text.find_index(header, "Relative To miR")
    mss_column = text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      
      intersection = True
      
      for i in range(0, sample_columns):
        if tokens[sample_column + i] == text.NA:
          intersection = False
          break
        
      if not intersection:
        continue
        
      type = tokens[type_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      location = tokens[location_column]
      centromere = tokens[centromere_column]
      
      if tokens[entrez_column] != text.NA:
        # Preference is given to peaks within a gene
        entrezes = tokens[entrez_column].split(";")
        symbols = tokens[symbol_column].split(";")
        refseqs = tokens[refseq_column].split(";")
        tsses = tokens[tss_column].split(";")
        type = tokens[type_column]
      else:
        entrezes = tokens[closest_entrez_column].split(";")
        symbols = tokens[closest_symbol_column].split(";")
        refseqs = tokens[closest_refseq_column].split(";")
        tsses = tokens[closest_tss_column].split(";")
        type = tokens[closest_type_column]
        
      
      for i in range(0, len(refseqs)):
        refseq = refseqs[i]
        entrez = entrezes[i]
        symbol = symbols[i]
        tss = tsses[i]
        
        if refseq != text.NA:
          self.refseqs.add(refseq)
          self.collapsed_entrezes[refseq] = entrez
          self.collapsed_symbols[refseq] = symbol
          self.collapsed_p[refseq].append(p)
          self.collapsed_scores[refseq].append(score)
          self.collapsed_locations[refseq].append(location)
          self.collapsed_tss[refseq].append(tss)
          self.collapsed_types[refseq].append(type)
          self.collapsed_centromeres[refseq].append(centromere)
        
      mir = tokens[mir_column]
    
      if mir != text.NA:
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
    

class GeneOrientatedRegionCores(genes.GeneOrientatedPeaks):
  """
  For annotating peaks that lie within any of the three way intersections
  """
  
  def __init__(self):
    super().__init__("Region")
    
  
  def load_annotations(self, file):
    """
    Override to cope with checking for true intersections
    """
    
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = text.find_index(header, headings.LOCATION)
    entrez_column = text.find_index(header, headings.ENTREZ_ID)
    refseq_column = text.find_index(header, headings.REFSEQ_ID)
    symbol_column = text.find_index(header, headings.GENE_SYMBOL)
    type_column = text.find_index(header, "Relative To Gene")
    p_column = text.find_index(header, headings.P_VALUE)
    score_column = text.find_index(header, headings.SCORE)
    tss_column = text.find_index(header, headings.TSS_DISTANCE)
    centromere_column = text.find_index(header, headings.CENTROMERE)
    sample_column = text.find_index(header, headings.SAMPLE)
    
    mir_column = text.find_index(header, headings.MIR_SYMBOL)
    mir_type_column = text.find_index(header, "Relative To miR")
    mss_column = text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      type = tokens[type_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      location = tokens[location_column]
      centromere = tokens[centromere_column]
      
      intersection = False
      
      # One of the intersections must be true
      if (tokens[sample_column] != "n/a" and tokens[sample_column + 1] != "n/a") or \
        (tokens[sample_column] != "n/a" and tokens[sample_column + 2] != "n/a") or \
        (tokens[sample_column + 1] != "n/a" and tokens[sample_column + 2] != "n/a"):
        intersection = True
        
      if not intersection:
        continue
      
      
      entrezes = tokens[entrez_column].split(";")
      symbols = tokens[symbol_column].split(";")
      refseqs = tokens[refseq_column].split(";")
      tsses = tokens[tss_column].split(";")
      
      for i in range(0, len(refseqs)):
        entrez = entrezes[i]
        symbol = symbols[i]
        refseq = refseqs[i]
        tss = tsses[i]
        
        if refseq != "n/a":
          self.refseqs.add(refseq)
          self.collapsed_entrezes[refseq] = entrez
          self.collapsed_symbols[refseq] = symbol
          self.collapsed_p[refseq].append(p)
          self.collapsed_scores[refseq].append(score)
          self.collapsed_locations[refseq].append(location)
          self.collapsed_tss[refseq].append(tss)
          self.collapsed_types[refseq].append(type)
          self.collapsed_centromeres[refseq].append(centromere)
        
      mir = tokens[mir_column]
    
      if mir != "n/a":
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
