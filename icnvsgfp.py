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

from . import expression


class RnaSeqGeneICNvsGFPExpression(expression.RnaSeqExpression):
  def __init__(self):
    super().__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq_icn_vs_gfp/")
    

class CustomICNvsGFPExpression(expression.CustomExpression):
  """
  Annotate up or down in MO cells.
  """
  def __init__(self, type):
    super().__init__("RNAseq ICNvsGFP", RnaSeqGeneICNvsGFPExpression())
