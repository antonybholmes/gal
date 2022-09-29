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

import re
from . import text

def get_sample_column_count(header:list[str]) -> int:
  """
  Returns the number of sample columns in a gene file.

  Args:
      header (list[str]): table column header

  Returns:
      int: index of first column containing "Sample"
  """
  ret = 0
  
  c = text.find_index(header, "Sample")
  
  for i in range(c, len(header)):
    if not re.match(r'^Sample.*', header[i]):
      break
    
    ret += 1
    
  return ret
