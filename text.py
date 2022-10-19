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

from typing import List
import numpy as np

NA = "n/a"
EMPTY_ARRAY = np.array([])


def get_header(f):
    return f.readline().strip().split("\t")


def find_index(tokens: List[str], text: str, offset: int = 0):
    """
    Find the first heading in list that matches some text.
    """

    idx = find_indices(tokens, text, offset)

    if idx.size > 0:
        return idx[0]
    else:
        return -1


def find_indices(tokens: List[str], text: str) -> int:
    """
    Find all the headings in list that matches some text.
    """

    ltokens = [x.lower() for x in tokens]
    lt = text.lower()

    idx = np.where([lt in x for x in ltokens])[0]

    if idx.size > 0:
        return idx
    else:
        return EMPTY_ARRAY


def empty_line(n: int) -> str:
    """
    Produce an empty line of n tab separated n/as
    Args:
        n:  number of columns
    Return
        Tab separated string of n n/as
    """

    return '\t'.join([NA] * n)
