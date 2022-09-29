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

import subprocess

def cmd(cmd):
    """
    Run a command line program and return the output as a string.

    Args:
      cmd     An array of strings to form an command line expression
              such as ['cat', 'file.txt']
    """
    return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read()


def cmd_lines(cmd):
    """
    Run a command line program and return the output as a list of
    string.

    Args:
      cmd     An array of strings to form an command line expression
              such as ['cat', 'file.txt']
    """
    return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.readlines()
