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

Copyright (C) 2024 Antony Holmes.
"""

import sqlite3
from typing import Union
from .genomic import Location

# the use of MIN(start, end) is because the coordinates are stranded
# so if strand is '-' then start is greater than end. For convention
# we report coordinates always on the forward strand so MIN(start, end)
# ensures the forward strand start is always used regardless of gene
# orientation

WITHIN_SQL = ("SELECT id, chr, start, end, strand, gene_id, gene_symbol, MIN(start, end) - ? "
              "FROM genes "
              "WHERE level = 1 AND chr = ? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) "
              "ORDER BY start ASC")

CLOSEST_SQL = ("SELECT id, chr, start, end, strand, gene_id, gene_symbol, stranded_start - ? "
               "FROM genes "
               "WHERE level=1 AND chr=? "
               "ORDER BY ABS(stranded_start - ?) ASC "
               "LIMIT ?")

class Loctogene:
    def __init__(self, db: str):
        self._db = db
        self._conn = sqlite3.connect(db)
        self._cursor = self._conn.cursor()

    def close(self):
        self._cursor.close()
        self._conn.close()

    def get_genes_within(self, loc: Location) -> list[dict[str, Union[str, int]]]:
        mid = int((loc.start + loc.end) / 2)

        self._conn.set_trace_callback(print)
        rows = self._cursor.execute(
            WITHIN_SQL,
            [mid, loc.chr, loc.start, loc.start, loc.end, loc.end],
        ).fetchall()
     
        return Loctogene._format_rows(rows)

    def get_closest_genes(
        self, loc: Location, n: int = 10
    ) -> list[dict[str, Union[str, int]]]:
        mid = int((loc.start + loc.end) / 2)

        rows = self._cursor.execute(
            CLOSEST_SQL,
            [mid, loc.chr, mid, n],
        ).fetchall()

        return Loctogene._format_rows(rows)

    @staticmethod
    def _format_rows(rows: list[any]) -> list[dict[str, Union[str, int]]]:
        # swap start and end around based on strand so that coordinates
        # are always returned as if on forward strand
 
        return [
            {
                "id": row[0],
                "chr": row[1],
                "start": row[2] if row[4] == "+" else row[3],
                "end": row[3] if row[4] == "+" else row[2],
                "strand": row[4],
                "gene_id": row[5],
                "gene_symbol": row[6],
                "dist": row[7],
            }
            for row in rows
        ]
