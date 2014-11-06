# This file is part of sequence-aligner.
# Copyright (C) 2014 Christopher Kyle Horton <chorton@ltu.edu>

# sequence-aligner is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# sequence-aligner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with sequence-aligner. If not, see <http://www.gnu.org/licenses/>.


# MCS 5603 Intro to Bioinformatics, Fall 2014
# Christopher Kyle Horton (000516274), chorton@ltu.edu
# Last modified: 11/4/2014

from scoring_matrix import ScoringMatrix

def print_bottom_border(columns):
    '''Prints the bottom border of a table row.'''
    output_row = "-+"
    for column in range(0, columns):
        output_row += "--+"
    print output_row

def print_matrix(sm):
    '''Prints the given ScoringMatrix to the terminal.'''
    # Sequence on top
    output_row = " |  |"
    top_sequence = sm.get_top_sequence()
    left_sequence = sm.get_left_sequence()
    for c in top_sequence:
        output_row += " " + c + "|"
    print output_row
    print_bottom_border(sm.get_columns())
    # All other rows
    for row in range(0, sm.get_rows()):
        # Top half of row
        output_row = " |"
        for column in range(0, sm.get_columns()):
            bl = sm.get_backlinks(row, column)
            output_row += "\\" if bl["diagonal"] else " "
            output_row += "^" if bl["up"] else " "
            output_row += "|"
        print output_row
        # Bottom half of row
        output_row = left_sequence[row - 1] if row >= 1 else " "
        output_row += "|"
        for column in range(0, sm.get_columns()):
            bl = sm.get_backlinks(row, column)
            output_row += "<" if bl["left"] else " "
            output_row += str(sm.get_score(row, column))
            output_row += "|"
        print output_row
        # Bottom border of row
        print_bottom_border(sm.get_columns())
