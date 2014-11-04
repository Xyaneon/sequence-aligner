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

match_score = 1
gap_score = -1
mismatch_score = 0
terminal_gap_score = -1

def initialize_edges(sm):
    '''Sets up the top and left edges of the provided ScoringMatrix. This is
    the first step in the dynamic programming algorithm.

    Based on pseudocode provided on page 54 of our textbook.'''
    sm.set_score(0, 0, 0)
    for i in range(1, sm.get_rows()):
        sm.set_score(i, 0, sm.get_score(i - 1, 0) + gap_score)
    for j in range(1, sm.get_cols()):
        sm.set_score(0, i, sm.get_score(0, i - 1) + gap_score)

def fill_matrix(sm):
    '''Uses dynamic programming to fill out a provided ScoringMatrix after
    the edges have been initialized.

    Based on pseudocode provided on page 54 of our textbook.'''
    initialize_edges(sm)
    for i in range(1, sm.get_rows()):
        for j in range(1, sm.get_cols()):
            score1 = 0
            if sm.match(i, j):
                score1 = sm.get_score(i - 1, j - 1) + match_score
            else:
                score1 = sm.get_score(i - 1, j - 1) + mismatch_score
            score2 = sm.get_score(i, j - 1) + gap_score
            score3 = sm.get_score(i - 1, j) + gap_score
            sm.set_score(i, j, max(score1, score2, score3))

if __name__ == "__main__":
    # Unit test
    import terminal_output
    sm = ScoringMatrix("CGCA", "CACGTAT")
    sm.fill_matrix()
    terminal_output.print_matrix(sm)
