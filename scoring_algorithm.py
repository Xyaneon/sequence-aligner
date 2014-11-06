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
# Last modified: 11/6/2014

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
        sm.add_up_backlink(i, 0)
    for i in range(1, sm.get_columns()):
        sm.set_score(0, i, sm.get_score(0, i - 1) + gap_score)
        sm.add_left_backlink(0, i)

def fill_matrix(sm):
    '''Uses dynamic programming to fill out a provided ScoringMatrix after
    the edges have been initialized.

    Based on pseudocode provided on page 54 of our textbook.'''
    initialize_edges(sm)
    for i in range(1, sm.get_rows()):
        for j in range(1, sm.get_columns()):
            # Calculate scores
            score_diagonal = 0
            if sm.match(i, j):
                score_diagonal = sm.get_score(i - 1, j - 1) + match_score
            else:
                score_diagonal = sm.get_score(i - 1, j - 1) + mismatch_score
            score_left = sm.get_score(i, j - 1) + gap_score
            score_up = sm.get_score(i - 1, j) + gap_score
            max_score = max(score_diagonal, score_left, score_up)
            sm.set_score(i, j, max_score)
            # Establish backlink(s)
            if max_score == score_diagonal:
                sm.add_diagonal_backlink(i, j)
            if max_score == score_left:
                sm.add_left_backlink(i, j)
            if max_score == score_up:
                sm.add_up_backlink(i, j)

def get_alignments(sm):
    '''Returns a list of the alignments generated from the scoring matrix.

    Port of code from global-grid2.rb.'''
    fill_matrix(sm)
    # Find the optimal alignment paths starting from the lower right corner.
    seq = (sm.get_top_sequence(), sm.get_left_sequence())
    print seq
    last_col = len(seq[0]) - 1
    last_row = len(seq[1]) - 1
    todo_list = [[last_row,last_col,"",""]] # Entry (row,col,string0,string1)
    done_list = [] # Entry (str0,str1)
    backlink_used = [[False for x in range(len(seq[0]) + 1)]
                      for x in range(len(seq[1]) + 1)]
    print backlink_used
    while todo_list:
        row, col, str0, str1 = todo_list.pop()
        print (row, col, str0, str1)
        backlinks = sm.get_backlinks(row, col)
        print backlinks
        if True in backlinks.values(): # If some back-link exists.
            print "({}, {}) backlinks used".format(row, col)
            backlink_used[row][col] = True # Mark linked cells as used as we go.
            if backlinks["diagonal"]:
                todo_list.append([row - 1, col - 1, seq[0][col] + str0,
                                 seq[1][row] + str1])
            if backlinks["up"]:
                todo_list.append([row - 1, col, '_' + str0,
                                 seq[1][row] + str1])
            if backlinks["left"]:
                todo_list.append([row, col - 1, seq[0][col] + str0,
                                 '_' + str1])
        else:
            print "No backlinks found in ({}, {})\n".format(row, col)
            done_list.append([str0,str1])
    # Clean up unused backlinks
    for row in range(len(seq[1])):
        print backlink_used[row]
        for col in range(len(seq[0])):
            if not backlink_used[row][col]:
                sm.remove_backlinks(row, col)
    return done_list

if __name__ == "__main__":
    # Unit test
    import terminal_output
    sm = ScoringMatrix("CGCA", "CACGTAT")
    fill_matrix(sm)
    terminal_output.print_matrix(sm)
    terminal_output.print_alignments(get_alignments(sm))
