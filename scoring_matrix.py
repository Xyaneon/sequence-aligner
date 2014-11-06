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

class ScoringMatrixCell:
    '''A class implementing individual cells within the scoring matrix.'''
    def __init__(self, score=0, up=False, diagonal=False, left=False):
        '''Initializes this ScoringMatrixCell instance.

        If the arguments are omitted, a cell with no backlinks and a score of
        0 will be created by default.'''
        self.score = score
        self.up = up
        self.diagonal = diagonal
        self.left = left

    def get_score(self):
        '''Returns this cell's current score.'''
        return self.score

    def get_backlinks(self):
        '''Returns a dictionary containing this cell's current backlinks as
        Boolean values.'''
        return {"up": self.up, "diagonal": self.diagonal, "left": self.left}

    def set_score(self, score):
        '''Updates this cell's score to the new one supplied.'''
        self.score = score

    def add_up_backlink(self):
        '''Adds an up backlink to this cell.'''
        self.up = True

    def add_diagonal_backlink(self):
        '''Adds a diagonal backlink to this cell.'''
        self.diagonal = True

    def add_left_backlink(self):
        '''Adds a left backlink to this cell.'''
        self.left = True

    def remove_backlinks(self):
        '''Removes all backlinks from this cell.'''
        self.up = False
        self.diagonal = False
        self.left = False

class ScoringMatrix:
    '''A class implementing a scoring matrix.'''
    def __init__(self, sequence1, sequence2):
        '''Initializes a new scoring matrix with the supplied sequences.

        sequence1 is the sequence displayed along the left side of the scoring
        matrix, and sequence2 is displayed along the top.'''
        seql = (len(sequence1) + 1, len(sequence2) + 1)
        self.matrix = [[ScoringMatrixCell() for x in range(seql[1])]
                       for x in range(seql[0])]
        self.left_sequence = sequence1
        self.top_sequence = sequence2

    def get_top_sequence(self):
        '''Returns the sequence along the top edge of the matrix.'''
        return self.top_sequence

    def get_left_sequence(self):
        '''Returns the sequence along the left edge of the matrix.'''
        return self.left_sequence

    def get_rows(self):
        '''Returns the number of rows in this scoring matrix.

        This should be equal to the length of the left sequence, plus one.'''
        return len(self.left_sequence) + 1

    def get_columns(self):
        '''Returns the number of columns in this scoring matrix.

        This should be equal to the length of the top sequence, plus one.'''
        return len(self.top_sequence) + 1

    def get_score(self, row, column):
        '''Gets the current score at the specified row and column.'''
        try:
            return self.matrix[row][column].get_score()
        except IndexError:
            print "IndexError in get_score({0!s}, {1!s})".format(row, column)
            exit(1)

    def set_score(self, row, column, score):
        '''Sets the score at the specified row and column.'''
        try:
            self.matrix[row][column].set_score(score)
        except IndexError:
            print "IndexError in set_score({0!s}, {1!s})".format(row, column)
            exit(1)

    def get_backlinks(self, row, column):
        '''Returns the current backlinks at the specified row and column in
        the form of a dictionary of Boolean values for each direction.'''
        try:
            return self.matrix[row][column].get_backlinks()
        except IndexError:
            print "IndexError in get_backlinks({0!s}, {1!s})".format(row, column)
            exit(1)

    def add_up_backlink(self, row, column):
        '''Adds an up backlink at the specified row and column.'''
        try:
            self.matrix[row][column].add_up_backlink()
        except IndexError:
            print "IndexError in add_up_backlink({0!s}, {1!s})".format(row, column)
            exit(1)

    def add_diagonal_backlink(self, row, column):
        '''Adds a diagonal backlink at the specified row and column.'''
        try:
            self.matrix[row][column].add_diagonal_backlink()
        except IndexError:
            print "IndexError in add_diagonal_backlink({0!s}, {1!s})".format(row, column)
            exit(1)

    def add_left_backlink(self, row, column):
        '''Adds a left backlink at the specified row and column.'''
        try:
            self.matrix[row][column].add_left_backlink()
        except IndexError:
            print "IndexError in add_left_backlink({0!s}, {1!s})".format(row, column)
            exit(1)

    def remove_backlinks(self, row, column):
        '''Removes all backlinks at the specified row and column.'''
        try:
            self.matrix[row][column].remove_backlinks()
        except IndexError:
            print "IndexError in remove_backlinks({0!s}, {1!s})".format(row, column)
            exit(1)

    def match(self, row, col):
        '''Returns True iff the corresponding characters in the sequences are
        equal.

        Note that the row and column of the matrix must be supplied, not the
        character's indexes within their sequences directly.'''
        try:
            return self.left_sequence[row - 1] == self.top_sequence[col - 1]
        except IndexError:
            print "IndexError in match({0!s}, {1!s})".format(row, col)
            exit(1)

if __name__ == '__main__':
    '''Unit test for this module.'''

    def test_failed(test_number, fail_string):
        print "Test {0!s} failed: {1}".format(test_number, fail_string)
        exit(1)

    def test_passed(test_number):
        print "Test {0!s} passed.".format(test_number)

    # Test 1: ScoringMatrixCell initialization.
    test_num = 1
    test1 = ScoringMatrixCell()
    score1 = test1.get_score()
    if score1 != 0:
        fail_message = "Default score is {0!s} instead of 0".format(score1)
        test_failed(test_num, fail_message)
    backlinks1 = test1.get_backlinks()
    if backlinks1["up"]:
        fail_message = "Up backlink shouldn't exist"
        test_failed(test_num, fail_message)
    if backlinks1["diagonal"]:
        fail_message = "Diagonal backlink shouldn't exist"
        test_failed(test_num, fail_message)
    if backlinks1["left"]:
        fail_message = "Left backlink shouldn't exist"
        test_failed(test_num, fail_message)
    test_passed(test_num)

    # Test 2: Alternative ScoringMatrixCell initialization and deep copy check.
    test_num += 1
    test2 = ScoringMatrixCell(1, True, True, True)
    score2 = test2.get_score()
    backlinks2 = test2.get_backlinks()
    if score1 != 0:
        fail_message = "Default score is {0!s} instead of 0".format(score1)
        test_failed(test_num, fail_message)
    backlinks1 = test1.get_backlinks()
    if backlinks1["up"]:
        fail_message = "Up backlink shouldn't exist"
        test_failed(test_num, fail_message)
    if backlinks1["diagonal"]:
        fail_message = "Diagonal backlink shouldn't exist"
        test_failed(test_num, fail_message)
    if backlinks1["left"]:
        fail_message = "Left backlink shouldn't exist"
        test_failed(test_num, fail_message)
    if score2 != 1:
        fail_message = "Set score is {0!s} instead of 1".format(score1)
        test_failed(test_num, fail_message)
    if not backlinks2["up"]:
        fail_message = "Up backlink should exist"
        test_failed(test_num, fail_message)
    if not backlinks2["diagonal"]:
        fail_message = "Diagonal backlink should exist"
        test_failed(test_num, fail_message)
    if not backlinks2["left"]:
        fail_message = "Left backlink should exist"
        test_failed(test_num, fail_message)
    test_passed(test_num)

    # Test 3: ScoringMatrix initialization.
    test_num += 1
    seq1 = "CGCA"
    seq2 = "CACGTAT"
    test3 = ScoringMatrix(seq1, seq2)
    num_rows = test3.get_rows()
    num_cols = test3.get_columns()
    if num_rows != len(seq1) + 1:
        fail_message = "Incorrect number of rows ({0!s} instead of {1!s})".format(num_rows, len(seq1) + 1)
    if num_cols != len(seq2) + 1:
        fail_message = "Incorrect number of columns ({0!s} instead of {1!s})".format(num_cols, len(seq2) + 1)
    for row in range(num_rows):
        for column in range(num_cols):
            score3 = test3.get_score(row, column)
            backlinks3 = test3.get_backlinks(row, column)
            if score3 != 0:
                fail_message = "Default score for row {0!s}, column {1!s} is {2!s} instead of 0".format(row, column, score3)
                test_failed(test_num, fail_message)
            if backlinks3["up"]:
                fail_message = "Up backlink shouldn't exist for row {0!s}, column {1!s}".format(row, column)
                test_failed(test_num, fail_message)
            if backlinks3["diagonal"]:
                fail_message = "Diagonal backlink shouldn't exist for row {0!s}, column {1!s}".format(row, column)
                test_failed(test_num, fail_message)
            if backlinks3["left"]:
                fail_message = "Left backlink shouldn't exist for row {0!s}, column {1!s}".format(row, column)
                test_failed(test_num, fail_message)
    test_passed(test_num)

    # Test 4: Ensure deep copy of each cell in ScoringMatrix initialization.
    test_num += 1
    for row in range(num_rows):
        for column in range(num_cols):
            test3.set_score(row, column, row * column)
    for row in range(num_rows):
        for column in range(num_cols):
            if test3.get_score(row, column) != row * column:
                fail_message = "Deep copy test failed."
                test_failed(test_num, fail_message)
    test_passed(test_num)
    print "All {0!s} test cases for scoring_matrix.py passed.".format(test_num)
