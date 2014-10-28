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
# Last modified: 10/28/2014

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

class ScoringMatrix:
    '''A class implementing a scoring matrix.'''
    def __init__(self, sequence_length1, sequence_length2):
        '''Initializes a new scoring matrix with the supplied sequence
        lengths.'''
        seql = (sequence_length1 + 1, sequence_length2 + 1)
        self.matrix = [[ScoringMatrixCell() for x in range(seql[0])]
                       for x in range(seql[1])]

    def get_score(self, row, column):
        '''Gets the current score at the specified row and column.'''
        return self.matrix[row][column].get_score()

    def set_score(self, row, column, score):
        '''Sets the score at the specified row and column.'''
        self.matrix[row][column].set_score(score)

    def get_backlinks(self, row, column):
        '''Returns the current backlinks at the specified row and column in
        the form of a dictionary of Boolean values for each direction.'''
        self.matrix[row][column].get_backlinks()

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
