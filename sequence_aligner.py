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

import argparse
from copy import deepcopy
import os.path
import webbrowser

from scoring_matrix import ScoringMatrix
from scoring_algorithm import get_alignments
from terminal_output import print_matrix, print_alignments
from html_output import write_html

version = "v0.0.0"
desc = "sequence-aligner " + version
desc += "\nFinds semi-global alignments between FASTA sequences."
infile_help="""
Reads in the sequence from the given file path if the file exists.
Otherwise, treats this as a sequence string to align.
"""

#============================================================================
# Main program code
#============================================================================

# Set up commandline argument parser
parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=desc
            )
parser.add_argument("sequence1", help=infile_help)
parser.add_argument("sequence2", help=infile_help)
args = parser.parse_args()

# Read in sequences from FASTA files, if they exist
# Ignore first line since that's just header info, not part of the sequence
path1, path2 = args.sequence1, args.sequence2
sequence1 = sequence2 = ""
if os.path.exists(path1):
    try:
        with open(path1, 'r') as infile1_reading:
            lines1 = infile1_reading.readlines()[1:]
        for line in lines1:
            sequence1 += line.upper().strip()
    except IOError:
        print "Error: could not open first file:", path1
        exit(1)
else:
    sequence1 = path1
if os.path.exists(path2):
    try:
        with open(path2, 'r') as infile2_reading:
            lines2 = infile2_reading.readlines()[1:]
        for line in lines2:
            sequence2 += line.upper().strip()
    except IOError:
        print "Error: could not open second file:", path2
        exit(1)
else:
    sequence2 = path2

# Ensure sequence1 is always the longer one.
if len(sequence1) > len(sequence2):
    sequence1, sequence2 = sequence2, sequence1

sm = ScoringMatrix(sequence1, sequence2)
alignments = get_alignments(sm)
print_matrix(sm)
print_alignments(deepcopy(alignments))
html_file = write_html(sm, alignments)
print "Output written to", html_file
if raw_input("Open HTML output in your web browser (y/n)? ").lower() == 'y':
    webbrowser.get().open(html_file)
exit(0)
