sequence-aligner
================

Performs semi-global alignments on FASTA sequences. Homework 3 for Dr. Miller's Intro to Bioinformatics class.

## Usage ##

    python sequence-aligner [-h] [-g] sequence1 sequence2

The `-h` option will show help. The `-g` option will perform a global alignment instead of the default semi-global alignment. `sequence1` and `sequence2` can each be either a literal sequence string or a FASTA filename.

This program will output both its dynamic programming table (with unused backlinks cleaned up) and possible alignments found, both in the terminal and in an HTML5 file which might be nicer to look at. The program will even ask you if you want to view the latter in your browser after it finishes.

## License ##
GNU GPLv3
