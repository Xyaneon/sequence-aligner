sequence-aligner
================

Performs semi-global alignments on FASTA sequences. Homework 3 for Dr. Miller's Intro to Bioinformatics class.

## Usage ##

    python sequence-aligner [-h] [-g] [-v] sequence1 sequence2

The `-h` option will show help. The `-g` option will perform a global alignment instead of the default semi-global alignment. The `-v` option will show the HTML5 output file automatically in your web browser without asking you at the end. `sequence1` and `sequence2` can each be either a literal sequence string or a FASTA filename.

This program will output both its dynamic programming table (with unused backlinks cleaned up) and possible alignments found, both in the terminal and in an HTML5 file which might be nicer to look at. The program will even ask you if you want to view the latter in your browser after it finishes, if you haven't already told it to do so via the `-v` option.

The generated HTML5 file will also show a button you can use to print the page when you view it in a browser. Handy!

## License ##
GNU GPLv3
