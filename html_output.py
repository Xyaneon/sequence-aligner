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

# This module is based on code from draw-grib.rb, located at:
# http://medicalopensource.net/mcs5603/recap.html
# This is basically a port of Ruby code written by Professor Miller, with some
# modifications. For example, this implementation only supports HTML5 output.

HEADER = '''<!DOCTYPE html>
<html lang="en">
<head>
    <!-- saved from url=(0014)about:internet -->
    <meta charset="utf-8">
    <title>{0}</title>
    <script>
        function arrow(dc,x,y,degrees) {{
            dc.beginPath();
            dc.save();
            dc.translate(x,y);
            dc.rotate(-Math.PI * 2 * degrees / 360.0);
            dc.moveTo(0,{1!s});
            dc.lineTo(0,-{1!s});
            dc.lineTo(-{2!s},-{1!s} + {2!s});
            dc.moveTo(0,-{1!s});
            dc.lineTo({2!s},-{1!s} + {2!s});
            dc.stroke();
            dc.restore();
        }}
        window.onload = function() {{
        var canvas = document.getElementById("drawingCanvas");
        var dc = canvas.getContext("2d");
        dc.font = '10pt Helvetica';
        dc.textAlign = 'center';
'''
BODYTOP = '''        }};
    </script>
</head>
<body>
    <style type="text/css" media="print">
    .printbutton {{
      visibility: hidden;
      display: none;
    }}
    </style>
    <h1>{} Alignment</h1>
    <script>
        document.write("<input type='button' " +
        "onClick='window.print()' " +
        "class='printbutton' " +
        "value='Print'/>");
    </script>
    <h2>Dynamic programming table</h2>
'''
CANVAS = '  <canvas id="drawingCanvas" width="{0!s}" height="{1!s}"></canvas>'
BODYBOTTOM = '''</body>
</html>
'''

CELL = 40

xmax = ymax = 0

# Functions

def cell_score(row, col, score):
    arg = (score, CELL/2 + (col+1) * CELL, CELL*2/3 + CELL * (row + 1))
    return "dc.fillText('{0!s}',{1!s},{2!s});\n".format(arg[0], arg[1], arg[2])

def cell_left(row, col):
    arg = ((col+1) * CELL, CELL/3 + CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},90);\n".format(arg[0], arg[1])

def cell_top(row, col):
    arg = (CELL/3 + (col+1) * CELL, CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},0);\n".format(arg[0], arg[1])

def cell_diag(row, col):
    arg = ((col+1) * CELL, CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},45);\n".format(arg[0], arg[1])

def draw_cell_grid(f, seq):
    f.write("dc.beginPath();\n")
    # Vertical cell lines:
    for i in range(1, len(seq[0]) + 3):
        f.write("dc.moveTo({0!s},{1!s});\n".format(i * CELL, CELL))
        f.write("dc.lineTo({0!s},{1!s});\n".format(i * CELL, CELL * (len(seq[1]) + 2)))
    # Horizontal cell lines:
    for i in range(1, len(seq[1]) + 3):
        f.write("dc.moveTo({0!s},{1!s});\n".format(CELL, i * CELL))
        f.write("dc.lineTo({0!s},{1!s});\n".format(CELL * (len(seq[0]) + 2), i * CELL))
    f.write("dc.stroke();\n")
    # Sequences:
    for i in range(0, len(seq[0])):
        arg = (seq[0][i], CELL/2 + (2+i) * CELL, CELL * 2/3)
        f.write("dc.fillText('{0}',{1!s},{2!s});\n".format(arg[0], arg[1], arg[2]))
    for i in range(0, len(seq[1])):
        arg = (seq[1][i], CELL/2, CELL * 2/3 + (i+2) * CELL)
        f.write("dc.fillText('{0}',{1!s},{2!s});\n".format(arg[0], arg[1], arg[2]))

def cell_fill(f, row, col, sm):
    links = sm.get_backlinks(row, col)
    f.write(cell_score(row, col, sm.get_score(row, col)))
    if links["left"]:
        f.write(cell_left(row, col))
    if links["up"]:
        f.write(cell_top(row, col))
    if links["diagonal"]:
        f.write(cell_diag(row, col))

def draw_grid(f, sm):
    seq = [sm.get_top_sequence(), sm.get_left_sequence()]
    title = "-".join((seq[0], seq[1]))
    xmax = CELL * (len(seq[0]) + 2)
    ymax = CELL * (len(seq[1]) + 2)
    draw_cell_grid(f, seq)
    for row in range(0, len(seq[1]) + 1):
        for col in range(0, len(seq[0]) + 1):
            cell_fill(f, row, col, sm)
    return (xmax, ymax)

def write_alignments(f, alignments):
    '''Writes the alignmentsfound to the HTML output.'''
    i = 1
    while alignments:
        seq = alignments.pop()
        f.write("  <h2>Alignment #{}</h2>\n".format(i))
        for s in seq:
            f.write("  <code>{}</code><br />\n".format(s))
        i += 1

# Main function:
def write_html(sm, alignments, alignment_is_global=False):
    '''Puts together the HTML file for the table and alignments.'''
    seq = [sm.get_top_sequence(), sm.get_left_sequence()]
    align_type = "Global" if alignment_is_global else "Semi-Global"
    #title = "-".join((seq[0], seq[1], align_type.lower()))
    title= "output"
    filename = title + ".html"
    with open(filename, "w") as f:
        f.write(HEADER.format(title, CELL/6, CELL/10))
        xmax, ymax = draw_grid(f, sm)
        f.write(BODYTOP.format(align_type))
        f.write(CANVAS.format(xmax + CELL, ymax + CELL))
        write_alignments(f, alignments)
        f.write(BODYBOTTOM)
    return filename

if __name__ == "__main__":
    # Unit testing
    from scoring_algorithm import fill_matrix, get_alignments
    sm = ScoringMatrix("CGCA", "CACGTAT")
    alignments = get_alignments(sm)
    html_file = write_html(sm, alignments)
