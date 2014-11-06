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

HEADER = '''
<!DOCTYPE html>
<html lang="en">
<head>
  <!-- saved from url=(0014)about:internet -->
  <meta charset="utf-8">
  <title>{0}</title>
  <script>
    function arrow(dc,x,y,degrees) {
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
    }
    window.onload = function() {
    var canvas = document.getElementById("drawingCanvas");
    var dc = canvas.getContext("2d");
    dc.font = '10pt Helvetica';
    dc.textAlign = 'center';
'''
FOOTER = '''
};
  </script>
</head>
<body>
  <canvas id="drawingCanvas" width="{0!s}" height="{1!s}"></canvas>
</body>
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
    return "arrow(dc,{0!s},{1!s},90);".format(arg[0], arg[1])

def cell_top(row, col):
    arg = (CELL/3 + (col+1) * CELL, CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},0);".format(arg[0], arg[1])

def cell_diag(row, col):
    arg = ((col+1) * CELL, CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},45);".format(arg[0], arg[1])

def draw_cell_grid(f, seq):
    f.write("dc.beginPath();")
    # Vertical cell lines:
    for i in range(1, seq[0].size + 2):
        f.write("dc.moveTo({0!s},{1!s});".format(i * CELL, CELL))
        f.write("dc.lineTo({0!s},{1!s});".format(i * CELL, ymax))
    # Horizontal cell lines:
    for i in range(1, seq[1].size + 2):
        f.write("dc.moveTo({0!s},{1!s});".format(CELL, i * CELL))
        f.write("dc.lineTo({0!s},{1!s});".format(CELL * (seq[0].size+2), i * CELL))
    f.write("dc.stroke();")
    # Sequences:
    for i in range(0, seq[0].length):
        arg = (seq[0][i], CELL/2 + (2+i) * CELL, CELL * 2/3)
        f.write("dc.fillText('{0}',{1!s},{2!s});".format(arg[0], arg[1], arg[2]))
    for i in range(0, seq[1].length):
        arg = (seq[1][i], CELL/2, CELL * 2/3 + (i+2) * CELL)
        f.write("dc.fillText('{0}',{1!s},{2!s});".format(arg[0], arg[1], arg[2]))

def cell_fill(f, row, col, sm):
    links = sm.get_backlinks()
    f.write(cell_score(row, col, sm.get_score(row, col)))
    if links["left"]:
        f.write(cell_left(row, col))
    if links["top"]:
        f.write(cell_top(row, col))
    if links["diagonal"]:
        f.write(cell_diag(row, col))

# Main function:
def draw_grid(sm):
    #seq = seq.map{|s| s[1..-1]}
    seq = [sm.get_top_sequence(), sm.get_left_sequence()]
    title = "-".join(seq[0], seq[1])
    xmax = cell * (seq[0].size+2)
    ymax = cell * (seq[1].size+2)
    with open(title + ".html", "w") as f:
        f.write(HEADER.format(title, CELL/6, CELL/10))
        draw_cell_grid(f, seq)
        for row in range(0, seq[1].length):
            for col in range(0, seq[0].length):
                cell_fill(f, row, col, sm)
        f.write(FOOTER.format(xmax + CELL, ymax + CELL))

if __name__ == "__main__":
    # Unit testing
    from scoring_algorithm import fill_matrix
    sm = ScoringMatrix("CGCA", "CACGTAT")
    fill_matrix(sm)
    draw_grid(sm)
