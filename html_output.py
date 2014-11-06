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
# Last modified: 11/5/2014

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
  <title>%s</title>
  <script>
    function arrow(dc,x,y,degrees) {
      dc.beginPath();
      dc.save();
      dc.translate(x,y);
      dc.rotate(-Math.PI * 2 * degrees / 360.0);
      dc.moveTo(0,#{CELL/6});
      dc.lineTo(0,-#{CELL/6});
      dc.lineTo(-#{CELL/10},-#{CELL/6} + #{CELL/10});
      dc.moveTo(0,-#{CELL/6});
      dc.lineTo(#{CELL/10},-#{CELL/6} + #{CELL/10});
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
  <canvas id="drawingCanvas" width="%d" height="%d"></canvas>
</body>
</html>
'''

CELL = 40

# Functions

def cell_score(row, col, score):
    arg = (score, CELL/2 + (col+1) * CELL, CELL*2/3 + CELL * (row + 1))
    return "dc.fillText('{0!s}',{1!s},{2!s});\n".format(arg[0], arg[1], arg[2])

def cell_left(row, col):
    arg = ((col+1) * CELL, CELL/3 + CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},90);".format(arg[0], arg[1])

def cell_top(row, col):
    arg = (CELL/3 + (col+1) * CELL}, CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},0);".format(arg[0], arg[1])

def cell_diag(row, col):
    arg = ((col+1) * CELL}, CELL * (row + 1))
    return "arrow(dc,{0!s},{1!s},45);".format(arg[0], arg[1])

def draw_cell_grid(f, seq):
    f.write("dc.beginPath();")
    # Vertical cell lines:
    (1..seq[0].size+2).each do |i|
      f.write("dc.moveTo(#{i * $CELL},#{$CELL});")
      f.write("dc.lineTo(#{i * $CELL},#{$ymax});")
    end
    # Horizontal cell lines:
    (1..seq[1].size+2).each do |i|
      f.write("dc.moveTo(#{$CELL},#{i * $CELL});")
      f.write("dc.lineTo(#{$CELL * (seq[0].size+2)},#{i * $CELL});")
    end
    f.write("dc.stroke();")
    # Sequences:
    seq[0].each_with_index do |e,i|
      f.write("dc.fillText('#{e}',#{$CELL/2 + (2+i) * $CELL},#{$CELL * 2/3});")
    end
    seq[1].each_with_index do |e,i|
      f.write("dc.fillText('#{e}',#{$CELL/2},#{$CELL * 2/3 + (i+2) * $CELL});")
    end

def cell_fill(f, row, col, score, backlink):
    links = backlink[row][col].split("")
    f.write(cell_score(display,row,col,score[row][col]))
    f.write(cell_left(display,row,col) if links.include?('l'))
    f.write(cell_top(display,row,col) if links.include?('t'))
    f.write(cell_diag(display,row,col) if links.include?('d'))

# Main function:
def draw_grid(seq, score, backlink):
    seq = seq.map{|s| s[1..-1]}
    title = "-".join(seq[0], seq[1])
    $xmax = cell * (seq[0].size+2)
    $ymax = cell * (seq[1].size+2)
    with open(title + Suffix[display], "w") as f:
        f.write(HEADER % title)
        draw_cell_grid(f, seq)
        for row in range(0, seq[1].length):
            for col in range(0, seq[0].length):
                cell_fill(f, row, col, score, backlink)
        f.write(FOOTER % [$xmax + CELL, $ymax + CELL])
