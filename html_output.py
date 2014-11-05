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
      dc.moveTo(0,#{Cell[:html5]/6});
      dc.lineTo(0,-#{Cell[:html5]/6});
      dc.lineTo(-#{Cell[:html5]/10},-#{Cell[:html5]/6} + #{Cell[:html5]/10});
      dc.moveTo(0,-#{Cell[:html5]/6});
      dc.lineTo(#{Cell[:html5]/10},-#{Cell[:html5]/6} + #{Cell[:html5]/10});
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

# Functions

def cell_score(row, col, score):
    return "dc.fillText('#{score}'," + "#{$cell/2 + (col+1) * $cell},#{$cell*2/3 + $cell * (row + 1)});\n"

def cell_left(row, col):
    return "arrow(dc,#{(col+1) * $cell},#{$cell/3 + $cell * (row + 1)},90);"

def cell_top(row, col):
    return "arrow(dc,#{$cell/3 + (col+1) * $cell},#{$cell * (row + 1)},0);"

def cell_diag(row, col):
    return "arrow(dc,#{(col+1) * $cell},#{$cell * (row + 1)},45);"

def draw_cell_grid(f, seq):
    f.puts "dc.beginPath();"
    # Vertical cell lines:
    (1..seq[0].size+2).each do |i|
      f.puts "dc.moveTo(#{i * $cell},#{$cell});"
      f.puts "dc.lineTo(#{i * $cell},#{$ymax});"
    end
    # Horizontal cell lines:
    (1..seq[1].size+2).each do |i|
      f.puts "dc.moveTo(#{$cell},#{i * $cell});"
      f.puts "dc.lineTo(#{$cell * (seq[0].size+2)},#{i * $cell});"
    end
    f.puts "dc.stroke();"
    # Sequences:
    seq[0].each_with_index do |e,i|
      f.puts "dc.fillText('#{e}',#{$cell/2 + (2+i) * $cell},#{$cell * 2/3});"
    end
    seq[1].each_with_index do |e,i|
      f.puts "dc.fillText('#{e}',#{$cell/2},#{$cell * 2/3 + (i+2) * $cell});"
    end

def cell_fill(f, row, col, score, backlink):
    links = backlink[row][col].split("")
    f.puts cell_score(display,row,col,score[row][col])
    f.puts cell_left(display,row,col) if links.include?('l')
    f.puts cell_top(display,row,col) if links.include?('t')
    f.puts cell_diag(display,row,col) if links.include?('d')

# Main function:
def draw_grid(seq, score, backlink):
    cell = 40
    seq = seq.map{|s| s[1..-1]}
    title = seq[0].join + "-" + seq[1].join
    $xmax = cell * (seq[0].size+2)
    $ymax = cell * (seq[1].size+2)
    with open(title + Suffix[display], "w") as f:
        f.write(HEADER % title)
        draw_cell_grid(f, seq)
        for row in range(0, seq[1].length):
            for col in range(0, seq[0].length):
                cell_fill(f, row, col, score, backlink)
        f.write(FOOTER % [$xmax + cell, $ymax + cell])
