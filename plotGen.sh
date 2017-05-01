gnuplot plots.gnu
latex plot.tex
dvipdf plot.dvi plot.pdf

mv plot.pdf imgs/RENDER.pdf
rm plot.tex plot.aux plot.log plot-inc.eps plot.dvi
