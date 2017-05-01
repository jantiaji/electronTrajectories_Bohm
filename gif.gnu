reset

	#set terminal epslatex size 3.5,2.62 standalone color colortext #header "\\usepackage{graphicx}"
	#set output 'plot.tex'
	
	set terminal gif size 1920,1080 animate delay 2 font 'times' 30
	set output "animate.gif"

	set autoscale
	#set lmargin 6
	#set rmargin 3
	#set xlabel 'x [cm]' offset 0
	set xlabel 'y [{/Symbol m}m]' offset -3
	set ylabel '|{/Symbol Y}_{A}+{/Symbol Y}_{B}|^{2}'
	#set ylabel 'Concurrence' offset 2
	#set format x '\scriptsize $%.1f$'
	#set format y '\scriptsize $%.2f$'
	
	set key at 4.2,0.9
    #set key spacing 0.8
    #set key samplen 2
    #set key box 3
    #set key width -4
    #set key height 0.3
    #set key reverse Left
    
    set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
	set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
    #show grid
    
    set yr [0:1]   
	set xr [-6:6] 
	#set zr [0:1]
	
	#set xr [-12:12]
	#set yr [0:10]
    #set xtics (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)	
    
    
    nx=500
    ny=500
    
    nx=nx+1
    ny=ny+1
    #set hidden3d
	#set dgrid3d ny,nx
	#set pm3d 
	#unset colorbox
	#set cbrange [0:1]
	#set isosample 80
	#set view 90,00
	# GIF generator
	
	do for [t=0:500] {# plot sprintf('punto1_%d.dat', i) using 1:4; pause 0.5 
	
    #splot	'psi2_1D_whole.dat' every :::t::t using 2:3:4 title sprintf("t = %.3f ns", t*0.004) with line lw 4 dashtype 1 lc 1    
    plot	'psi2_1D_whole.dat' every :::t::t using 3:4 title sprintf("t = %.3f ns", t*0.004) with line lw 4 dashtype 1,\
     		#'psi2_1D_half.dat' every :::t::t using 3:4 notitle '\scriptsize $|\Psi_{A}|+|\Psi_{B}|$' with line lw 4 dashtype 1
    
    #		#sprintf('data/concurrence.dat') every ::i::i using ($1*1e3):3 notitle 'E^{+}' with points pointtype 13 lw 2 pointsize 3 lc 4;
	}
	
	
	
	
			#sprintf('concurrence.dat') using ($1*1e3):2 title '\scriptsize Experiment' with points lw 3 dashtype 3, \
			#
	
	
