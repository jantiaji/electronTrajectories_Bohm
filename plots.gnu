reset

	set terminal epslatex size 3.5,2.62 standalone color colortext #header "\\usepackage{graphicx}"
	set encoding iso_8859_1  
	set output 'plot.tex'

	#set autoscale
	#set lmargin 6
	#set rmargin 3
	#set ylabel 'Concurrence' offset 2
	#set format x '\scriptsize $%.1f$'
	#set format y '\scriptsize $%.2f$'
	
	#set key at 2.8,2.3
    #set key spacing 0.8
    #set key samplen 2
    set key box 3
    set key width 3
    #set key height 0.3
    #set key reverse Left
    
    set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
	set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
    #show grid
    
    #set yr [-6:6]   
	#set xr [-1:32.0] 
	#set zr [0:25]
	
	#set xr [-6:6]
	#set yr [0:0.2]
    #set xtics (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)	
    
    
    nx=200
    ny=200
    
    nx=nx+1
    ny=ny+1
    #set hidden3d
	#set dgrid3d ny,nx
	#set pm3d 
	#set cbrange [0:1]
	#set isosample 80
	#set view 90,00
	# GIF generator
	
	#do for [t=0:500] {# plot sprintf('punto1_%d.dat', i) using 1:4; pause 0.5 
	
	tt=500
	yy=150
	
	# Labels square
	#set label '\small t=2 ns' at 5.5,0.13
	#set label '\small t=1 ns' at -1,0.315
	#set label '\small t=0.2 ns' at 1,0.65
	#set label '\small t=0.02 ns' at 1,2.2
	#set label '\small t=0 ns' at 1.1,0.9
	
	# Labels gaussian
	set label '\small t=2 ns' at -3.0,0.18
	#set label '\small t=1 ns' at -1.5,0.315
	#set label '\small t=0.2 ns' at -0.6,0.45
	#set label '\small t=0.02 ns' at -0.35,1.34
	#set label '\small t=0 ns' at 1.1,0.9
	
	# Functions to plot
	set samples 10000
	h=1.054571*(10**(-34))
	m=9.109*(10**(-31))
	d=1*10**(-6)
	a=1/(2*(0.09*10**(-6))**2)
	A=(2*a/pi)**(1/4)
	v=1.8*10**(8)
	k=h*a/m
	t=0.36/(1*v)
	
	g(x,t)=2*(A**2)*exp(-(2*a*(x**2+(d/2)**2)/(1+4*(k*t)**2)))*(cosh(4*a*(d/2)*x/(1+4*(k*t)**2))+cos(8*a*(d/2)*x*(k*t)/(1+4*(k*t)**2)))/(sqrt(1+4*(k*t)**2))

	if(1==1){
	set xlabel '$y$ [$\mu$m]'
	set ylabel '$|\Psi_{A}+\Psi_{B}|^{2}$' offset 1.5
	set xr [-10:10]
	set yr[0:0.2]
    plot	'psi2gauss.dat' every :::tt::tt using 3:4 title '\scriptsize Numeric' with line lw 4 dashtype 1,\
    		 g(x*10**(-6),t) title '\scriptsize Analytic' with line lw 4
    }
    
    ## Wavefunction
    if(1==2){
    set xlabel '$y$ [$\mu$m]'
	set ylabel '$|\Psi_{A}+\Psi_{B}|^{2}$' offset 1.5
    #set xr [0:30]
	#set yr [0:30]
    	plot 'psi2.dat' every :::00::00 using 2:5 notitle '\scriptsize Trajectories' with line lw 4 dashtype 1,\
    }
    
    
    ## Test Gradient
    if(1==2){
    
    set key at 2.8,1.3
    #set key spacing 0.8
    #set key above
    #set key samplen 2
    #set key box 3
    #set key width 3
    #set key height 0.3
    #set key reverse Left
    
    set xlabel '$y$ [$\mu$m]'
	set ylabel '$f(y)$' offset 1.5
    #set xr [0:30]
	set yr [-1.5:1.5]
    plot	'testGrad.dat' using 1:2 title '\scriptsize $\sin(x)$' with line lw 4 dashtype 1, \
    		'testGrad.dat' using 1:3 title '\scriptsize $\nabla \sin(x)$' with line lw 4 dashtype 1
    }
    
    
    ## Trajectories
    if(1==2){
       
    set arrow from 0,0.59 to 0,5 nohead front lw 5
    set arrow from 0,-0.59 to 0,-5 nohead front lw 5
    set arrow from 0,0.41 to 0,-0.41 nohead front lw 5
	
    set key nobox
    set xlabel '$x$ [cm]'
	set ylabel '$y$ [$\mu$m]' offset 0
    
    set xr [-5:35]
	set yr [-5:5]
	
	set contour base
	set view map
	
	set palette rgbformula -7,-7,2
	set cbrange [0:0.5]
	
   	plot	'traj.dat' every :::00::1000 using 2:3:(0) notitle with line lw 0 dashtype 1, \
    		'trajPrev.dat' every :::00::1000 using (-$2):3:(0) notitle with line lw 0 dashtype 1 lc 1
    		 
   	
   	#set isosample 2
    #splot	'psi2.dat' using 2:3:4:4 with pm3d notitle
    
	#set hidden3d
	#set dgrid3d ny,nx
	#set contour
	#unset colorbox
	#set cbrange [0:1]

	
	#set parametric
	#set view 0,0,1
	#unset surface
	#set cntrparam levels 5
	#set dgrid3d
	
	#do for [tt=0:500] {
    	#splot 'psi2.dat' using 2:3:4:4 with pm3d notitle
    	#}
    }
    
    
	
