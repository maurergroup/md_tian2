#!/usr/bin/gnuplot


### SPECULAR ###

set term png enhanced size 1920,1920 font "Computer Modern, 40"

set lmargin at screen 0.05
set rmargin at screen 0.85
set bmargin at screen 0.1
set tmargin at screen 0.9

set pm3d map interpolate 0,0
unset key

set output "analysis/3d-polar_wp.png"
set multiplot

# plot the heatmap
set isosamples 500

unset border
set xtics nomirror 0.25
unset ytics
set autoscale

set angles degrees
r = 6
set xlabel "-E_f/E_i"
#set urange[0:1]   # radius
#set vrange[0:360] # angle
set xrange[-1.1:0]
set yrange[0:1.1]
set key top left Left reverse samplen -1.2 width -2.9 box lw 5
set colorbox user origin 0.9,0.1 size 0.03,0.8
plot "analysis/component_fast.txt" u ($1*cos($2+90)):($1*sin($2+90)) t "Barrier not crossed" w p ps 3 lw 2,\
     "analysis/component_slow_single.txt" u ($1*cos($2+90)):($1*sin($2+90)) t "Barrier crossed / S"w p ps 3 lw 2 ,\
     "analysis/component_slow_multi.txt"  u ($1*cos($2+90)):($1*sin($2+90)) t "Barrier crossed / M" w p ps 3 lw 2
set object 1 circle at first 0,0 size first 1.0 fc rgb "black" lw 5
set object 2 circle at first 0,0 size first 0.75 fc rgb "black" lw 5
set object 3 circle at first 0,0 size first 0.5 fc rgb "black" lw 5
set object 4 circle at first 0,0 size first 0.25 fc rgb "black" lw 5

set arrow from 0                     ,1                     to 0,0 nohead lw 5
set arrow from -(sqrt(3)-1)/2/sqrt(2),(sqrt(3)+1)/2/sqrt(2) to 0,0 nohead lw 5
set arrow from -0.5                  , sqrt(3)/2            to 0,0 nohead lw 5
set arrow from -cos(45)              , sin(45)              to 0,0 nohead lw 5
set arrow from -sqrt(3)/2            ,0.5                   to 0,0 nohead lw 5
set arrow from -(sqrt(3)+1)/2/sqrt(2),(sqrt(3)-1)/2/sqrt(2) to 0,0 nohead lw 5
set arrow from -1                    ,0                     to 0,0 nohead lw 5

set angles degrees
set label "15°" at first 1.05*cos(105), first 1.05*sin(105) front textcolor rgb "black" center
set label "30°" at first 1.05*cos(120), first 1.05*sin(120) front textcolor rgb "black" center
set label "45°" at first 1.05*cos(135), first 1.05*sin(135) front textcolor rgb "black" center
set label "60°" at first 1.05*cos(150), first 1.05*sin(150) front textcolor rgb "black" center
set label "75°" at first 1.05*cos(165), first 1.05*sin(165) front textcolor rgb "black" center

# now plot the polar grid only
set style line 11 lc rgb 'white' lw 2
unset pm3d
unset xtics 
unset xlabel

#plot (sqrt( (1E2*x)**2       + x**2)<1?-1E2*x    :NaN) lc rgb "black" lw 5 # 90°
#plot (sqrt((-(2-sqrt(3))*x)**2 + x**2)<1?-(2-sqrt(3))*x:NaN) lc rgb "black" lw 3	# 75°
#plot (sqrt((-(1/sqrt(3))*x)**2 + x**2)<1?-(1/sqrt(3))*x:NaN) lc rgb "black" lw 3	# 60°
#plot (sqrt( x**2               + x**2)<1?-x            :NaN) lc rgb "black" lw 3	# 45°
#plot (sqrt((-sqrt(3)*x)**2     + x**2)<1?-sqrt(3)*x    :NaN) lc rgb "black" lw 3 # 30°
#plot (sqrt((-(2+sqrt(3))*x)**2 + x**2)<1?-(2+sqrt(3))*x:NaN) lc rgb "black" lw 3	# 15°
#plot (sqrt( -(1E-10*x)**2      + x**2)<1?-1E-10*x       :NaN) lc rgb "black" lw 5 # 0°
plot NaN notitle
unset multiplot

