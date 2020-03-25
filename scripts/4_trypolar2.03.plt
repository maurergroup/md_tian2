#!/usr/bin/gnuplot

# intention: plot the out-of-plane scattering

E_INC = 1.92


### SPECULAR ###

set term png enhanced size 1920,1920 font "Computer Modern, 40"

set lmargin at screen 0.05
set rmargin at screen 0.85
set bmargin at screen 0.1
set tmargin at screen 0.9

set pm3d map interpolate 0,0
unset key

set output "analysis/3d-polar.png"
set multiplot

# plot the heatmap
set isosamples 500

unset border
set xtics nomirror 0.25
unset ytics
set autoscale

set angles degrees
r = 6
#set urange[0:1]   # radius
#set vrange[0:360] # angle
set xrange[-1.1:0]
set yrange[0:1.1]
set colorbox user origin 0.9,0.1 size 0.03,0.8
splot "analysis/ang_res_eloss.txt" u ((E_INC-$1)/E_INC*cos($2+90)):((E_INC-$1)/E_INC*sin($2+90)):3
set object 1 circle at first 0,0 size first 1.0 fc rgb "white" lw 5
set object 2 circle at first 0,0 size first 0.75 fc rgb "white" lw 5
set object 3 circle at first 0,0 size first 0.5 fc rgb "white" lw 5
set object 4 circle at first 0,0 size first 0.25 fc rgb "white" lw 5

set angles degrees
set label "15°" at first 1.1*cos(105), first 1.1*sin(105) front textcolor rgb "black" center
set label "30°" at first 1.1*cos(120), first 1.1*sin(120) front textcolor rgb "black" center
set label "45°" at first 1.1*cos(135), first 1.1*sin(135) front textcolor rgb "black" center
set label "60°" at first 1.1*cos(150), first 1.1*sin(150) front textcolor rgb "black" center
set label "75°" at first 1.1*cos(165), first 1.1*sin(165) front textcolor rgb "black" center

# now plot the polar grid only
set style line 11 lc rgb 'white' lw 2
unset pm3d
unset xtics 

plot -10**8 * x     lc rgb "black" lw 5 # 90°
plot -(2-sqrt(3))*x lc rgb "white" lw 3	# 75°
plot -(1/sqrt(3))*x lc rgb "white" lw 3	# 60°
plot -x             lc rgb "white" lw 3	# 45°
plot -sqrt(3)*x     lc rgb "white" lw 3 # 30°
plot -(2+sqrt(3))*x lc rgb "white" lw 3	# 15°
plot 0              lc rgb "black" lw 5 # 0°
unset multiplot


### OVER ALL AZIMUTH ANGLES ###

set output "analysis/3d-polar_azi_integrated.png"
set multiplot

set pm3d map interpolate 0,0
unset key
unset label

# plot the heatmap
set isosamples 500

unset border
set xtics nomirror 0.25
unset ytics
set autoscale

set xrange[-1.1:0]
set yrange[0:1.1]
set colorbox user origin 0.9,0.1 size 0.03,0.8
splot "analysis/polar_scatt_azi_int.txt" u ((E_INC-$1)/E_INC*cos($2+90)):((E_INC-$1)/E_INC*sin($2+90)):3

set label "15°" at first 1.1*cos(105), first 1.1*sin(105) front textcolor rgb "black" center
set label "30°" at first 1.1*cos(120), first 1.1*sin(120) front textcolor rgb "black" center
set label "45°" at first 1.1*cos(135), first 1.1*sin(135) front textcolor rgb "black" center
set label "60°" at first 1.1*cos(150), first 1.1*sin(150) front textcolor rgb "black" center
set label "75°" at first 1.1*cos(165), first 1.1*sin(165) front textcolor rgb "black" center

unset pm3d
unset xtics
plot -10**8 * x     lc rgb "black" lw 5 # 90°
plot -(2-sqrt(3))*x lc rgb "white" lw 3 # 75°
plot -(1/sqrt(3))*x lc rgb "white" lw 3 # 60°
plot -x             lc rgb "white" lw 3 # 45°
plot -sqrt(3)*x     lc rgb "white" lw 3 # 30°
plot -(2+sqrt(3))*x lc rgb "white" lw 3 # 15°
plot 0              lc rgb "black" lw 5 # 0°

unset multiplot
