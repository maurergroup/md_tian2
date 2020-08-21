#!/usr/bin/gnuplot

# intention: plot all beforehand created data files to analyze the scattering

E_INC = 1.92

set term png enhanced size 1920,1080 font "Computer Modern, 30"

set xlabel "Number of bounces"
set ylabel "Probability density"
set output "analysis/bounces.png"
plot "analysis/bounces.txt" u 1:3 t "scattered"   w l lw 5,\
     "analysis/bounces.txt" u 1:4 t "absorbed"    w l lw 5,\
     "analysis/bounces.txt" u 1:5 t "transmitted" w l lw 5,\
     "analysis/bounces.txt" u 1:2 t "all"         w l lw 5 lt -1


set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "Probability density/eV^{-1}"
set output "analysis/eloss.png"
plot "analysis/eloss.txt" u 1:3 t "single bounce" w l lw 5,\
     "analysis/eloss.txt" u 1:4 t "double bounce" w l lw 5,\
     "analysis/eloss.txt" u 1:5 t "multi bounce"  w l lw 5,\
     "analysis/eloss.txt" u 1:2 t "all"           w l lw 5 lt -1


set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "Probability density/eV^{-1}"
set output "analysis/spec_eloss.png"
plot "analysis/spec_eloss.txt" u 1:3 t "single bounce" w l lw 5,\
     "analysis/spec_eloss.txt" u 1:4 t "double bounce" w l lw 5,\
     "analysis/spec_eloss.txt" u 1:5 t "multi bounce"  w l lw 5,\
     "analysis/spec_eloss.txt" u 1:2 t "all"           w l lw 5 lt -1


set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "Probability density/eV^{-1}"
set output "analysis/in_plane_eloss_1d.png"
plot "analysis/in_plane_eloss.txt" u 1:3 t "single bounce" w l lw 5,\
     "analysis/in_plane_eloss.txt" u 1:4 t "double bounce" w l lw 5,\
     "analysis/in_plane_eloss.txt" u 1:5 t "multi bounce"  w l lw 5,\
     "analysis/in_plane_eloss.txt" u 1:2 t "all"           w l lw 5 lt -1


set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "Probability density/eV^{-1}"
set output "analysis/eloss_to_ehps.png"
plot "analysis/eloss_to_ehps.txt" notitle w l lw 5


set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "Probability density/eV^{-1}"
set output "analysis/eloss_to_ehps_spec.png"
plot "analysis/eloss_to_ehps_spec.txt" notitle w l lw 5


set xlabel "{/cmti10 r}_{sp}/Å"
set ylabel "Probability density/Å^{-1}"
set output "analysis/min_proj_surf_dist.png"
plot "analysis/ps_dist.txt" notitle w l lw 5


set xlabel "Polar scattering angle"
set ylabel "counts"
set title "In plane scattering polar angle distribution"
set output "analysis/ang_res_occ.png"
plot "analysis/ang_res_occurrence.txt" w l lw 5 notitle


set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "Scattering angle"
set cblabel "counts"
set pm3d map interpolate 0,0
set yrange [] reverse writeback
set title "In plane scattering angle vs energy loss"
set output "analysis/in_plane_eloss.png"
splot "analysis/ang_res_eloss.txt" notitle



reset
set xlabel "Number of bounces"
set ylabel "E_{loss}/eV"
set pm3d map interpolate 1,1
set output "analysis/bounces_vs_eloss.png"
splot "analysis/bounces_vs_eloss.txt" u ($1+0.5):2:3 notitle

set term png enhanced size 2000,1200 font "Computer Modern, 30"
set pm3d map interpolate 0,0
set xlabel "{/Symbol D} out of plane / °"
set ylabel "{/Symbol D} in plane / °"
set output "analysis/spherical_rel.png"
splot "analysis/rel_spherical_symmetry.txt" notitle

set xlabel "azimuth / °"
set ylabel "polar / °"
set output "analysis/spherical_abs.png"
splot "analysis/abs_spherical_symmetry.txt" u 1:2:3  notitle

#set xlabel "t_{int}/fs"
#set ylabel "E_{loss}/eV"
#set output "analysis/eloss_tint_all.png"
#splot "analysis/eloss_tint_all.txt" u 1:2:3 notitle

#set output "analysis/eloss_tint_inplane.png"
#splot "analysis/eloss_tint_ip.txt" u 1:2:3 notitle

set xlabel "{/cmti10 E}_{loss}/eV"
set ylabel "{/cmti10 r}_{sp}/Å"
set output "analysis/eloss_psd.png"
set title "Total energy loss vs closest proj-surf approach"
splot "analysis/eloss_psd.txt" notitle

set output "analysis/eloss_psd_in_plane.png"
set title "In plane energy loss vs closest proj-surf approach"
splot "analysis/eloss_psd_in_plane.txt" notitle


set xlabel "{/cmti10 r}_{sp}/Å"
set ylabel "Scattering angle"
set yrange [] reverse writeback
set output "analysis/polar_psd.png"
set title "In plane scattering angle vs closest proj-surf approach"
splot "analysis/polar_psd.txt" notitle


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

# points plot
reset
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
