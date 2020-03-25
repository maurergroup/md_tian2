#!/usr/bin/gnuplot

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
