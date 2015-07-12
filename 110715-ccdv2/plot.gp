set terminal postscript enhanced
set output "deccd-g.eps"
#
set xlabel "g"
set ylabel "{/Symbol D}E_{CCD}"
#
plot "lccd.out" u 1:2 w l lt -1 notitle
#
#
#
set output "tccd-g.eps"
#
set xlabel "g"
set ylabel "t_@{ij}^{ab}"
set key top left
set key font ",10"
#
plot "lccd.out" u 1:3 w l lt -1 title "t_@{12}^{56}",\
     "lccd.out" u 1:4 w l lt  0 title "t_@{34}^{56}",\
     "lccd.out" u 1:5 w l lt  1 title "t_@{12}^{78}",\
     "lccd.out" u 1:6 w l lt  2 title "t_@{34}^{78}"
#
set output
set terminal unknown
reset
