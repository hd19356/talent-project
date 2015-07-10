set terminal postscript enhanced eps
#
set output "exdiag.eps"
#
set xlabel "g"
set ylabel "g.s. energy"
#
plot "diag.out" u 1:3 w l lt  0 title "reference energy",\
     "diag.out" u 1:2 w l lt -1 title "exact FCI"

#
#
#
set output "mbpte2.eps"
#
set xlabel "g"
set ylabel "correlation energy"
#
plot "mbpt.out" u 1:2 w l lt -1 title "{/Symbol D}E^{(2)}",\
     "diag.out" u 1:(($2)-($3)) w l lt 0 title "FCI"
#
#
#
set output "gs-g.eps"
#
set xlabel "g"
set ylabel "g.s. structure"
set yrange [-0.1:1.1]
set key center right
#
plot "diag.out" u 1:(($4)**2) w l lt -1 title "4h state",\
     "diag.out" u 1:(($5)**2) w l lt  0 title "2h-2p state",\
     "diag.out" u 1:(($6)**2) w l lt  1 title "2h-2p state",\
     "diag.out" u 1:(($7)**2) w l lt  2 title "2h-2p state",\
     "diag.out" u 1:(($8)**2) w l lt  3 title "2h-2p state",\
     "diag.out" u 1:(($9)**2) w l lt  4 title "4p state"
#
#
#
set output
set terminal unknown
reset
