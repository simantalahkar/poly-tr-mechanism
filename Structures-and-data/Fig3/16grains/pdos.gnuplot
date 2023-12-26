set term post enha colo 20
set out "pdos.eps"

set xlabel "frequency (THz)"
set ylabel "Phonon DOS"
unset key
plot "pdos.dat" u 1:2 w l
