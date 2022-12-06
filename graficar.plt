set term gif animate delay 2
set output "viaje.gif"

stats "trayectoria_nave.txt"

set style data points

set xr[-1.5:1.5]
set yr[-1.5:1.5]
set object circle at first 0,0 radius 0.05
do for [i=0:floor(1000000.0/997)-1] {
    plot "trayectoria_nave.txt" u 1:2 every ::i::i t "Nave" pt 25 ps 1 lc "blue", "trayectoria_luna.txt" u 1:2 every ::i::i t "Luna" pt 7 ps 2 lc "dark-grey"
}
