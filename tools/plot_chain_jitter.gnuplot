# Recorded C++ SDK trajectories; plot only, no new physics calculation.
# From repo root: gnuplot tools/plot_chain_jitter.gnuplot
# Publication wrapper: sh tools/plot_chain_jitter.sh (also builds the diagram and normalizes SVG).
# Optional -e "DATA='...'; OUT='...'; FORMAT='png'" for another data/output directory.
if (!exists("DATA")) DATA = "docs/book/evidence/chain-jitter-20261001"
if (!exists("OUT")) OUT = "docs/images/chain-jitter-20261001"
if (!exists("FORMAT")) FORMAT = "svg"
set encoding utf8
if (FORMAT eq "png") {set terminal pngcairo size 1100,780 font "DejaVu Sans,11"} else {set terminal svg size 1100,780 font "DejaVu Sans,11" noenhanced}
set datafile separator comma
set grid back lc rgb "#dddddd"
set border 3
set tics nomirror
set key opaque top right
set style line 1 lc rgb "#2166ac" lw 1.5
set style line 2 lc rgb "#d6604d" lw 1.5
set style line 3 lc rgb "#1b7837" lw 1.5 dt 2
set xrange [10:30]
set xlabel "Время от запуска, с"
set ylabel "T, мкДж"
set output OUT."/rest-kinetic.".FORMAT
set multiplot layout 2,1 title "Нетронутая цепь: остаточное движение без сна — 2026-10-01"
set title "Поздний интервал 10–30 с"
set yrange [0:500]
plot DATA."/rest.csv" using 2:($3*1e6) with lines ls 1 title "SDK, 12 итераций", \
     DATA."/rp3-rest.csv" using 2:($3*1e6) with lines ls 2 title "Кандидат RP3D, 12 итераций", \
     DATA."/rp3-iter24-rest.csv" using 2:($3*1e6) with lines ls 3 title "Кандидат RP3D, 24 итерации"
set title "Увеличение интервала 20–30 с; те же CSV"
set xrange [20:30]
set yrange [0:130]
plot DATA."/rest.csv" using 2:($3*1e6) with lines ls 1 notitle, \
     DATA."/rp3-rest.csv" using 2:($3*1e6) with lines ls 2 notitle, \
     DATA."/rp3-iter24-rest.csv" using 2:($3*1e6) with lines ls 3 notitle
unset multiplot
unset output

set xrange [10:30]
set yrange [0:180]
set ylabel "max(ΔE, 0), мкДж / кадр"
set output OUT."/positive-energy.".FORMAT
set multiplot layout 2,1 title "Положительные приращения полной энергии; отрицательная часть исключена явно"
set title "Без захвата — отсчёт от запуска, сон выключен"
plot DATA."/rest.csv" using 2:($6>0?$6*1e6:0) with impulses ls 1 title "SDK, 12 итераций", \
     DATA."/rp3-rest.csv" using 2:($6>0?$6*1e6:0) with impulses ls 2 title "Кандидат RP3D, 12 итераций", \
     DATA."/rp3-iter24-rest.csv" using 2:($6>0?$6*1e6:0) with impulses ls 3 title "Кандидат RP3D, 24 итерации"
set title "После стандартного захвата — кольцо отпущено на 4-й секунде"
set label 1 "SDK 12 и кандидат RP3D 24: положительных ΔE после 10 с нет в этом прогоне" at graph 0.03,0.75 front
plot DATA."/pull.csv" using 2:($6>0?$6*1e6:0) with impulses ls 1 notitle, \
     DATA."/rp3-pull.csv" using 2:($6>0?$6*1e6:0) with impulses ls 2 notitle, \
     DATA."/rp3-iter24-pull.csv" using 2:($6>0?$6*1e6:0) with impulses ls 3 notitle
unset label 1
unset multiplot
unset output

# Sum pre-serialization deltas. Do not subtract rounded total_J or double-count the counterfactual.
set autoscale x
set autoscale y
F = DATA."/counterfactual-785.csv"
stats F using (strcol(2) eq "gravity-kick" ? $15 : 1/0) nooutput
G = STATS_sum
stats F using (strcol(2) eq "warm-start" ? $15 : 1/0) nooutput
W = STATS_sum
stats F using (strcol(2) eq "contact-iteration" ? $15 : 1/0) nooutput
C = STATS_sum
stats F using (strcol(2) eq "shock" ? $15 : 1/0) nooutput
S = STATS_sum
stats F using (strcol(2) eq "pose-drift" ? $15 : 1/0) nooutput
D = STATS_sum
stats F using (strcol(2) ne "split-counterfactual" ? $13 : 1/0) nooutput
T = STATS_sum
stats F using (strcol(2) ne "split-counterfactual" ? $14 : 1/0) nooutput
U = STATS_sum
stats F using (strcol(2) ne "split-counterfactual" ? $15 : 1/0) nooutput
E = STATS_sum
stats F using (strcol(2) eq "split-counterfactual" ? $14 : 1/0) nooutput
B = STATS_sum
set print $Stages
print sprintf('"Гравитация",%.12g', G*1e3)
print sprintf('"Warm start",%.12g', W*1e3)
print sprintf('"Контакты",%.12g', C*1e3)
print sprintf('"Shock",%.12g', S*1e3)
print sprintf('"Перенос поз",%.12g', D*1e3)
unset print
set print $Net
print sprintf('"ΔT",%.12g', T*1e6)
print sprintf('"ΔU",%.12g', U*1e6)
print sprintf('"ΔE",%.12g', E*1e6)
unset print
set output OUT."/frame-785-budget.".FORMAT
set multiplot layout 1,2 title "SDK: измеренный баланс кадра 785 (13.0833 с); 40/40 подшагов совпали"
unset key
unset xlabel
set style fill solid 0.75 border -1
set boxwidth 0.65
set xrange [-0.7:4.7]
set yrange [-14:10]
set ylabel "Изменение E за стадию, мДж"
set title "Промежуточные изменения; warm start не является внешней работой"
set xtics rotate by -30
plot $Stages using 0:2:xticlabels(1) with boxes lc rgb "#2166ac"
set xrange [-0.7:2.7]
set yrange [-40:240]
set ylabel "Изменение за весь кадр, мкДж"
set title "Итог: рост через потенциальную энергию"
set xtics norotate
set label 1 sprintf("Локальная разница split / без split:\nΔU = +%.5f мДж\nне прибавляется повторно к балансу", B*1e3) at graph 0.03,0.95 front
plot $Net using 0:2:xticlabels(1) with boxes lc rgb "#d6604d"
unset label 1
unset multiplot
unset output
set print OUT."/plot-values.txt"
print sprintf("frame785_stage_delta_E_J=%.12g", E)
print sprintf("frame785_split_local_delta_U_J=%.12g", B)
unset print
