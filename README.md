# Turnover-of-sex-mutation

The simulation was performed in Perl scripts.

simulateZW.pl: N W haplotypes and 3N Z haplotypes

simulateW.pl: 4N W haplotypes

simulateZ.pl: 4N Z haploytpes

Command

i: size of N

j: lable of repeat

x: mutation rate

y: recombination rate

simulateZW.pl ${i} ${j} ${x} ${y} &>ZWsimulate-5-3_${i}_${j}.log &

simulateZ.pl ${i} ${j} ${x} ${y} &>Zsimulate-5-3_${i}_${j}.log &

simulateW.pl ${i} ${j} ${x} ${y} &>Wsimulate-5-3_${i}_${j}.log &


