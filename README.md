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

for i in {50,100,200,500,1000}
do
        for j in {1..10}
        do
                /DATA6/leon/Tov-diversity/simulate/simulateZW.pl ${i} ${j} 1E-5 1E-3 &>ZWsimulate-5-3_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateZW.pl ${i} ${j} 1E-3 1E-5 &>ZWsimulate-3-5_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateZW.pl ${i} ${j} 1E-8 1E-6 &>ZWsimulate-8-6_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateZW.pl ${i} ${j} 1E-6 1E-8 &>ZWsimulate-6-8_${i}_${j}.log &
        done
done

for i in {50,100,200,500,1000}
do
        for j in {1..10}
        do
                /DATA6/leon/Tov-diversity/simulate/simulateZ.pl ${i} ${j} 1E-5 1E-3 &>Zsimulate-5-3_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateZ.pl ${i} ${j} 1E-3 1E-5 &>Zsimulate-3-5_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateZ.pl ${i} ${j} 1E-8 1E-6 &>Zsimulate-8-6_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateZ.pl ${i} ${j} 1E-6 1E-8 &>Zsimulate-6-8_${i}_${j}.log &
        done
done


for i in {50,100,200,500,1000}
do
        for j in {1..10}
        do
                /DATA6/leon/Tov-diversity/simulate/simulateW.pl ${i} ${j} 1E-5 1E-3 &>Wsimulate-5-3_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateW.pl ${i} ${j} 1E-3 1E-5 &>Wsimulate-3-5_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateW.pl ${i} ${j} 1E-8 1E-6 &>Wsimulate-8-6_${i}_${j}.log &
                /DATA6/leon/Tov-diversity/simulate/simulateW.pl ${i} ${j} 1E-6 1E-8 &>Wsimulate-6-8_${i}_${j}.log &
        done
done


