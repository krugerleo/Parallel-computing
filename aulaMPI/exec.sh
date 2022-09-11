#!/bin/sh
# tee -a file.txt >/dev/null
now="$(date +'%d-%m-%Y-%H-%M')"

#!/bin/bash
echo "Size (fileB,fileA,uniqAB);Serial Score;Serial LCS Time;Parallel Score;Parallel LCS Time;Execution" >> ./resultados/result$now.txt
for x in `seq 30`
do
    mpirun kruger entradas/file1.in entradas/file2.in entradas/uniqAB.in | tee -a ./resultados/result$now.txt >/dev/null
    echo "$x" >> ./resultados/result$now.txt
done
for y in `seq 30`
do
    mpirun kruger entradas/fileC.in entradas/fileD.in entradas/uniqAB.in | tee -a ./resultados/result$now.txt >/dev/null
    echo "$y" >> ./resultados/result$now.txt
done
for z in `seq 30`
do
    mpirun kruger entradas/fileA.in entradas/fileB.in entradas/uniqAB.in | tee -a ./resultados/result$now.txt >/dev/null
    echo "$z" >> ./resultados/result$now.txt
done
# for z in `seq 30`
# do
#     ./parallel 34 -wh | tee -a ./resultados/result$now.txt >/dev/null
#     echo "$z" >> ./resultados/result$now.txt
# done