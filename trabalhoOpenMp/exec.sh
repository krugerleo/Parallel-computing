#!/bin/sh
# tee -a file.txt >/dev/null
now="$(date +'%d-%m-%Y-%H-%M')"

#!/bin/bash
echo "Size (fileB,fileA,uniqAB);Serial Score;Serial LCS Time;Parallel Score;Parallel LCS Time;Execution" >> ./resultados/result$now.txt
for x in `seq 20`
do
    ./parallel CD -wh | tee -a ./resultados/result$now.txt >/dev/null
    echo "$x" >> ./resultados/result$now.txt
done
for y in `seq 20`
do
    ./parallel 12 -wh | tee -a ./resultados/result$now.txt >/dev/null
    echo "$y" >> ./resultados/result$now.txt
done
for z in `seq 2`
do
    ./parallel AB -wh | tee -a ./resultados/result$now.txt >/dev/null
    echo "$z" >> ./resultados/result$now.txt
done