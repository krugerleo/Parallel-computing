#!/bin/sh
# tee -a file.txt >/dev/null
now="$(date +'%d-%m-%Y-%H-%M')"

#!/bin/bash
echo "Size (SIZE);Parallel Time;Execution" >> ./resultados/result$now.txt
for x in `seq 20`
do
    mpirun -np 6 kruger | tee -a ./resultados/result$now.txt >/dev/null
    echo "$x" >> ./resultados/result$now.txt
done
