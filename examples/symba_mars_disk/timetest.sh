#!/bin/bash
outname="timetest.csv"
echo "N cores, wall time (s)" > $outname
for value in {1..12}
do
   export OMP_NUM_THREADS=$value
   walltime="$(/usr/bin/time -f %e ./swiftest_driver symba param.in 2>&1 > /dev/null)"
   echo "$value,$walltime" >> $outname
done