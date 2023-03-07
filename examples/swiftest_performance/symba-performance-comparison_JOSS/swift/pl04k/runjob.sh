#!/bin/zsh -l
#SBATCH -A standby
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=swift-4k
#SBATCH -o %x.out
#SBATCH -e %x.err

cd $SLURM_SUBMIT_DIR
filename="${SLURM_SUBMIT_DIR}/swift-4k-timehist.log"
rm ${filename}
echo "N cores,wall time(s)" > ${filename}

export KMP_STACKSIZE=2G
count=1
while [ $count -lt 2 ]
do
  n=`printf %03d $count`
  export OMP_NUM_THREADS=${count}
  walltime="$(/usr/bin/time -f %e ./swift_symba5 < start.in 2>&1 > term$n.out)"
  echo "$count,$walltime" >> $filename
  count=`expr $count + 1`
done

