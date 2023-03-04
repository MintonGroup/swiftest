#!/bin/zsh -l
#SBATCH -A standby
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=swifter-16k
#SBATCH --mem-per-cpu=3072
#SBATCH -o %x.out
#SBATCH -e %x.err

cd $SLURM_SUBMIT_DIR
filename="${SLURM_SUBMIT_DIR}/swifter-16k-timehist.log"
rm ${filename}
echo "N cores,wall time(s)" > ${filename}

export KMP_STACKSIZE=2G
count=1
while [ $count -lt 25 ]
do
  n=`printf %03d $count`
  export OMP_NUM_THREADS=${count}
  ./swifter_symba_omp < start.in > term${n}.out
  cat term${n}.out | grep complete | awk -v NTHREADS=$OMP_NUM_THREADS '{print NTHREADS,",",$7}' >> ${filename}
  count=`expr $count + 1`
done

