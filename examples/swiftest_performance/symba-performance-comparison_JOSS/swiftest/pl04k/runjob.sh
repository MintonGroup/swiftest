#!/bin/zsh -l
#SBATCH -A standby
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=swiftest-4k
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3072
#SBATCH -o %x.out
#SBATCH -e %x.err

cd $SLURM_SUBMIT_DIR
filename="${SLURM_SUBMIT_DIR}/swiftest-4k-timehist.log"
rm ${filename}
echo "N cores,wall time(s)" > ${filename}

count=1
while [ $count -lt 25 ]
do
  n=`printf %03d $count`
  export OMP_NUM_THREADS=${count}
  ./swiftest_driver symba param.in > term${n}.out
  grep Integration term${n}.out | tail -n1 | awk -v NTHREADS=$OMP_NUM_THREADS '{print NTHREADS,",",$6}' | sed 's/.$//' >> ${filename} 
  count=`expr $count + 1`
done

