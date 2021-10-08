#!/bin/zsh -l
#SBATCH -A standby
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=swiftest-2k
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=MaxMemPerCPU
#SBATCH -o %x.out
#SBATCH -e %x.err

cd $SLURM_SUBMIT_DIR
filename="${SLURM_SUBMIT_DIR}/swiftest-2k-timehist.log"
rm ${filename}
echo "N cores,wall time(s)" > ${filename}

count=1
while [ $count -lt 129 ]
do
  export OMP_NUM_THREADS=${count}
  ./swiftest_driver symba param.in > term${count}.out
  grep Integration term${count}.out | tail -n1 | awk -v NTHREADS=$OMP_NUM_THREADS '{print NTHREADS,",",$6}' | sed 's/.$//' >> ${filename} 
  mv interaction_timer.log interaction_timer${count}.log
  mv encounter_check_timer.log encounter_check_timer${count}.log
  count=`expr $count + 1`
done

