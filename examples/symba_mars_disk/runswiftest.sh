#!/bin/zsh -l
#SBATCH -A daminton
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=14-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --job-name=high_high
#SBATCH --no-requeue
#SBATCH -o %x.out
#SBATCH -e %x.err
cd $SLURM_SUBMIT_DIR

module load utilities monitor
#module load gcc/10.2.0
#netcdf-fortran/4.5.3

# track per-code CPU load
monitor cpu percent --all-cores >cpu-percent.log &
CPU_PID=$!

# track memory usage
monitor cpu memory >cpu-memory.log &
MEM_PID=$!


export OMP_NUM_THREADS=12
export KMP_STACKSIZE=2G
./swiftest_driver symba param.in

# shut down the resource monitors
kill -s INT $CPU_PID $MEM_PID

