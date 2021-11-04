#!/bin/bash -l
#SBATCH --job-name=CUDA_boost
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:01 # time (D-HH:MM)
#SBATCH --mem-per-cpu=4G # memory (MB)
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo 'running with OMP_NUM_THREADS =' $OMP_NUM_THREADS
echo 'running with MKL_NUM_THREADS =' $MKL_NUM_THREADS
echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "running with OMP_NUM_THREADS= $OMP_NUM_THREADS "
echo "running with SLURM_TASKS_PER_NODE= $SLURM_TASKS_PER_NODE "

if [ ! -f bin/main_gpu ] ; then
   echo "unable to find requested file"
   echo "you probably need to compile code"
   exit 2
fi

time ./bin/main_cpu "taylor_green.h5"
time ./bin/main_gpu "taylor_green.h5"
