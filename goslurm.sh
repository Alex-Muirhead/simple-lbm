#!/bin/bash -l
#SBATCH --job-name=CUDA_boost
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:10 # time (D-HH:MM)
#SBATCH --mem-per-cpu=4G # memory (MB)
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load hdf5

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo 'running from directory =' $CURRENT_DIR

if [ ! -f bin/main_gpu ] ; then
   echo "unable to find requested file"
   echo "you probably need to compile code"
   exit 2
fi

> $CURRENT_DIR/cpu_performance.txt
for i in {0..3}; do
   time ./bin/main_cpu "$CURRENT_DIR/config.h5" "$CURRENT_DIR/cpu_output.h5" >> $CURRENT_DIR/cpu_performance.txt
done

> $CURRENT_DIR/gpu_performance.txt
for i in {0..3}; do
   time ./bin/main_gpu "$CURRENT_DIR/config.h5" "$CURRENT_DIR/gpu_output.h5" >> $CURRENT_DIR/gpu_performance.txt
done