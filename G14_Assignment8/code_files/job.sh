#!/bin/bash
#SBATCH --job-name=hpc_assign8
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=output_%j.txt

# Load MPI module (adjust name if needed on your cluster)
module load openmpi

# Set OMP threads = cores per node
export OMP_NUM_THREADS=16

echo "Running MPI+OpenMP job"
echo "Nodes: $SLURM_NNODES | Tasks: $SLURM_NTASKS | OMP Threads: $OMP_NUM_THREADS"

mpirun -np 4 ./main_mpi.out input.bin
