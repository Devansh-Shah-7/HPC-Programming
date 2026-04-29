# HPC Assignment 8 — MPI + OpenMP Hybrid Parallelization

## Files
- `main.cpp`           → OpenMP only version (for local testing)
- `main_mpi.cpp`       → MPI + OpenMP hybrid version (for cluster)
- `utils.cpp`          → OpenMP utility functions
- `utils_mpi.cpp`      → MPI + OpenMP utility functions
- `utils.h`            → Header for all utility functions
- `init.cpp / init.h`  → Particle initialization and I/O
- `input_file_maker.cpp` → Generates binary input files
- `job.sh`             → SLURM job submission script

---

## Step 1: Compile everything

### OpenMP version (local machine):
```bash
g++ main.cpp utils.cpp init.cpp -lm -O2 -fopenmp -o main_omp.out
```

### MPI + OpenMP version:
```bash
mpicxx main_mpi.cpp utils_mpi.cpp init.cpp -lm -O2 -fopenmp -o main_mpi.out
```

### Input file maker:
```bash
g++ input_file_maker.cpp -o input_maker.out
```

---

## Step 2: Generate input files

```bash
./input_maker.out
# Config a: NX=250  NY=100  Points=900000    Iter=10
# Config b: NX=250  NY=100  Points=5000000   Iter=10
# Config c: NX=500  NY=200  Points=3600000   Iter=10
# Config d: NX=500  NY=200  Points=20000000  Iter=10
# Config e: NX=1000 NY=400  Points=14000000  Iter=10
```

---

## Step 3: Run OpenMP locally

```bash
OMP_NUM_THREADS=1  ./main_omp.out input.bin
OMP_NUM_THREADS=2  ./main_omp.out input.bin
OMP_NUM_THREADS=4  ./main_omp.out input.bin
OMP_NUM_THREADS=8  ./main_omp.out input.bin
OMP_NUM_THREADS=12 ./main_omp.out input.bin
```

---

## Step 4: Run MPI locally (for testing)

```bash
mpirun -np 1 ./main_mpi.out input.bin
mpirun -np 2 ./main_mpi.out input.bin
mpirun -np 4 ./main_mpi.out input.bin
```

---

## Step 5: Copy to cluster

Run this on YOUR LOCAL machine:
```bash
scp -r hpc_assign8_cluster/ username@gics1.daiict.ac.in:~/
```

---

## Step 6: On the cluster

```bash
ssh username@gics1.daiict.ac.in
cd hpc_assign8_cluster

# Compile
mpicxx main_mpi.cpp utils_mpi.cpp init.cpp -lm -O2 -fopenmp -o main_mpi.out

# Generate input files on cluster
g++ input_file_maker.cpp -o input_maker.out
./input_maker.out

# Submit job
sbatch job.sh

# Check status
squeue -u username

# View output when done
cat output_*.txt
```

---

## Verify correctness

```bash
./main_omp.out Test_input.bin
diff Mesh.out Test_Mesh.out   # should print nothing
```
