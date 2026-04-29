#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include "init.h"
#include "utils.h"

// Global variables
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {

    // ── MPI Init ─────────────────────────────────────────────────────────
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) printf("Usage: %s <input_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    // ── Rank 0 reads header ───────────────────────────────────────────────
    if (rank == 0) {
        FILE *file = fopen(argv[1], "rb");
        if (!file) {
            printf("Error opening input file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        (void)fread(&NX,         sizeof(int), 1, file);
        (void)fread(&NY,         sizeof(int), 1, file);
        (void)fread(&NUM_Points, sizeof(int), 1, file);
        (void)fread(&Maxiter,    sizeof(int), 1, file);
        fclose(file);
    }

    // ── Broadcast parameters to all ranks ────────────────────────────────
    MPI_Bcast(&NX,         1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NY,         1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NUM_Points, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Maxiter,    1, MPI_INT, 0, MPI_COMM_WORLD);

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    if (rank == 0) {
        printf("Grid: %d x %d | Points: %d | Iterations: %d\n",
               NX, NY, NUM_Points, Maxiter);
        printf("MPI ranks: %d | OMP threads per rank: %d\n",
               size, omp_get_max_threads());
    }

    // ── Particle decomposition ───────────────────────────────────────────
    int local_start = (rank * NUM_Points) / size;
    int local_end   = ((rank + 1) * NUM_Points) / size;
    int local_count = local_end - local_start;

    Points *points    = (Points *)calloc(NUM_Points,  sizeof(Points));
    Points *local_pts = (Points *)calloc(local_count, sizeof(Points));

    double *mesh_value       = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
    double *local_mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));

    double total_interp = 0.0, total_mover = 0.0;

    // ── Open file on rank 0 for iteration reading ─────────────────────────
    FILE *file = NULL;
    if (rank == 0) {
        file = fopen(argv[1], "rb");
        fseek(file, 4 * sizeof(int), SEEK_SET);
    }

    for (int iter = 0; iter < Maxiter; iter++) {

        // Rank 0 reads all points
        if (rank == 0) {
            for (int i = 0; i < NUM_Points; i++) {
                (void)fread(&points[i].x, sizeof(double), 1, file);
                (void)fread(&points[i].y, sizeof(double), 1, file);
                points[i].active = 1;
            }
        }

        // Broadcast all points to every rank
        MPI_Bcast(points, NUM_Points * sizeof(Points), MPI_BYTE,
                  0, MPI_COMM_WORLD);

        // Copy this rank's chunk to local_pts
        for (int i = 0; i < local_count; i++)
            local_pts[i] = points[local_start + i];

        // ── Interpolation ─────────────────────────────────────────────────
        double t0 = MPI_Wtime();
        interpolation_mpi(local_mesh_value, local_pts, local_count);
        double t1 = MPI_Wtime();
        total_interp += t1 - t0;

        // ── Reduce all local meshes into global mesh ───────────────────────
        MPI_Allreduce(local_mesh_value, mesh_value,
                      GRID_X * GRID_Y, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);

        // ── Normalize ─────────────────────────────────────────────────────
        double vmin, vmax;
        normalize_mesh(mesh_value, &vmin, &vmax);

        // ── Mover ─────────────────────────────────────────────────────────
        double t2 = MPI_Wtime();
        mover_mpi(mesh_value, local_pts, local_count);
        double t3 = MPI_Wtime();
        total_mover += t3 - t2;

        // ── Denormalize ───────────────────────────────────────────────────
        denormalize_mesh(mesh_value, vmin, vmax);

        // Count active particles
        int local_active = 0;
        for (int p = 0; p < local_count; p++)
            if (local_pts[p].active) local_active++;

        int global_active = 0;
        MPI_Reduce(&local_active, &global_active, 1, MPI_INT,
                   MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("Iter %2d | interp: %.4fs | mover: %.4fs | active: %d/%d\n",
                   iter, t1-t0, t3-t2, global_active, NUM_Points);
        }
    }

    // ── Save mesh (rank 0 only) ───────────────────────────────────────────
    if (rank == 0) {
        save_mesh(mesh_value);
        printf("\n========== Timing Summary (Rank 0) ==========\n");
        printf("Total interpolation time : %.6f s\n", total_interp);
        printf("Total mover time         : %.6f s\n", total_mover);
        printf("Total pipeline time      : %.6f s\n", total_interp + total_mover);
        printf("=============================================\n");
        fclose(file);
    }

    free(points);
    free(local_pts);
    free(mesh_value);
    free(local_mesh_value);

    MPI_Finalize();
    return 0;
}
