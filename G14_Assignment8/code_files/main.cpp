#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

// Global variables
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {

    if (argc != 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "rb");
    if (!file) { printf("Error opening input file\n"); exit(1); }

    (void)fread(&NX,         sizeof(int), 1, file);
    (void)fread(&NY,         sizeof(int), 1, file);
    (void)fread(&NUM_Points, sizeof(int), 1, file);
    (void)fread(&Maxiter,    sizeof(int), 1, file);

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    printf("Grid: %d x %d | Points: %d | Iterations: %d | Threads: %d\n",
           NX, NY, NUM_Points, Maxiter, omp_get_max_threads());

    double *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points     = (Points *)calloc(NUM_Points, sizeof(Points));

    double total_interp = 0.0, total_mover = 0.0;

    for (int iter = 0; iter < Maxiter; iter++) {

        read_points(file, points);

        double t0 = omp_get_wtime();
        interpolation(mesh_value, points);
        double t1 = omp_get_wtime();
        total_interp += t1 - t0;

        double vmin, vmax;
        normalize_mesh(mesh_value, &vmin, &vmax);

        double t2 = omp_get_wtime();
        mover(mesh_value, points);
        double t3 = omp_get_wtime();
        total_mover += t3 - t2;

        denormalize_mesh(mesh_value, vmin, vmax);

        int active = 0;
        for (int p = 0; p < NUM_Points; p++)
            if (points[p].active) active++;

        printf("Iter %2d | interp: %.4fs | mover: %.4fs | active: %d/%d\n",
               iter, t1-t0, t3-t2, active, NUM_Points);
    }

    save_mesh(mesh_value);

    printf("\n========== Timing Summary ==========\n");
    printf("Total interpolation time : %.6f s\n", total_interp);
    printf("Total mover time         : %.6f s\n", total_mover);
    printf("Total pipeline time      : %.6f s\n", total_interp + total_mover);
    printf("=====================================\n");

    free(mesh_value);
    free(points);
    fclose(file);
    return 0;
}
