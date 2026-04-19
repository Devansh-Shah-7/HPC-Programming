#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "init.h"
#include "utils.h"

// Global variables (same names as main.cpp)
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {

    if (argc < 2 || argc > 3) {
        printf("Usage: %s <input_file> [num_threads]\n", argv[0]);
        return 1;
    }

    /* Optional thread count argument */
    if (argc == 3) {
        int nthreads = atoi(argv[2]);
        omp_set_num_threads(nthreads);
    }

    printf("Running with %d OpenMP thread(s)\n", omp_get_max_threads());

    FILE *file = fopen(argv[1], "rb");
    if (!file) { printf("Error opening input file\n"); return 1; }

    fread(&NX,        sizeof(int), 1, file);
    fread(&NY,        sizeof(int), 1, file);
    fread(&NUM_Points,sizeof(int), 1, file);
    fread(&Maxiter,   sizeof(int), 1, file);

    GRID_X = NX + 1;  GRID_Y = NY + 1;
    dx = 1.0 / NX;    dy = 1.0 / NY;

    double  *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
    Points  *points     = (Points *)calloc(NUM_Points,       sizeof(Points));

    double total_time = 0.0;

    for (int iter = 0; iter < Maxiter; iter++) {
        read_points(file, points);

        double t0 = omp_get_wtime();
        interpolation_omp(mesh_value, points);
        double t1 = omp_get_wtime();

        total_time += t1 - t0;
    }

    save_mesh(mesh_value);
    printf("Total interpolation time (parallel) = %lf seconds\n", total_time);

    free(mesh_value);
    free(points);
    fclose(file);
    return 0;
}
