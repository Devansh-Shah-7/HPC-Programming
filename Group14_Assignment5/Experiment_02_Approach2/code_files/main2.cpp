#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

// Global variables
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {

    // ── Experiment 2 Settings ──
    NX       = 500;           // Fixed grid
    NY       = 200;
    Maxiter  = 10;
    NUM_Points = 14000000;    // 14 million particles
    int approach = 2;         // 1 = Deferred, 2 = Immediate (change this)
    int num_threads = 16;      // 1, 2, 4, 8, 16 (change this)

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    omp_set_num_threads(num_threads);

    double *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points     = (Points *)calloc(NUM_Points, sizeof(Points));

    // Initialize ONCE outside the loop
    initializepoints(points);

    double total_interp = 0.0, total_mover = 0.0;

    printf("Iter\tInterp\t\tMover\t\tTotal\t\tDeleted\n");

    for (int iter = 0; iter < Maxiter; iter++) {
        long deleted = 0;

        double t0 = omp_get_wtime();
        interpolation(mesh_value, points);
        double t1 = omp_get_wtime();

        if (approach == 1)
            mover_deferred(points, dx, dy, &deleted);
        else
            mover_immediate(points, dx, dy, &deleted);
        double t2 = omp_get_wtime();

        double ti = t1 - t0;
        double tm = t2 - t1;
        total_interp += ti;
        total_mover  += tm;

        printf("%d\t%.4f\t\t%.4f\t\t%.4f\t\t%ld\n",
               iter+1, ti, tm, ti+tm, deleted);
    }

    printf("\nTOTAL interp: %.4f s\n", total_interp);
    printf("TOTAL mover:  %.4f s\n",  total_mover);
    printf("TOTAL:        %.4f s\n",  total_interp + total_mover);

    // CSV line for easy copy-paste into spreadsheet
    printf("\nCSV: %d,%d,%d,%d,%d,%.6f,%.6f,%.6f\n",
           NX, NY, NUM_Points, approach, num_threads,
           total_interp, total_mover, total_interp + total_mover);

    free(mesh_value);
    free(points);
    return 0;
}