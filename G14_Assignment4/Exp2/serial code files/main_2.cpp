#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {

    // CHANGE THESE FOR EACH RUN
    NX = 1000;   // Config 1: 250 | Config 2: 500 | Config 3: 1000
    NY = 400;   // Config 1: 100 | Config 2: 200 | Config 3: 400

    Maxiter = 10;
    NUM_Points = 100000000;  // Fixed at 10^8 for all configs
    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    Points *points = (Points *) malloc(NUM_Points * sizeof(Points));
    double *mesh_value = (double *) calloc(GRID_X * GRID_Y, sizeof(double));

    initializepoints(points);  // Only ONCE

    double total_interp_time = 0.0;
    for (int iter = 0; iter < Maxiter; iter++) {
        clock_t s = clock();
        interpolation(mesh_value, points);
        clock_t e = clock();
        total_interp_time += (double)(e - s) / CLOCKS_PER_SEC;
    }

    printf("Config (NX=%d, NY=%d): Total Interp Time = %lf\n", NX, NY, total_interp_time);

    free(points);
    free(mesh_value);
    return 0;
}