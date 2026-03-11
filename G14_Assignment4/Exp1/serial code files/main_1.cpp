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

    // ---- CHANGE THESE FOR EACH CONFIG RUN ----
    NX = 1000;   // Config 1: 250 | Config 2: 500 | Config 3: 1000
    NY = 400;   // Config 1: 100 | Config 2: 200 | Config 3: 400
    // ------------------------------------------

    Maxiter = 10;
    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    // ---- PASTE THE EXPERIMENT 1 LOOP HERE ----
    int particle_counts[] = {100, 10000, 1000000, 100000000 };

    printf("NUM_Points\tTotal_Interp_Time\n");

    for (int pc = 0; pc < 5; pc++) {
        NUM_Points = particle_counts[pc];
        Points *points = (Points *) malloc(NUM_Points * sizeof(Points));
        double *mesh_value = (double *) calloc(GRID_X * GRID_Y, sizeof(double));

        double total_interp_time = 0.0;
        for (int iter = 0; iter < Maxiter; iter++) {
            initializepoints(points);

            clock_t s = clock();
            interpolation(mesh_value, points);
            clock_t e = clock();
            total_interp_time += (double)(e - s) / CLOCKS_PER_SEC;
        }
        printf("%d\t%lf\n", NUM_Points, total_interp_time);

        free(points);
        free(mesh_value);
    }
    // ------------------------------------------

    return 0;
}