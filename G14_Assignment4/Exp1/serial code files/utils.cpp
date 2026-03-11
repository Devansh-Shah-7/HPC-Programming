#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "utils.h"

// Interpolation (Serial Code)
void interpolation(double *mesh_value, Points *points) {
    // Zero out the mesh first
    memset(mesh_value, 0, GRID_X * GRID_Y * sizeof(double));

    for (int i = 0; i < NUM_Points; i++) {
        // Find lower-left grid index
        int ix = (int)(points[i].x / dx);
        int iy = (int)(points[i].y / dy);

        // Clamp to valid range
        if (ix >= NX) ix = NX - 1;
        if (iy >= NY) iy = NY - 1;

        // Fractional offsets (weights)
        double fx = (points[i].x / dx) - ix;
        double fy = (points[i].y / dy) - iy;

        // Bilinear scatter to 4 corners
        mesh_value[ iy      * GRID_X + ix    ] += (1.0 - fx) * (1.0 - fy);
        mesh_value[ iy      * GRID_X + (ix+1)] +=        fx  * (1.0 - fy);
        mesh_value[(iy+1)   * GRID_X + ix    ] += (1.0 - fx) *        fy;
        mesh_value[(iy+1)   * GRID_X + (ix+1)] +=        fx  *        fy;
    }
}

// Stochastic Mover (Serial Code) 
void mover_serial(Points *points, double deltaX, double deltaY) {
    for (int i = 0; i < NUM_Points; i++) {
        double newx, newy;
        do {
            double rx = ((double)rand() / RAND_MAX) * 2.0 - 1.0; // [-1, 1]
            double ry = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            newx = points[i].x + rx * deltaX;
            newy = points[i].y + ry * deltaY;
        } while (newx < 0.0 || newx > 1.0 || newy < 0.0 || newy > 1.0);

        points[i].x = newx;
        points[i].y = newy;
    }
}

// Stochastic Mover (Parallel Code) 
void mover_parallel(Points *points, double deltaX, double deltaY) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < NUM_Points; i++) {
        unsigned int seed = omp_get_thread_num() * 1234 + i;
        double newx, newy;
        do {
            double rx = ((double)rand_r(&seed) / RAND_MAX) * 2.0 - 1.0;
            double ry = ((double)rand_r(&seed) / RAND_MAX) * 2.0 - 1.0;
            newx = points[i].x + rx * deltaX;
            newy = points[i].y + ry * deltaY;
        } while (newx < 0.0 || newx > 1.0 || newy < 0.0 || newy > 1.0);

        points[i].x = newx;
        points[i].y = newy;
    }
}

// Write mesh to file
void save_mesh(double *mesh_value) {

    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) {
        printf("Error creating Mesh.out\n");
        exit(1);
    }

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++) {
            fprintf(fd, "%lf ", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }

    fclose(fd);
}