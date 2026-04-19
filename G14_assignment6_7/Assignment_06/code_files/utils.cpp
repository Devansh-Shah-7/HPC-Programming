#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "utils.h"

/*
 * Parallel bilinear scatter interpolation using mesh privatisation.
 * Strategy: each thread maintains a private copy of the mesh,
 * scatters its assigned particles, then reduces into the global mesh.
 * This eliminates all race conditions without atomic operations.
 */
void interpolation(double *mesh_value, Points *points) {
    int total_nodes = GRID_X * GRID_Y;
    memset(mesh_value, 0, total_nodes * sizeof(double));

    int max_threads = omp_get_max_threads();
    double **private_mesh = (double **)malloc(max_threads * sizeof(double *));
    for (int t = 0; t < max_threads; t++)
        private_mesh[t] = (double *)calloc(total_nodes, sizeof(double));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double *local_mesh = private_mesh[tid];

        #pragma omp for schedule(static)
        for (int p = 0; p < NUM_Points; p++) {
            double x = points[p].x;
            double y = points[p].y;
            int col = (int)(x / dx);
            int row = (int)(y / dy);
            if (col >= NX) col = NX - 1;
            if (row >= NY) row = NY - 1;
            double lx = x - col * dx;
            double ly = y - row * dy;
            double w00 = (dx - lx) * (dy - ly);
            double w10 = lx        * (dy - ly);
            double w01 = (dx - lx) * ly;
            double w11 = lx        * ly;
            local_mesh[ row      * GRID_X + col    ] += w00;
            local_mesh[ row      * GRID_X + col + 1] += w10;
            local_mesh[(row + 1) * GRID_X + col    ] += w01;
            local_mesh[(row + 1) * GRID_X + col + 1] += w11;
        }

        // Reduction phase: each thread adds its slab of local mesh to global
        #pragma omp for schedule(static)
        for (int i = 0; i < total_nodes; i++) {
            double sum = 0.0;
            for (int t = 0; t < max_threads; t++) sum += private_mesh[t][i];
            mesh_value[i] = sum;
        }
    }

    for (int t = 0; t < max_threads; t++) free(private_mesh[t]);
    free(private_mesh);
}

void save_mesh(double *mesh_value) {
    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) { printf("Error creating Mesh.out\n"); exit(1); }
    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++)
            fprintf(fd, "%lf ", mesh_value[i * GRID_X + j]);
        fprintf(fd, "\n");
    }
    fclose(fd);
}
