#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "utils.h"

// ── MPI Interpolation: local particles → local mesh ──────────────────────────
void interpolation_mpi(double *local_mesh, Points *local_pts, int count) {

    int total = GRID_X * GRID_Y;
    memset(local_mesh, 0, total * sizeof(double));

    #pragma omp parallel
    {
        double *thread_mesh = (double *)calloc(total, sizeof(double));

        #pragma omp for schedule(static)
        for (int p = 0; p < count; p++) {

            if (!local_pts[p].active) continue;

            double x = local_pts[p].x;
            double y = local_pts[p].y;

            int i = (int)(x / dx);
            int j = (int)(y / dy);

            if (i >= NX) i = NX - 1;
            if (j >= NY) j = NY - 1;

            double lx = x - i * dx;
            double ly = y - j * dy;

            double w_ij   = (dx - lx) * (dy - ly);
            double w_i1j  = lx        * (dy - ly);
            double w_ij1  = (dx - lx) * ly;
            double w_i1j1 = lx        * ly;

            thread_mesh[ j    * GRID_X + i     ] += w_ij;
            thread_mesh[ j    * GRID_X + (i+1) ] += w_i1j;
            thread_mesh[(j+1) * GRID_X + i     ] += w_ij1;
            thread_mesh[(j+1) * GRID_X + (i+1) ] += w_i1j1;
        }

        #pragma omp critical
        {
            for (int k = 0; k < total; k++)
                local_mesh[k] += thread_mesh[k];
        }

        free(thread_mesh);
    }
}

// ── MPI Mover: global mesh → local particles ─────────────────────────────────
void mover_mpi(double *mesh_value, Points *local_pts, int count) {

    #pragma omp parallel for schedule(static)
    for (int p = 0; p < count; p++) {

        if (!local_pts[p].active) continue;

        double x = local_pts[p].x;
        double y = local_pts[p].y;

        int i = (int)(x / dx);
        int j = (int)(y / dy);

        if (i >= NX) i = NX - 1;
        if (j >= NY) j = NY - 1;

        double lx = x - i * dx;
        double ly = y - j * dy;

        double w_ij   = (dx - lx) * (dy - ly);
        double w_i1j  = lx        * (dy - ly);
        double w_ij1  = (dx - lx) * ly;
        double w_i1j1 = lx        * ly;

        double Fi =  w_ij   * mesh_value[ j    * GRID_X + i     ]
                   + w_i1j  * mesh_value[ j    * GRID_X + (i+1) ]
                   + w_ij1  * mesh_value[(j+1) * GRID_X + i     ]
                   + w_i1j1 * mesh_value[(j+1) * GRID_X + (i+1) ];

        double cell_area = dx * dy;
        if (cell_area > 1e-15) Fi /= cell_area;

        double x_new = x + Fi * dx;
        double y_new = y + Fi * dy;

        if (x_new < 0.0 || x_new > 1.0 || y_new < 0.0 || y_new > 1.0)
            local_pts[p].active = 0;
        else {
            local_pts[p].x = x_new;
            local_pts[p].y = y_new;
        }
    }
}

// ── Normalize mesh to [-1, 1] ────────────────────────────────────────────────
void normalize_mesh(double *mesh_value, double *out_min, double *out_max) {

    double vmin =  1e18, vmax = -1e18;

    #pragma omp parallel for reduction(min:vmin) reduction(max:vmax)
    for (int k = 0; k < GRID_X * GRID_Y; k++) {
        if (mesh_value[k] < vmin) vmin = mesh_value[k];
        if (mesh_value[k] > vmax) vmax = mesh_value[k];
    }

    double range = vmax - vmin;
    if (range < 1e-15) range = 1.0;

    #pragma omp parallel for schedule(static)
    for (int k = 0; k < GRID_X * GRID_Y; k++)
        mesh_value[k] = 2.0 * (mesh_value[k] - vmin) / range - 1.0;

    *out_min = vmin;
    *out_max = vmax;
}

// ── Denormalize mesh ─────────────────────────────────────────────────────────
void denormalize_mesh(double *mesh_value, double vmin, double vmax) {

    double range = vmax - vmin;
    if (range < 1e-15) range = 1.0;

    #pragma omp parallel for schedule(static)
    for (int k = 0; k < GRID_X * GRID_Y; k++)
        mesh_value[k] = (mesh_value[k] + 1.0) / 2.0 * range + vmin;
}

// ── Save mesh to file ────────────────────────────────────────────────────────
void save_mesh(double *mesh_value) {

    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) { printf("Error creating Mesh.out\n"); exit(1); }

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++) {
            if (j > 0) fprintf(fd, " ");
            fprintf(fd, "%lf", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }
    fclose(fd);
}
