#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "utils.h"

/*
 * OpenMP parallel bilinear scatter-mover pipeline.
 *
 * SCATTER: privatisation per thread, then #pragma omp critical reduction.
 *   - Each thread writes to its own private mesh (no races).
 *   - After all scatter is done (implicit barrier after omp for), each thread
 *     adds its full private mesh to global inside a critical section.
 *   - Critical overhead = nthreads × O(GRID_SIZE) work – tiny vs scatter.
 *
 * NORMALIZE / DENORM: parallel loops, parallel max-reduction.
 * MOVER: embarrassingly parallel (grid read-only, each particle owns itself).
 */
void interpolation_omp(double *mesh_value, Points *points) {

    int total = GRID_X * GRID_Y;

    /* ---- 1. Reset mesh -------------------------------------------------- */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < total; i++)
        mesh_value[i] = 0.0;

    /* ---- 2. SCATTER: point -> mesh (privatisation + critical reduction) -- */
    #pragma omp parallel
    {
        double *local_mesh = (double *) calloc(total, sizeof(double));
        if (!local_mesh) { fprintf(stderr,"OOM\n"); exit(1); }

        /* Each thread scatters its share of particles (implicit barrier after). */
        #pragma omp for schedule(static)
        for (int p = 0; p < NUM_Points; p++) {
            if (!points[p].active) continue;

            double x = points[p].x;
            double y = points[p].y;

            int ix = (int)(x / dx);
            int iy = (int)(y / dy);
            if (ix >= NX) ix = NX - 1;
            if (iy >= NY) iy = NY - 1;

            double lx = x - ix * dx;
            double ly = y - iy * dy;

            double w00 = (dx - lx) * (dy - ly);
            double w10 =        lx * (dy - ly);
            double w01 = (dx - lx) *        ly;
            double w11 =        lx *        ly;

            local_mesh[ iy      * GRID_X + ix    ] += w00;
            local_mesh[ iy      * GRID_X + ix + 1] += w10;
            local_mesh[(iy + 1) * GRID_X + ix    ] += w01;
            local_mesh[(iy + 1) * GRID_X + ix + 1] += w11;
        }
        /* Implicit barrier after the omp for: all scatter done. */

        /* Add this thread's full local mesh to global (serialised, but O(GRID) not O(N)). */
        #pragma omp critical
        {
            for (int i = 0; i < total; i++)
                mesh_value[i] += local_mesh[i];
        }

        free(local_mesh);

        /* Wait for all threads to finish reducing before normalise. */
        #pragma omp barrier
    }

    /* ---- 3. NORMALIZE to [0, 1] ----------------------------------------- */
    double max_val = 0.0;
    #pragma omp parallel for reduction(max:max_val) schedule(static)
    for (int i = 0; i < total; i++)
        if (mesh_value[i] > max_val) max_val = mesh_value[i];

    double inv_max = (max_val > 1e-15) ? 1.0 / max_val : 0.0;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < total; i++)
        mesh_value[i] *= inv_max;

    /* ---- 4. MOVER: mesh -> point (embarrassingly parallel) --------------- */
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int p = 0; p < NUM_Points; p++) {
        if (!points[p].active) continue;

        double x = points[p].x;
        double y = points[p].y;

        int ix = (int)(x / dx);
        int iy = (int)(y / dy);
        if (ix >= NX) ix = NX - 1;
        if (iy >= NY) iy = NY - 1;

        double lx = x - ix * dx;
        double ly = y - iy * dy;

        double w00 = (dx - lx) * (dy - ly);
        double w10 =        lx * (dy - ly);
        double w01 = (dx - lx) *        ly;
        double w11 =        lx *        ly;

        double Fi = w00 * mesh_value[ iy      * GRID_X + ix    ]
                  + w10 * mesh_value[ iy      * GRID_X + ix + 1]
                  + w01 * mesh_value[(iy + 1) * GRID_X + ix    ]
                  + w11 * mesh_value[(iy + 1) * GRID_X + ix + 1];

        double new_x = x + Fi * dx;
        double new_y = y + Fi * dy;

        if (new_x < 0.0 || new_x > 1.0 || new_y < 0.0 || new_y > 1.0)
            points[p].active = 0;
        else { points[p].x = new_x; points[p].y = new_y; }
    }

    /* ---- 5. DENORMALIZE -------------------------------------------------- */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < total; i++)
        mesh_value[i] *= max_val;
}
