#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"

/*
 * Bilinear scatter-mover pipeline (serial):
 *   1. Reset mesh to zero.
 *   2. SCATTER  : each active particle distributes fi=1 to 4 surrounding
 *                 grid nodes using bilinear weights.
 *   3. NORMALIZE: divide every cell by the maximum value so the field
 *                 lies in [0, 1] before driving the mover.
 *   4. MOVER    : gather the normalised field back to each particle and
 *                 update its position.  Particles that leave [0,1]x[0,1]
 *                 are marked inactive and excluded from future iterations.
 *   5. DENORM   : multiply by max_val to restore original field magnitudes.
 *      (The denormalised mesh is what save_mesh() writes to disk.)
 */
void interpolation(double *mesh_value, Points *points) {

    int total = GRID_X * GRID_Y;

    /* ------------------------------------------------------------------ */
    /* 1. Reset mesh                                                        */
    /* ------------------------------------------------------------------ */
    memset(mesh_value, 0, total * sizeof(double));

    /* ------------------------------------------------------------------ */
    /* 2. SCATTER  point → mesh  (bilinear, fi = 1 for every particle)     */
    /* ------------------------------------------------------------------ */
    for (int p = 0; p < NUM_Points; p++) {
        if (!points[p].active) continue;

        double x = points[p].x;
        double y = points[p].y;

        /* Cell indices */
        int ix = (int)(x / dx);
        int iy = (int)(y / dy);
        if (ix >= NX) ix = NX - 1;   /* clamp boundary particles */
        if (iy >= NY) iy = NY - 1;

        /* Sub-cell offsets */
        double lx = x - ix * dx;
        double ly = y - iy * dy;

        /* Bilinear weights (area of the opposite rectangle / cell area) */
        double w00 = (dx - lx) * (dy - ly);   /* bottom-left  (ix  , iy  ) */
        double w10 =        lx * (dy - ly);   /* bottom-right (ix+1, iy  ) */
        double w01 = (dx - lx) *        ly;   /* top-left     (ix  , iy+1) */
        double w11 =        lx *        ly;   /* top-right    (ix+1, iy+1) */

        mesh_value[ iy      * GRID_X + ix    ] += w00;
        mesh_value[ iy      * GRID_X + ix + 1] += w10;
        mesh_value[(iy + 1) * GRID_X + ix    ] += w01;
        mesh_value[(iy + 1) * GRID_X + ix + 1] += w11;
    }

    /* ------------------------------------------------------------------ */
    /* 3. NORMALIZE mesh to [0, 1]                                         */
    /* ------------------------------------------------------------------ */
    double max_val = 0.0;
    for (int i = 0; i < total; i++)
        if (mesh_value[i] > max_val) max_val = mesh_value[i];

    double inv_max = (max_val > 1e-15) ? 1.0 / max_val : 0.0;
    for (int i = 0; i < total; i++)
        mesh_value[i] *= inv_max;

    /* ------------------------------------------------------------------ */
    /* 4. MOVER  mesh → point  (gather + position update)                 */
    /* ------------------------------------------------------------------ */
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

        /* Interpolated field value at particle location */
        double Fi = w00 * mesh_value[ iy      * GRID_X + ix    ]
                  + w10 * mesh_value[ iy      * GRID_X + ix + 1]
                  + w01 * mesh_value[(iy + 1) * GRID_X + ix    ]
                  + w11 * mesh_value[(iy + 1) * GRID_X + ix + 1];

        /* Update particle position */
        double new_x = x + Fi * dx;
        double new_y = y + Fi * dy;

        if (new_x < 0.0 || new_x > 1.0 || new_y < 0.0 || new_y > 1.0) {
            points[p].active = 0;   /* left domain → deactivate */
        } else {
            points[p].x = new_x;
            points[p].y = new_y;
        }
    }

    /* ------------------------------------------------------------------ */
    /* 5. DENORMALIZE  (restore original field magnitudes for output)      */
    /* ------------------------------------------------------------------ */
    for (int i = 0; i < total; i++)
        mesh_value[i] *= max_val;
}

/* ------------------------------------------------------------------ */
/* Write mesh to file                                                   */
/* ------------------------------------------------------------------ */
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
