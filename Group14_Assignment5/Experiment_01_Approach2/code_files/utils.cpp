#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "utils.h"

/* ── Bilinear (Cloud-In-Cell) Interpolation ────────────────── */
void interpolation(double *mesh_value, Points *points) {
    memset(mesh_value, 0, GRID_X * GRID_Y * sizeof(double));

    for (int i = 0; i < NUM_Points; i++) {
        if (!points[i].active) continue;

        // Find lower-left grid node
        int ix = (int)(points[i].x / dx);
        int iy = (int)(points[i].y / dy);

        // Clamp to valid range
        if (ix >= NX) ix = NX - 1;
        if (iy >= NY) iy = NY - 1;

        // Fractional distances
        double fx = (points[i].x - ix * dx) / dx;
        double fy = (points[i].y - iy * dy) / dy;

        // Bilinear weights to 4 surrounding nodes
        mesh_value[ iy      * GRID_X + ix    ] += (1-fx)*(1-fy);
        mesh_value[ iy      * GRID_X + ix + 1] += (  fx)*(1-fy);
        mesh_value[(iy + 1) * GRID_X + ix    ] += (1-fx)*(  fy);
        mesh_value[(iy + 1) * GRID_X + ix + 1] += (  fx)*(  fy);
    }
}

/* ── Serial Mover (Assignment 4 — particles stay in domain) ── */
void mover_serial(Points *points, double deltaX, double deltaY) {
    for (int i = 0; i < NUM_Points; i++) {
        double rx = ((double)rand()/RAND_MAX * 2.0 - 1.0) * deltaX;
        double ry = ((double)rand()/RAND_MAX * 2.0 - 1.0) * deltaY;
        points[i].x += rx;
        points[i].y += ry;
        // Reflect/clamp — keep particle inside
        if (points[i].x < 0.0) points[i].x = 0.0;
        if (points[i].x > 1.0) points[i].x = 1.0;
        if (points[i].y < 0.0) points[i].y = 0.0;
        if (points[i].y > 1.0) points[i].y = 1.0;
    }
}

/* ── Parallel Mover (Assignment 4) ───────────────────────────*/
void mover_parallel(Points *points, double deltaX, double deltaY) {
    #pragma omp parallel
    {
        unsigned int seed = (unsigned int)(omp_get_thread_num() * 123456 + time(NULL));
        #pragma omp for schedule(static)
        for (int i = 0; i < NUM_Points; i++) {
            double rx = ((double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0) * deltaX;
            double ry = ((double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0) * deltaY;
            points[i].x += rx;
            points[i].y += ry;
            if (points[i].x < 0.0) points[i].x = 0.0;
            if (points[i].x > 1.0) points[i].x = 1.0;
            if (points[i].y < 0.0) points[i].y = 0.0;
            if (points[i].y > 1.0) points[i].y = 1.0;
        }
    }
}

/* ════════════════════════════════════════════════════════════
 * ASSIGNMENT 5 — APPROACH 1: DEFERRED INSERTION
 *
 *  Step 1: Move all particles. Mark out-of-domain ones inactive.
 *          Record their indices in a "void list".
 *  Step 2: After the loop, insert new particles into those slots.
 * ════════════════════════════════════════════════════════════ */
void mover_deferred(Points *points, double deltaX, double deltaY, long *deleted) {

    int max_threads = omp_get_max_threads();

    // Per-thread seed array
    unsigned int *seeds = (unsigned int *)malloc(max_threads * sizeof(unsigned int));
    for (int t = 0; t < max_threads; t++)
        seeds[t] = (unsigned int)(t * 999983 + time(NULL));

    // Per-thread void lists (indices of deleted particles)
    long **local_voids    = (long **)malloc(max_threads * sizeof(long *));
    long  *local_void_cnt = (long  *)calloc(max_threads, sizeof(long));
    long chunk = NUM_Points / max_threads + 1024;
    for (int t = 0; t < max_threads; t++)
        local_voids[t] = (long *)malloc(chunk * sizeof(long));

    /* ── Step 1: Parallel move + collect void indices ── */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        unsigned int seed = seeds[tid];

        #pragma omp for schedule(static)
        for (int i = 0; i < NUM_Points; i++) {
            double rx = ((double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0) * deltaX;
            double ry = ((double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0) * deltaY;
            points[i].x += rx;
            points[i].y += ry;

            if (points[i].x < 0.0 || points[i].x > 1.0 ||
                points[i].y < 0.0 || points[i].y > 1.0) {
                points[i].active = 0;
                local_voids[tid][local_void_cnt[tid]++] = i;
            }
        }
        seeds[tid] = seed;
    }

    /* ── Merge void lists into a single flat array ── */
    long total_voids = 0;
    for (int t = 0; t < max_threads; t++)
        total_voids += local_void_cnt[t];

    long *void_idx = (long *)malloc(total_voids * sizeof(long));
    long offset = 0;
    for (int t = 0; t < max_threads; t++) {
        memcpy(void_idx + offset, local_voids[t], local_void_cnt[t] * sizeof(long));
        offset += local_void_cnt[t];
        free(local_voids[t]);
    }
    free(local_voids);
    free(local_void_cnt);

    /* ── Step 2: Parallel insertion into void slots ── */
    #pragma omp parallel for schedule(static)
    for (long k = 0; k < total_voids; k++) {
        int tid = omp_get_thread_num();
        long idx = void_idx[k];
        points[idx].x      = (double)rand_r(&seeds[tid]) / RAND_MAX;
        points[idx].y      = (double)rand_r(&seeds[tid]) / RAND_MAX;
        points[idx].active = 1;
    }

    *deleted = total_voids;

    free(void_idx);
    free(seeds);
}

/* ════════════════════════════════════════════════════════════
 * ASSIGNMENT 5 — APPROACH 2: IMMEDIATE REPLACEMENT
 *
 *  For each particle: move it.
 *  If it exits the domain → immediately replace it in-place
 *  with a new random particle. No separate insertion step.
 * ════════════════════════════════════════════════════════════ */
void mover_immediate(Points *points, double deltaX, double deltaY, long *deleted) {

    int max_threads = omp_get_max_threads();
    unsigned int *seeds = (unsigned int *)malloc(max_threads * sizeof(unsigned int));
    for (int t = 0; t < max_threads; t++)
        seeds[t] = (unsigned int)(t * 777777 + time(NULL));

    long total_deleted = 0;

    #pragma omp parallel reduction(+:total_deleted)
    {
        int tid = omp_get_thread_num();
        unsigned int seed = seeds[tid];

        #pragma omp for schedule(static)
        for (int i = 0; i < NUM_Points; i++) {
            double rx = ((double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0) * deltaX;
            double ry = ((double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0) * deltaY;
            points[i].x += rx;
            points[i].y += ry;

            if (points[i].x < 0.0 || points[i].x > 1.0 ||
                points[i].y < 0.0 || points[i].y > 1.0) {
                // Replace immediately at same memory slot
                points[i].x = (double)rand_r(&seed) / RAND_MAX;
                points[i].y = (double)rand_r(&seed) / RAND_MAX;
                total_deleted++;
            }
            points[i].active = 1;
        }
        seeds[tid] = seed;
    }

    *deleted = total_deleted;
    free(seeds);
}

/* ── Save mesh to file ───────────────────────────────────── */
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