#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

// Global definitions
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    int num_threads = omp_get_max_threads();
    printf("=== Parallel Interpolation (Mesh Privatization) ===\n");
    printf("Threads: %d\n", num_threads);

    FILE *file = fopen(argv[1], "rb");
    if (!file) {
        printf("Error: Cannot open input file: %s\n", argv[1]);
        exit(1);
    }

    // Read binary input parameters
    fread(&NX, sizeof(int), 1, file);
    fread(&NY, sizeof(int), 1, file);
    fread(&NUM_Points, sizeof(int), 1, file);
    fread(&Maxiter, sizeof(int), 1, file);

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    printf("Grid cells: %d x %d (points: %d x %d)\n", NX, NY, GRID_X, GRID_Y);
    printf("Particles: %d, Iterations: %d\n", NUM_Points, Maxiter);

    // Allocate memory
    double *mesh_value = (double*)calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points = (Points*)malloc(NUM_Points * sizeof(Points));

    double total_time = 0.0;

    // Main iteration loop
    for (int iter = 0; iter < Maxiter; iter++) {
        read_points(file, points);
        
        double start = omp_get_wtime();
        interpolation(mesh_value, points);
        double end = omp_get_wtime();
        
        total_time += (end - start);
    }

    save_mesh(mesh_value);
    
    // Performance metrics
    double avg_time = total_time / Maxiter;
    double throughput = (double)NUM_Points * Maxiter / total_time / 1e6;
    
    printf("\n=== Performance Results ===\n");
    printf("Total interpolation time: %.6f seconds\n", total_time);
    printf("Average per iteration:    %.6f seconds\n", avg_time);
    printf("Throughput:               %.2f million particles/second\n", throughput);
    printf("Output saved to: Mesh.out\n");

    // Cleanup
    free(mesh_value);
    free(points);
    fclose(file);

    return 0;
}