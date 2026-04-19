#ifndef UTILS_H
#define UTILS_H
#include <time.h>
#include "init.h"

// Serial pipeline: scatter + normalize + mover + denormalize
void interpolation(double *mesh_value, Points *points);

// OpenMP parallel pipeline (defined in utils_omp.cpp)
void interpolation_omp(double *mesh_value, Points *points);

// Shared output
void save_mesh(double *mesh_value);

#endif
