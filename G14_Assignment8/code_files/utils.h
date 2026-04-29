#ifndef UTILS_H
#define UTILS_H
#include "init.h"

// Serial / OpenMP functions
void interpolation(double *mesh_value, Points *points);
void normalize_mesh(double *mesh_value, double *out_min, double *out_max);
void denormalize_mesh(double *mesh_value, double vmin, double vmax);
void mover(double *mesh_value, Points *points);
void save_mesh(double *mesh_value);

// MPI + OpenMP functions
void interpolation_mpi(double *local_mesh, Points *local_pts, int count);
void mover_mpi(double *mesh_value, Points *local_pts, int count);

#endif
