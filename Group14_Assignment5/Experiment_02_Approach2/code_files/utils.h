#ifndef UTILS_H
#define UTILS_H
#include <time.h>
#include "init.h"

void interpolation(double *mesh_value, Points *points);
void mover_serial(Points *points, double deltaX, double deltaY);
void mover_parallel(Points *points, double deltaX, double deltaY);

// NEW for Assignment 5
void mover_deferred(Points *points, double deltaX, double deltaY, long *deleted);
void mover_immediate(Points *points, double deltaX, double deltaY, long *deleted);

void save_mesh(double *mesh_value);
#endif