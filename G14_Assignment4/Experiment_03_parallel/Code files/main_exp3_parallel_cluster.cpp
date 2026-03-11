#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main()
{
    NX=1000;
    NY=400;

    NUM_Points=14000000;
    Maxiter=10;

    GRID_X=NX+1;
    GRID_Y=NY+1;

    dx=1.0/NX;
    dy=1.0/NY;

    omp_set_num_threads(4);

    double *mesh_value=
    (double*)calloc(GRID_X*GRID_Y,sizeof(double));

    Points *points=
    (Points*)calloc(NUM_Points,sizeof(Points));

    initializepoints(points);

    printf("Iter,Interp,Mover,Total\n");

    for(int iter=0;iter<Maxiter;iter++)
    {
        clock_t s1=clock();
        interpolation(mesh_value,points);
        clock_t e1=clock();

        clock_t s2=clock();
        mover_parallel(points,dx,dy);
        clock_t e2=clock();

        double interp=(double)(e1-s1)/CLOCKS_PER_SEC;
        double mover=(double)(e2-s2)/CLOCKS_PER_SEC;

        printf("%d,%lf,%lf,%lf\n",
        iter+1,interp,mover,interp+mover);
    }

    free(mesh_value);
    free(points);

    return 0;
}