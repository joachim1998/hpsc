#ifndef _IMPLICIT_H_
#define _IMPLICIT_H_

void implicitEuler(int myrank, int np, Parameters parameters, Map map, double* h,
                  double delta_xh, double delta_yh, MPI_Status status);

#endif