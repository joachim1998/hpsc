#ifndef _EXPLICIT_H_
#define _EXPLICIT_H_

void explicitEuler(int myrank, int np, Parameters parameters, Map map, double* h, double delta_xh,
                  double delta_yh, MPI_Status status);

#endif