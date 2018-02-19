#ifndef GAUSS_H
#define GAUSS_H

int init_matrix (double *matr, double *rhs, int size);

int gauss_solver (double *matr, double *x, double *rhs, int size);

#endif // GAUSS_H
