#include "gauss.h"
#include "math.h"

int init_matrix(double *matr, double *rhs, int size)
{
    int i = 0, j = 0;

    for (i = 0; i < size; ++i)
        for (j = 0; j < size; ++j)
            matr[j + i * size] = abs (i + j + 1);
    for (i = 0; i < size; ++i)
        rhs[i] = i + j + 1;
}

int gauss_solver(double *matr, double *x, double *rhs, int size)
{
    int i = 0, j = 0, k = 0;
    double diag = 0.0;
    double sum = 0.0;

    for (i = 0; i < size; ++i)
    {
        diag = matr[i + i * size];
        rhs [i] /= diag;
        for (j = i + 1; j < size; ++j)
        {
            matr[j + i * size] /= diag;
            for (k = i + 1; k < size; ++k)
            {
                matr[k + j * size] -= matr[k + i * size] * matr[i + j * size];
            }
            rhs[j] -= rhs[j] * matr[i + j * size];
        }
    }

    for (i = size - 1; i >=0; --i)
    {
        for (j = size - 1; j > i; --j)
        {
            sum += matr[j + i * size] * x[j];
        }
        x[i] = x[i] - sum;
        sum = 0.0;
    }
    return 0;
}
