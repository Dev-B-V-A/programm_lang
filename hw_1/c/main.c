// c-implementation gauss solver for hw_1 programm_languages

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int parse_args (int argc, char *argv[]);
void print (double *buf, int col, int row);
int init_matrix (double *matr, double *rhs, int size);
int gauss_solver (double *matr, double *x, double *rhs, int size);

int main (int argc, char *argv[])
{
    double *workspace = NULL;
    double *matr = NULL;
    double *x = NULL;
    double *rhs = NULL;
    int size = 0;

    size = parse_args(argc, argv);
    workspace = (double *)malloc (sizeof (double) * (size + 2) * size);
    if (workspace == NULL)
    {
        printf ("Cannot allocate memmory: %d", sizeof (double) * (size + 2) * size);
        return -1;
    }
    matr = workspace;
    rhs = matr + size * size;
    x = rhs + size;

    init_matrix (matr, rhs, size);

    gauss_solver (matr, x, rhs, size);

    printf ("Solve:\n");
    print (x, size, 1);
    printf (">>>>>>>>>>>>>>\n");

    if (workspace)
      free (workspace);
}

int parse_args (int argc, char *argv[])
{
    int read_size = 0;
    if (argc != 2)
    {
        printf ("programm: %s  required 1 argument - size matrix!\n", argv[0]);
        return -1;
    }
    if (sscanf (argv[1], "%d", &read_size) < 1)
        return -1;

    return read_size;
}

void print(double *buf, int col, int row)
{
    int i = 0, j = 0;
    if (buf == NULL)
        return;

    for (i = 0; i < row; i++)
    {
       for (j = 0; j < col; j++)
           printf ("%.4lf  ", buf[j + col * i]);
        printf ("\n");
    }
    printf ("\n");
}

int init_matrix(double *matr, double *rhs, int size)
{
    int i = 0, j = 0;
    double sum = 0.0;
    for (i = 0; i < size; ++i)
        for (j = 0; j < size; ++j)
            matr[j + i * size] = abs (i - j) + 1;
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
            sum += matr[j + i * size];
        rhs[i] = sum;
        sum = 0.0;
    }
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
        for (k = i + 1; k < size; ++k)
            matr[k + i * size] /= diag;

        for (j = i + 1; j < size; ++j)
        {
            for (k = i + 1; k < size; ++k)
                matr[k + j * size] -= matr[k + i * size] * matr[i + j * size];

            rhs[j] -= rhs[i] * matr[i + j * size];
        }
    }

    for (i = size - 1; i >=0; --i)
    {
        for (j = size - 1; j > i; --j)
        {
            sum += matr[j + i * size] * x[j];
        }
        x[i] = rhs[i] - sum;
        sum = 0.0;
    }

    return 0;
}
