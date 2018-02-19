// c-implementation gauss solver for hw_1 programm_languages


#include "gauss.h"
#include "stdio.h"
#include "stdlib.h"

int parse_args (int argc, char *argv[]);
void print (double *buf, int col, int row = 1);

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

    init_matrix (matr, x, rhs, size);

    print (matr, size, size);
    print (rhs, size);

    gauss_solver (matr, x, rhs, size);

    print (x, size);

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

int print(double *buf, int col, int row)
{
    int i = 0, j = 0;
    if (buf == NULL)
        return;

    for (i = 0; i < row; i++)
    {
       for (j = 0; j < col; j++)
           printf ("%.4lf", buf[j + col * i]);
        printf ("\n");
    }
}
