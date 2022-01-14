/* gauss_seidel.c - Poisson problem in 3d
 
 */
#include <math.h>
#include <stdio.h>
#include <omp.h>

void gauss_seidel(double ***u, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0, count = 1;
    double delta = 2.0 / (double)(N + 1), u_upd;
    float distance = INFINITY;

    double s = 1.0 / 6.0;
    double u_before = 0.0;

    while (distance > threshold && num_iter < iter_max)
    {
        distance = 0;

        for (i = 1; i < N + 1; i++)
        {
            for (j = 1; j < N + 1; j++)
            {

                for (k = 1; k < N + 1; k++)
                {
                    u_upd = s * (u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1] + delta * delta * f[i][j][k]);

                    distance += (u_upd - u[i][j][k]) * (u_upd - u[i][j][k]);
                    u[i][j][k] = u_upd;
                }
            }
        }

        num_iter += 1;
    }
    printf("N: %d num_iter: %d ", N, num_iter);
}

void gauss_seidel_parallel(double ***u, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0, count = 1;
    double delta = 2.0 / (double)(N + 1);
    ;
    float distance = INFINITY;

    double s = 1.0 / 6.0;
    double u_before = 0.0;

    while (num_iter < iter_max)
    {
#pragma omp parallel shared(num_iter) private(j, k)
#pragma omp for ordered(2) schedule(static, 1)
        for (i = 1; i < N + 1; i++)
        {
            for (j = 1; j < N + 1; j++)
            {
#pragma omp ordered depend(sink                    \
                           : i - 1, j) depend(sink \
                                              : i, j - 1)
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = s * (u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1] + delta * delta * f[i][j][k]);
                }
#pragma omp ordered depend(source)
            }
        }

        num_iter += 1;
    }
    printf("%d,%d,", N, num_iter);
}
