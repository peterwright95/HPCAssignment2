/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdio.h>
#include <omp.h>

// N is the number of points per dimension
int jacobi(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0;
    double delta = 2 / (double)(N + 1);

    float distance = INFINITY;

    double s = 1.0 / 6.0;
    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        distance = 0;
        for (i = 1; i < N + 1; i++)
        {
            for (j = 1; j < N + 1; j++)
            {
                for (k = 1; k < N + 1; k++)
                {
                    u_upd[i][j][k] = s*(u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1] + delta * delta * f[i][j][k]);

                    // Update the distance between u and u_upd
                    distance += (u_upd[i][j][k] - u[i][j][k]) * (u_upd[i][j][k] - u[i][j][k]);
                }
            }
        }

        // Put u_upd into u to remake the calculations another time
        for (i = 1; i < N + 1; i++)
            for (j = 1; j < N + 1; j++)
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = u_upd[i][j][k];
                }

        num_iter += 1;
    }
    printf("%d %f %d ", N, distance, num_iter);
    return 0;
}


int jacobi_baseline(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0, threads;
    double delta = 2 / (double)(N + 1);

    float distance = INFINITY;

    threads = omp_get_num_threads();

    double s = 1.0 / 6.0;
    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        distance = 0;
        for (i = 1; i < N + 1; i++)
        {
            #pragma omp parallel for reduction(+: distance)
            for (j = 1; j < N + 1; j++)
            {
                for (k = 1; k < N + 1; k++)
                {
                    u_upd[i][j][k] = s*(u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1] + delta * delta * f[i][j][k]);

                    // Update the distance between u and u_upd
                    distance += (u_upd[i][j][k] - u[i][j][k]) * (u_upd[i][j][k] - u[i][j][k]);
                }
            }
        }

        // Put u_upd into u to remake the calculations another time
        for (i = 1; i < N + 1; i++)
            for (j = 1; j < N + 1; j++)
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = u_upd[i][j][k];
                }

        num_iter += 1;
    }
    printf("%d %f %d %d", N, distance, num_iter, threads);
    return threads;
}


int jacobi_collapse(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0, threads;
    double delta = 2 / (double)(N + 1);

    float distance = INFINITY;

    threads = omp_get_num_threads();

    double s = 1.0 / 6.0;
    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        distance = 0;
    		#pragma omp parallel shared(delta, u, u_upd, N, distance) private(i, j, k)
        for (i = 1; i < N + 1; i++)
        {
            #pragma omp for reduction(+: distance)
            for (j = 1; j < N + 1; j++)
            {
                for (k = 1; k < N + 1; k++)
                {
                    u_upd[i][j][k] = s*(u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1] + delta * delta * f[i][j][k]);

                    // Update the distance between u and u_upd
                    distance += (u_upd[i][j][k] - u[i][j][k]) * (u_upd[i][j][k] - u[i][j][k]);
                }
            }
        }

        // Put u_upd into u to remake the calculations another time
        for (i = 1; i < N + 1; i++)
            for (j = 1; j < N + 1; j++)
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = u_upd[i][j][k];
                }

        num_iter += 1;
    }
    printf("%d %f %d %d", N, distance, num_iter, threads);
    return threads;
}