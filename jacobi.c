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
    int i, j, k, num_iter = 0, threads = 0;
    double delta = 2 / (double)(N + 1);

    float distance = INFINITY;

    double s = 1.0 / 6.0;

    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        distance = 0;
        // #pragma omp parallel private(i, j, k) shared(distance, num_iter)
        #pragma omp parallel
        {
          // threads = omp_get_num_threads();
          #pragma omp for
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
        
        #pragma omp for
          for (i = 1; i < N + 1; i++)
              for (j = 1; j < N + 1; j++)
                  for (k = 1; k < N + 1; k++)
                  {
                      u[i][j][k] = u_upd[i][j][k];
                  }
      }
      num_iter += 1;
    }
    printf("%d %f %d %d ", N, distance, num_iter, threads);
    return threads;
}

int jacobi_reduction(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0, threads = 0;
    double delta = 2 / (double)(N + 1);

    float distance = INFINITY;

    double s = 1.0 / 6.0;

    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        distance = 0;
        // #pragma omp parallel private(i, j, k) shared(distance, num_iter)
        #pragma omp parallel
        {
          threads = omp_get_num_threads();
          #pragma omp for reduction(+: distance)
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
        
        #pragma omp for
          for (i = 1; i < N + 1; i++)
              for (j = 1; j < N + 1; j++)
                  for (k = 1; k < N + 1; k++)
                  {
                      u[i][j][k] = u_upd[i][j][k];
                  }
      }
      num_iter += 1;
    }
    printf("%d %f %d %d ", N, distance, num_iter, threads);
    return threads;
}

int jacobi_collapse(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = 0, threads;
    double delta = 2 / (double)(N + 1);

    float distance = INFINITY;

    // threads = omp_get_num_threads();

    double s = 1.0 / 6.0;

    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
      distance = 0;
      #pragma omp parallel
      {
      #pragma omp for reduction(+: distance)
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
      	// #pragma omp parallel for collapse(3) private(i, j, k) reduction(+: distance) schedule(static)
        #pragma omp for collapse(3) 
        for (i = 1; i < N + 1; i++)
        {
            for (j = 1; j < N + 1; j++)
            {
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = u_upd[i][j][k];
                }
            }
        }
      }
        num_iter += 1;
    }
    printf("%d %f %d %d ", N, distance, num_iter, threads);
    return threads;
}

int jacobi_ultimate(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold)
{
    int i, j, k, num_iter = -1, threads;
    double delta = 2 / (double)(N + 1);


    threads = omp_get_num_threads();

    float total_distance = INFINITY;
    double s = 1.0 / 6.0;

    float partial_distance = 0.0;

    
    #pragma omp parallel private(partial_distance, i , j ,k) shared(num_iter, total_distance)
    while (total_distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        #pragma omp single 
        {
        total_distance = 0.0;
        partial_distance = 0.0;
        num_iter++;
        }
        #pragma omp for
        for (i = 1; i < N + 1; i++)
        {
            for (j = 1; j < N + 1; j++)
            {
                for (k = 1; k < N + 1; k++)
                {
                    u_upd[i][j][k] = s*(u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1] + delta * delta * f[i][j][k]);

                    // Update the distance between u and u_upd
                    partial_distance += (u_upd[i][j][k] - u[i][j][k]) * (u_upd[i][j][k] - u[i][j][k]);
                }
            }
        }

        // Put u_upd into u to remake the calculations another time
        #pragma omp for
        for (i = 1; i < N + 1; i++)
        {
            for (j = 1; j < N + 1; j++)
            {
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = u_upd[i][j][k];
                }
            }
        }

        #pragma omp critical
        {
            total_distance += partial_distance;
        }
        #pragma omp barrier
    }
    printf("%d %f %d ", N, total_distance, num_iter);
    return threads;
}

