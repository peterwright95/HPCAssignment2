/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void
gauss_seidel(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold) {
    int i, j, k, num_iter = 0;
    double delta = 2 / (double)(N + 1);;

    float distance = INFINITY;

    double s = 1.0 / 6.0;
    while (distance > threshold && num_iter < iter_max)
    {
        // Put the distance between u and u_upd to 0 now that the loop is started
        distance = 0;
        for (i = 1; i < N + 1; i++)
            for (j = 1; j < N + 1; j++)
                for (k = 1; k < N + 1; k++)
                {
                    u_upd[i][j][k] = u_upd[i - 1][j][k] + u_upd[i + 1][j][k] + u_upd[i][j - 1][k] + u_upd[i][j + 1][k] + u_upd[i][j][k - 1] + u_upd[i][j][k + 1] + delta * delta * f[i][j][k];
                    u_upd[i][j][k] *= s;

                    // Update the distance between u and u_upd
                    distance += (u_upd[i][j][k] - u[i][j][k]) * (u_upd[i][j][k] - u[i][j][k]);
                }

        // Put u_upd into u to remake the calculations another time
        for (i = 1; i < N + 1; i++)
            for (j = 1; j < N + 1; j++)
                for (k = 1; k < N + 1; k++)
                {
                    u[i][j][k] = u_upd[i][j][k];
                }
        
        if (num_iter < 10){
            printf("Iteration : %d\nDistance is: %f\n\n",num_iter,distance);
        }

        num_iter += 1;
    }
}

