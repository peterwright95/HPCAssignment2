/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

    int 	N = N_DEFAULT, M;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL, ***f = NULL, ***u_upd = NULL;
    int i, j, k;


    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // output type
    }
    M = N + 2;
    // allocate memory for u, u_upd and f
    if ( (u = d_malloc_3d(M, M, M)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( (u_upd = d_malloc_3d(M, M, M)) == NULL ) {
        perror("array u_upd: allocation failed");
        exit(-1);
    }
    if ( (f = d_malloc_3d(M, M, M)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
    }
    

    //Set initial values for u
    for(i = 0; i < M; i++){
        for(j = 0; j < M; j++){
            for(k = 0; k < M; k++){
                if ((i == 0) || (i == M) || (j == M) || (k == 0) || (k == M)){
                    u[i][j][k] = 20;
                }
                else if (j == 0){
                    u[i][j][k] = 0;
                }
                else{
                    u[i][j][k] = start_T;
                }

            }
        }
    }

    //Set values for f
    for(i = 0; i < M; i++){
        for(j = 0; j < M; j++){
            for(k = 0; k < M; k++){
                if ((i <= M * 3 / 16) && (j <= M / 4) && (k <= M / 2) && (k >= M * 2 / 16)){
                    f[i][j][k] = 200;
                }
                else{
                    f[i][j][k] = 0;
                }

            }
        }
    }

    //Obtain solution u using the Jacobi method (not sure if N or M, depends on implementation by Jasmine)
    #ifdef _JACOBI
    jacobi(u, u_upd, f, N, iter_max, tolerance);
    #endif

    #ifdef _GAUSS_SEIDEL
    gauss_seidel(u, u_upd, f, N, iter_max, tolerance);
    #endif

    // dump  results if wanted 
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free(u);

    return(0);
}
