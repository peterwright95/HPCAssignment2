/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBI_H
#define _JACOBI_H

int jacobi(double ***, double ***, double ***, int, int, double);
int jacobi_baseline(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold);
int jacobi_collapse(double ***u, double ***u_upd, double ***f, int N, int iter_max, double threshold);
#endif
