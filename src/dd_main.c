/* C source code is found in dgemm_example.c */

#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "mkl.h"
#include "control.h"


void checkConsistency(Matrix *A, Matrix *B, Matrix *C, Matrix *D)
{
  // A must be square
  if (A->m != A->n)
    goto FAIL;
  
  // A and B must have the same number of rows
  if (A->m != B->m)
    goto FAIL;

  // A, C and D must have the same number of columns
  if (A->n != C->n || A->n != D->n)
    goto FAIL;
  
  return;

 FAIL:
  printf("Inconsistent matrix sizes\n");
  exit(1);
  return;
}

int main(int argc, char **argv)
{
  // MPI-related variables
  int rank, size;
  MPI_Request requests[4];
  MPI_Status status[4];
  double time1, time2;
  double *buf0, *buf1, *buf2, *buf3;
  int iter, max_iter;

  double beta, dt, h, betai;
  int Nt, Nti;
  Matrix *A, *B, *C, *D;
  Problem *problem;
  double *f, *g, *yhat, *g2;
  double a1, a2;
  double q, pp, qq;
  double *ygrid;
  int i, j, k, ndofi;
  double dtmp;   // temporary double
  double *tmpsol, *refsol;
  FILE *fid;
  Param *param;
  int length;
  char *fname;
  double *yhat_glob;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (i=0; i < 4; i++) {
    requests[i] = MPI_REQUEST_NULL;
  }

  //  max_iter = 140;

  // argv[1] contains local path
  length = strlen(argv[1]);
  fname = (char *)malloc((length+30)*sizeof(char));
  strcpy(fname, argv[1]);

  // Read parameters from file
  param = (Param *)malloc(sizeof(Param));
  strcat(fname, "/params.txt");
  readParams_par(fname, param, 8);

  /* Time horizon and time steps */
  beta = param->beta;
  Nt = param->Nt;
  dt = beta/Nt;
  h = dt;

  // Parameters for global problem
  a1 = param->a1;
  a2 = param->a2;
  q = param->q;
  // Best for four subdomains
  pp = param->pp;
  qq = param->qq;

  // Parameters for GMRES
  max_iter = param->max_iter;

  if (rank==0) {
    printf("=======================\n");
    printf("    Parameter Values   \n");
    printf("=======================\n");
    printf("# Subdomains = %d\n",size);
    printf("beta = %.5G\n",beta);
    printf("Nt = %d\n",Nt);
    printf("a1 = %.5G\n",a1);
    printf("a2 = %.5G\n",a2);
    printf("q = %.5G\n",q);
    printf("pp = %.5G\n",pp);
    printf("qq = %.5G\n",qq);
    printf("max_iter = %d\n",max_iter);
    printf("------------------------\n");
  }


  printf("Reading matrices...\n");
  /* Read in matrices */
  fname[length] = '\0';
  strcat(fname, "/A.txt");
  A = readMatrix_par(fname);
  fname[length] = '\0';
  strcat(fname, "/B.txt");
  B = readMatrix_par(fname);
  fname[length] = '\0';
  strcat(fname, "/C.txt");
  C = readMatrix_par(fname);
  fname[length] = '\0';
  strcat(fname, "/D.txt");
  D = readMatrix_par(fname);
  if (rank < size-1) {
    D = createIdentityMatrix(A->n);
  }

  printf("Checking consistency between sizes...\n");
  /* Check consistency between sizes */
  checkConsistency(A, B, C, D);

  // Local time horizons
  Nti = (rank+1)*Nt/size - rank*Nt/size;
  betai = dt * Nti;
  ndofi = A->n + B->n*Nti;
  tmpsol = (double *)malloc(ndofi * sizeof(double));
  refsol = (double *)malloc(ndofi * sizeof(double));


  // Read reference solution
  fname[length] = '\0';
  strcat(fname, "/sol.txt");
  readRefSol_par(fname, refsol, rank, size, Nt, A->n, B->n);


  printf("Allocating data...\n");
  /* Set up the rest of the problem */
  f  = NULL;   /* No inhomogeneous forcing term */
  g  = (double *)malloc(2* A->n * sizeof(double));
  g2 = g+A->n;
  //  g2 = (double *)malloc(D->m * sizeof(double));
  //  yhat = (double *)malloc(C->m * Nti * sizeof(double));
  buf0 = (double *)malloc(A->n * sizeof(double));
  buf1 = (double *)malloc(A->n * sizeof(double));

  printf("Initializing data... \n");
  // g = 0
  memset(g, 0, 2* A->n * sizeof(double));


  // Set up observation vector
  yhat_glob = (double *)malloc(C->m*Nt*sizeof(double));

  fname[length] = '\0';
  strcat(fname, "/yhat.txt");
  readDoubleVector_par(fname, yhat_glob, C->m, 0, Nt);

  yhat = &yhat_glob[C->m*rank*Nt/size];

  // g2 = last vector of yhat
  memcpy(g2, &(yhat_glob[(Nt-1)*C->m]), D->m);  


  /*

  // g2 = exp(-16*(y-0.5)^2)
  dtmp = h/2.0 - 0.5;
  for (i=0; i < D->m; i++) {
    g2[i] = exp(-16.0*dtmp*dtmp);
    dtmp += h;
  }

  // yhat = Nt copies of g2
  for (i=0; i < Nti; i++) {
    memcpy(&(yhat[i*(C->m)]), g2, C->m * sizeof(double));
  }

  */

  // g2 is an input trace except for the last subdomain, reset to zero
  if (rank < size-1) {
    memset(g2, 0, A->n * sizeof(double));
  }

  // Parameters for local subproblem
  if (rank < size-1) {
    a2 = pp;   // optimized parameter
  }
  if (rank > 0) {
    q = qq;    // optimized parameter
  }

  // For technical reasons, we set up g2 to be premultiplied by a2
  for (i=0; i < D->m; i++) {
    g2[i] *= a2;
  }

  time2 = MPI_Wtime();
  printf("Creating problem...\n");
  problem = createProblem(Nti, betai, A, B, C, D,
			  f, g, yhat, g2, a1, a2, q);


  printf("Starting DD iterations...\n");
  for (iter = 0; iter < max_iter; iter++) {
    printf("Rank %d, iteration %d:\n",rank, iter);
    /* Make sure initial and final conditions are ready */
    if (iter > 0) {
      MPI_Waitall(2, &requests[0], &status[0]);

      /* Switch buffers */
      if (rank > 0) {
	g = problem->g;
	problem->g = buf0;
	buf0 = g;
      }
      if (rank < size-1){
	g2 = problem->g2;
	problem->g2 = buf1;
	buf1 = g2;
      }
    }
    
    /* Ready to get new initial and final data */
    if (iter < max_iter-1) {
      if (rank > 0) {
	MPI_Irecv(buf0, A->n, MPI_DOUBLE, rank-1, iter, MPI_COMM_WORLD, &requests[0]);
      }
      if (rank < size-1) {
	MPI_Irecv(buf1, A->n, MPI_DOUBLE, rank+1, iter, MPI_COMM_WORLD, &requests[1]);
      }
    }

    /* Do local solve */
    solveProblem(problem);
    
    // Check solution
    memcpy(&(tmpsol[A->n]), &(refsol[A->n]), Nti*(B->n)*sizeof(double));
    cblas_daxpy(Nti*(B->n), -1.0, &(problem->sol[A->n]), 1, &(tmpsol[A->n]), 1);
    i = cblas_idamax(Nti*(B->n), &(tmpsol[A->n]), 1);
    printf("#%d %d %10.5G\n",iter,rank,tmpsol[A->n+i]);


    /* Send data to neighbours */
    if (iter < max_iter-1) {
      /* Calculate outbound traces */
      residual(tmpsol, problem->sol, problem);
      computeOutboundTraces(rank==0, rank==size-1, pp, qq, problem);

      MPI_Waitall(2, &requests[2], &status[2]);
      if (rank > 0) {
	MPI_Isend(problem->g, A->n, MPI_DOUBLE, rank-1, iter, MPI_COMM_WORLD, &requests[2]);
      }
      if (rank < size-1) {
	MPI_Isend(problem->g2, A->n, MPI_DOUBLE, rank+1, iter, MPI_COMM_WORLD, &requests[3]);
      }
    }
  }
  time1 = MPI_Wtime();
  printf("Iteration time for rank %d = %10.5G\n",rank, time1-time2);
  return (0);
}
