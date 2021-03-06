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
  double beta, dt, h;
  int Nt;
  Matrix *A, *B, *C, *D;
  Problem *problem;
  double *f, *g, *yhat, *g2;
  double a1, a2;
  double q;
  double *ygrid;
  int i, ndof;
  double dtmp;   // temporary double
  FILE *fid;
  Param *param;
  char *fname;
  int length;
  double *r;

  // argv[1] contains local path, argv[2] contains prefix
  if(argc < 3)
  {
    printf("Please specify directory in which to run, and prefix. Aborting...\n");
    exit(3);
  }
  length = strlen(argv[1]) + 1 + strlen(argv[2]); // path/prefix
  fname = (char *)malloc((length+30)*sizeof(char));
  strcpy(fname, argv[1]);
  strcat(fname, "/");
  strcat(fname, argv[2]);

  // Read parameters from file
  param = (Param *)malloc(sizeof(Param));
  strcat(fname, "s2_params.txt");
  readParams(fname, param, 5);

  /* Time horizon and time steps */
  beta = param->beta;
  Nt = param->Nt;
  dt = beta/Nt;
  h = dt;

  // Parameters for global problem
  a1 = param->a1;
  a2 = param->a2;
  q = param->q;

  printf("Reading matrices...\n");
  /* Read in matrices */
  fname[length] = '\0';
  strcat(fname, "s2_A.txt");
  A = readMatrix(fname);
  fname[length] = '\0';
  strcat(fname, "s2_B.txt");
  B = readMatrix(fname);
  fname[length] = '\0';
  strcat(fname, "s2_C.txt");
  C = readMatrix(fname);
  fname[length] = '\0';
  strcat(fname, "s2_D.txt");
  D = readMatrix(fname);

  printf("Checking consistency between sizes...\n");
  /* Check consistency between sizes */
  checkConsistency(A, B, C, D);

  printf("Allocating data...\n");
  /* Set up the rest of the problem */
  f  = NULL;   /* No inhomogeneous forcing term */
  g  = (double *)malloc(A->n * sizeof(double));
  g2 = (double *)malloc(D->m * sizeof(double));
  yhat = (double *)malloc(C->m * Nt * sizeof(double));

  printf("Initializing data... \n");
  // g = 0
  memset(g, 0, A->n * sizeof(double));

  /*
  // g2 = exp(-16*(y-0.5)^2)
  dtmp = h/2.0 - 0.5;
  for (i=0; i < D->m; i++) {
    g2[i] = exp(-16.0*dtmp*dtmp);
    dtmp += h;
  }
  */
  printf("reading yhat...\n");
  fname[length] = '\0';
  strcat(fname, "s2_yhat.txt");
  readDoubleVector(fname, yhat, C->m*Nt);

  /*
  // yhat = Nt copies of g2
  for (i=0; i < Nt; i++) {
    memcpy(&(yhat[i*(C->m)]), g2, C->m * sizeof(double));
  }
  */

  // g2 = last vector of yhat
  memcpy(g2, &(yhat[(Nt-1)*C->m]), D->m);

  for (i=0; i < D->m; i++) {
    g2[i] *= a2;
  }

  printf("Creating problem...\n");
  problem = createProblem(Nt, beta, A, B, C, D,
			  f, g, yhat, g2, a1, a2, q);

  printf("Solving problem...\n");
  solveProblem(problem);
  ndof = A->n + B->n*Nt;
  r = (double *)malloc(ndof*sizeof(double));
  residual(r, problem->sol, problem);

  printf("Writing solution...\n");
  // Write vector
  fname[length] = '\0';
  strcat(fname, "s3_sol.txt");
  fid = fopen(fname,"w");
  for (i=0; i < ndof; i++) {
    fprintf(fid, "%.16e\n", problem->sol[i]);
  }
  fclose(fid);

  printf("Writing trajectory...\n");
  fname[length] = '\0';
  strcat(fname, "s3_yy.txt");
  fid = fopen(fname,"w");
  for (i=0; i < A->n*Nt; i++) {
    fprintf(fid, "%.16e\n", problem->yy[i]);
  }
  fclose(fid);

  printf("Done.\n");
  return 0;
}
