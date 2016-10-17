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
  double *tmpsol;
  FILE *fid;
  double *initCond, *finCond;
  char *fname;
  int length;

  // GMRES-related variables
  int gndof;                /* Number of degrees of freedom = 2*An*size*/
  int krylov;               /* max size of Krylov subspace */
  double *rhs, *rhs_copy;  /* Two copies of RHS (first one destroyed by GMRES */
  int ipar[128];           /* Workspaces for GMRES */
  double dpar[128];
  double *tmp;
  int tmp_size;
  int RCI_request;
  double *v_ptr, *w_ptr;
  int offset, it_count;
  double tol;
  Param *param;
  char *aux_str;
  double *yhat_glob;
  int count;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (i=0; i < 4; i++) {
    requests[i] = MPI_REQUEST_NULL;
  }

  // argv[1] contains local path
  length = strlen(argv[1]);
  fname = (char *)malloc((length+30)*sizeof(char));
  strcpy(fname, argv[1]);

  //  printf("Hello, this is node %d of %d.\n",rank,size);

  // Read parameters from file
  param = (Param *)malloc(sizeof(Param));
  strcat(fname, "/params.txt");
  readParams_par(fname, param, 10);
  //  printf("Outside: beta = %10.5G\n",param->beta);
  //  printf("Outside: krylov = %d\n",param->krylov);
  //  printf("I was supposed to read from %s\n",fname);

  /* Time horizon and time steps */
  beta = param->beta;
  //  beta = 4.0;
  Nt = param->Nt;
  //  Nt = 64;
  dt = beta/Nt;
  h = dt;

  // Parameters for global problem
  a1 = param->a1;
  //  a1 = 100.0;
  a2 = param->a2;
  //  a2 = 0.0;
  q = param->q;
  //  q = 0.0;
  // Best for four subdomains
  pp = param->pp;
  //  pp = 6.5;
  qq = param->qq;
  //  qq = 0.065;
  //  pp = 0.1;  // Best for eight subdomains
  //  qq = 0.1;

  // Parameters for GMRES
  krylov = param->krylov;
  //  krylov = 200;
  max_iter = param->max_iter;
  //  max_iter = 200;
  tol = param->tol;
  //  tol = 1.0e-6;

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
    printf("tol = %.5G\n",tol);
    printf("max_iter = %d\n",max_iter);
    printf("krylov = %d\n",krylov);
    printf("------------------------\n");
  }

  //  printf("Reading matrices...\n");
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
  //  printf("Checking consistency between sizes...\n");
  /* Check consistency between sizes */
  checkConsistency(A, B, C, D);

  // Local time horizons
  Nti = (rank+1)*Nt/size - rank*Nt/size;
  betai = dt * Nti;
  ndofi = A->n + B->n*Nti;
  tmpsol = (double *)malloc(ndofi * sizeof(double));

  //  printf("Allocating data...\n");
  /* Set up the rest of the problem */
  f  = NULL;   /* No inhomogeneous forcing term */

  g  = (double *)malloc(2*(A->n) * sizeof(double));    // g and g2 need to be contiguous
  g2 = g+A->n;                                         // g2 only needs D->m entries, but
                                                       // to make scatter/gather easier, we
                                                       // allocate A->n space
  //  yhat = (double *)malloc(C->m * Nti * sizeof(double));

  //  printf("Initializing data... \n");
  // g = 0
  memset(g, 0, 2*(A->n) * sizeof(double));
  /*
  // g2 = exp(-16*(y-0.5)^2)
  dtmp = h/2.0 - 0.5;
  for (i=0; i < D->m; i++) {
    g2[i] = exp(-16.0*dtmp*dtmp);
    dtmp += h;
  }

  // yhat = Nt copies of g2
  // Potential bug! yhat need not have the same length as g2
  for (i=0; i < Nti; i++) {
    memcpy(&(yhat[i*(C->m)]), g2, C->m * sizeof(double));
  }
  */

  yhat_glob = (double *)malloc(C->m*Nt*sizeof(double));

  fname[length] = '\0';
  strcat(fname, "/yhat.txt");
  readDoubleVector_par(fname, yhat_glob, C->m, 0, Nt);

  yhat = &yhat_glob[C->m*rank*Nt/size];

  // g2 = last vector of yhat
  memcpy(g2, &(yhat_glob[(Nt-1)*C->m]), D->m);  

  if (rank==0) {
    // Head node only
    // Save initial and final conditions, since they get destroyed all the time
    initCond = (double *)malloc(A->n * sizeof(double));
    memcpy(initCond, g, A->n * sizeof(double));

    finCond = (double *)malloc(A->n * sizeof(double));
    memset(finCond, 0, A->n * sizeof(double));
    cblas_daxpy(D->m, a2, g2, 1, finCond, 1);
  }

  // Parameters for local subproblem
  if (rank < size-1) {
    a2 = pp;   // optimized parameter
  }
  if (rank > 0) {
    q = qq;    // optimized parameter
  }

  //  printf("Creating problem...\n");
  problem = createProblem(Nti, betai, A, B, C, D,
			  f, g, yhat, g2, a1, a2, q);

  // For technical reasons, we set up g2 to be premultiplied by a2
  for (i=0; i < D->m; i++) {
    g2[i] *= a2;
  }

  printf("Starting DD iterations...\n");
  // GMRES version: Head node needs extra storage for gmres
  gndof = 2*(A->n)*size;
  if (krylov > gndof) {
    krylov = gndof;
  }
  if (rank==0) {
    // Head node only    
    // For GMRES
    time2 = MPI_Wtime();
    rhs = (double *)malloc(gndof*sizeof(double));
    rhs_copy = (double *)malloc(gndof*sizeof(double));
    tmp_size = (2*krylov+1)*gndof + krylov*(krylov+9)/2 + 1;
    tmp = (double *)malloc(tmp_size*sizeof(double));
    if (!(rhs && rhs_copy && tmp)) {
      printf("solveProblem: out of memory, ndof = %d, k = %d\n",gndof,krylov);
      exit(1);
    }

    // rhs = 0
    memset(rhs,0,gndof*sizeof(double));
  }

  // Send input traces to subdomains
  //  MPI_Scatter(rhs, 2*(A->n), MPI_DOUBLE, g, 2*(A->n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  memset(g, 0, 2*(A->n)*sizeof(double));

  // Run simulation and calculate outbound traces
  solveProblem(problem);
  residual(tmpsol, problem->sol, problem);
  computeOutboundTraces(rank==0, rank==size-1, pp, qq, problem);
  if (rank==size-1) {
    memset(g2+D->m,0,(A->n - D->m)*sizeof(double));
  }

  // Gather these new traces
  MPI_Gather(g, 2*(A->n), MPI_DOUBLE, rhs_copy, 2*(A->n), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Initialize GMRES
  if (rank == 0) {
    fname[length] = '\0';
    strcat(fname, "/err_gmres.txt");
    fid = fopen(fname,"w");
    // Set correct rhs with inhomogeneous g and g2
    // g = 0
    memcpy(rhs_copy, initCond, A->n*sizeof(double));
    // Set g2 to be the same as yhat
    memcpy(rhs_copy+(A->n)*(2*size-1), finCond, A->n*sizeof(double));

    // rhs = -r
    cblas_daxpy(gndof, -1.0, rhs_copy, 1, rhs, 1); 
    
    //memset(rhs_copy, 0, gndof*sizeof(double));
    dfgmres_init(&gndof, rhs_copy, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request!=0) {
      printf("GMRES initialization failed\n");
      exit(1);
    }

    /*---------------------------------------------------------------------------
      /* Set the desired parameters:
      /* LOGICAL parameters:
      /* do residual stopping test
      /* do not request for the user defined stopping test
      /* do the check of the norm of the next generated vector automatically
      /* DOUBLE PRECISION parameters
      /* set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
      /*---------------------------------------------------------------------------*/

    ipar[4]=max_iter;         // maximum number of iterations
    ipar[8]=1;           // perform residual stopping test
    ipar[9]=0;
    ipar[11]=1;          // No additional preconditioning
    ipar[14] = krylov;   // Also set Krylov subspace size to krylov
    dpar[0]=tol;
  
    /*---------------------------------------------------------------------------
      /* Check the correctness and consistency of the newly set parameters
      /*---------------------------------------------------------------------------*/
  
    dfgmres_check(&gndof, rhs_copy, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request!=0) {
      printf("GMRES initialization failed\n");
      exit(1);
    }
  }
  RCI_request = 1;
  while (RCI_request > 0) {
    if (rank == 0) {
      dfgmres(&gndof, rhs_copy, rhs, &RCI_request, ipar, dpar, tmp);
      //      printf("Global iteration %d: residual = %10.5G\n", ipar[3], dpar[4]);
      //      fprintf(fid,"err(%d) = %10.5G;\n",ipar[3]+1,dpar[4]);
    }
    MPI_Bcast(&RCI_request, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (RCI_request==1) {
	// Multiply by K
      if (rank == 0) {
	v_ptr = &(tmp[ipar[21]-1]);        /* Vector to be multiplied */
	w_ptr = &(tmp[ipar[22]-1]);        /* Result of mat-vec multiplication */
      }
      else {
	v_ptr = NULL; w_ptr = NULL;
      }

      MPI_Scatter(v_ptr, 2*A->n, MPI_DOUBLE, g, 2*A->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      // Run simulation and calculate outbound traces
      solveProblem(problem);
      residual(tmpsol, problem->sol, problem);
      computeOutboundTraces(rank==0, rank==size-1, pp, qq, problem);
      
      MPI_Gather(g, 2*A->n, MPI_DOUBLE, w_ptr, 2*A->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (rank == 0) {
	// Reset initial and final condition
	memcpy(w_ptr, initCond, A->n*sizeof(double));
	memcpy(w_ptr+A->n*(2*size-1), finCond, A->n*sizeof(double));
	
	// Subtract to get jump of traces
	offset = 0;
	cblas_daxpy(A->n, -1.0, v_ptr, 1, w_ptr, 1);
	offset += A->n;
	for (i=0; i < size-1; i++) {
	  cblas_daxpy(A->n, -1.0, v_ptr+offset, 1, w_ptr+offset+A->n, 1);
	  cblas_daxpy(A->n, -1.0, v_ptr+offset+A->n, 1, w_ptr+offset, 1);
	  offset += 2*A->n;
	}
	cblas_daxpy(A->n, -1.0, v_ptr+offset, 1, w_ptr+offset, 1);
      
	// Now subtract RHS to get homogeneous problem
	cblas_daxpy(gndof, 1.0, rhs, 1, w_ptr, 1);  /* w_ptr -= r; */
      }
    }
    else if (RCI_request !=0) {
      printf("GMRES failed, code = %d\n",RCI_request);
      exit(1);
    }
  }
  if (RCI_request==0) {
    // GMRES converged, get solution
    if (rank == 0) {
      dfgmres_get(&gndof, rhs_copy, rhs, &RCI_request, ipar, dpar, tmp, &it_count);    

      // toc
      time1 = MPI_Wtime();
      printf("Global GMRES converged in %d iterations\n",it_count);
      printf("Wall clock time = %10.5G\n",rank, time1-time2);
    }
    MPI_Scatter(rhs_copy, 2*A->n, MPI_DOUBLE, g, 2*A->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    solveProblem(problem);
    residual(tmpsol, problem->sol, problem);


  }
  else {
    printf("GMRES failed, code = %d\n",RCI_request);
    exit(1);
  }

  // Check solution
  fname[length] = '\0';
  strcat(fname, "/sol.txt");
  readRefSol_par(fname, tmpsol, rank, size, Nt, A->n, B->n);
  cblas_daxpy(Nti*(B->n), -1.0, &(problem->sol[A->n]), 1, &(tmpsol[A->n]), 1);
  i = cblas_idamax(Nti*(B->n), &(tmpsol[A->n]), 1);
  printf("Error for rank %d = %10.5G\n",rank,tmpsol[A->n+i]);

  if (rank==0) {
    fclose(fid);
  }
  /*  
  // Write solution
  count = 0;
  for (i=0; i < A->n*Nti; i++) {
    printf("out: %d %d %20.16g\n",rank,count++,problem->yy[i]);
  }
  */
  
  // Clean up workspaces 
  if (rank == 0) {
    free(rhs);
    free(rhs_copy);
    free(tmp);
  }

  return (0);
}


