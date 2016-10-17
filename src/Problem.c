#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "mkl.h"
#include "control.h"

/* --- Implementation starts here --- */

Problem *createProblem(double Nt, double beta, Matrix *A, Matrix *B,
		       Matrix *C, Matrix *D, double *f, double *g,
		       double *yhat, double *g2, double a1, double a2,
		       double q)
{
  Problem *problem = (Problem *)malloc(sizeof(Problem));
  problem->Nt = Nt;
  problem->beta = beta;
  problem->A = A;
  problem->B = B;
  problem->C = C;
  problem->D = D;
  problem->f = f;
  problem->g = g;
  problem->yhat = yhat;
  problem->g2 = g2;
  problem->a1 = a1;
  problem->a2 = a2;
  problem->q = q;
  problem->sol = NULL;
  problem->yy = NULL;
  problem->solver = NULL;
  return problem;
}



void freeProblem(Problem *problem)
{
  if (problem->sol) {
    free(problem->sol);
    problem->sol = NULL;
  }
  if (problem->yy) {
    free(problem->yy);
    problem->yy = NULL;
  }
  free(problem);
}

int solveProblem(Problem *problem)
{
  int ndof;                /* Number of degrees of freedom = size(u)+size(y0) */
  int k;                   /* max size of Krylov subspace */
  double *rhs, *rhs_copy;  /* Two copies of RHS (first one destroyed by GMRES */
  int ipar[128];           /* Workspaces for GMRES */
  double dpar[128];
  double *tmp;
  int tmp_size;
  int RCI_request;
  double *v_ptr, *w_ptr;
  int it_count;
  int counter;
  
  ndof = problem->A->n + problem->B->n * problem->Nt;
  k = 100;

  //  printf("Allocating vectors \n");
  /* Allocate vectors */
  if (!problem->sol) {
    problem->sol = (double *)malloc(ndof*sizeof(double));
    checkAlloc(problem->sol, "solveProblem", ndof*sizeof(double));
  }

  rhs = (double *)malloc(ndof*sizeof(double));
  rhs_copy = (double *)malloc(ndof*sizeof(double));
  tmp_size = (2*k+1)*ndof + k*(k+9)/2 + 1;
  tmp = (double *)malloc(tmp_size*sizeof(double));
  if (!(rhs && rhs_copy && tmp)) {
    printf("solveProblem: out of memory, ndof = %d, k = %d\n",ndof,k);
    exit(1);
  }
    
  // rhs = 0
  memset(rhs,0,ndof*sizeof(double));

  // Calculate RHS and store in rhs_copy
  residual(rhs_copy, rhs, problem);    /* rhs_copy = K*0 + r */

  // rhs = -r
  cblas_daxpy(ndof, -1.0, rhs_copy, 1, rhs, 1); 

  //  printf("Initializing GMRES\n");
  // Initialize GMRES solver to solve K*y = -r
  memset(problem->sol,0,ndof*sizeof(double));
  dfgmres_init(&ndof, problem->sol, rhs, &RCI_request, ipar, dpar, tmp);
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

  ipar[8]=1;
  ipar[9]=0;
  ipar[11]=1;
  ipar[14] = k;   // Also set Krylov subspace size to k
  dpar[0]=1.0E-8;
  
  /*---------------------------------------------------------------------------
  /* Check the correctness and consistency of the newly set parameters
  /*---------------------------------------------------------------------------*/
  
  dfgmres_check(&ndof, problem->sol, rhs, &RCI_request, ipar, dpar, tmp);
  if (RCI_request!=0) {
    printf("GMRES initialization failed\n");
    exit(1);
  }

  RCI_request = 1;
  counter = 0;
  while (RCI_request > 0) {
    dfgmres(&ndof, problem->sol, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request==1) {
      // printf("Residual at iteration %d: %10.4e\n",counter++,dpar[4]);
      // Multiply by K
      v_ptr = &(tmp[ipar[21]-1]);        /* Vector to be multiplied */
      w_ptr = &(tmp[ipar[22]-1]);        /* Result of mat-vec multiplication */
      residual(w_ptr, v_ptr, problem); /* w_ptr = K*v_ptr + r */
      cblas_daxpy(ndof, 1.0, rhs, 1, w_ptr, 1);  /* w_ptr -= r; */
    }
    else if (RCI_request !=0) {
      printf("GMRES failed, code = %d\n",RCI_request);
      exit(1);
    }
  }
  if (RCI_request==0) {
    // GMRES converged, get solution
    dfgmres_get(&ndof, problem->sol, rhs, &RCI_request, ipar, dpar, tmp, &it_count);        
    //    printf("GMRES converged in %d iterations\n",ipar[3]);
  }
  else {
    printf("GMRES failed, code = %d\n",RCI_request);
    exit(1);
  }

  // Clean up workspaces 
  free(rhs);
  free(rhs_copy);
  free(tmp);

  return(0);
}

int residual(double *r, double *x, Problem *problem)
{
  double *u;
  double dt;
  int yy_size, An, Bm, Bn, Cm, Cn, Dm, Dn;
  Matrix *M;
  int i, j;
  double *lam, *yi, *fi, *gtmp, *yy;
  int ndof, nRHS;
  int info;
  MKL_INT opt;
  double *lamsave;

  char notrans, trans;
  double c1, c2;
  char *matdescra;
  Solver *solver;

  //  printf("Initializing residual...\n");
  dt = (problem->beta)/(problem->Nt);
  An = problem->A->n;
  Bm = problem->B->m; Bn = problem->B->n;
  Cm = problem->C->m; Cn = problem->C->n;
  Dm = problem->D->m; Dn = problem->D->n;
  yy_size = An*(problem->Nt)+3*An+Dm;


  if (!problem->yy) {
    // Allocate space for solution and workspace
    problem->yy = (double *)malloc(yy_size*sizeof(double));
    checkAlloc(problem->yy, "residual", yy_size*sizeof(double));
  }

  yy = problem->yy;
  yi = yy+(An*(problem->Nt));
  lam  = yi + An;
  lamsave = lam + An;
  gtmp = lamsave + An;
  
  if (problem->solver==NULL) {
    //printf("Adjusting M...\n");
    M = IpcA(dt, problem->A);  
    //printf("Factorizing M...\n");
    solver = newSolver(M);
  }
  else {
    solver = problem->solver;
  }

  ndof = An + Bn * problem->Nt;

  //  printf("Set up time stepping...\n");
  memcpy(r, x, ndof * sizeof(double));
  memcpy(yi, x, An * sizeof(double));  // yi = y0;
  
  u = &(r[An]);
  fi = problem->f;
  for (i=0; i < problem->Nt; i++) {
    // yi += dt*B*u(:,i);
    notrans = 'N';
    c1 = dt;
    c2 = 1.0;
    matdescra = "GXXCXX";
    mkl_dcsrmv(&notrans, &Bm, &Bn, &c1, matdescra, problem->B->val, 
	       problem->B->colind, problem->B->rowptr, problem->B->rowptr+1, u, &c2, yi);
    u += Bn;

    // yi += dt*f(:,i);
    if (fi!=NULL) {
      cblas_daxpy(An, dt, fi, 1, yi, 1);
      fi += An;
    }

    // yy(:,i) = M\yi;
    opt = MKL_DSS_DEFAULTS;
    nRHS = 1;
    info = dss_solve_real(solver->handle, opt, yi, nRHS, &(yy[i*An]));
    if (info != MKL_DSS_SUCCESS) {
      printf("Solver returned error code %d\n", info);
      exit(1);
    }

    // yi = yy(:,i);
    memcpy(yi, &(yy[i*An]), An*sizeof(double));
  }

  // Careful: g2 has been redefined to contain a2 inside it
  // 
  // lln = D'*(g2-a2*D*yi) + a1*dt*C'(ynhat - C*yi);
  memcpy(gtmp, problem->g2, Dm*sizeof(double));
  notrans = 'N';
  c1 = -(problem->a2);
  c2 = 1.0;
  matdescra = "GXXCXX";
  mkl_dcsrmv(&notrans, &Dm, &Dn, &c1, matdescra,  problem->D->val,      // gtmp := g2 - a2*D*yi;
	      problem->D->colind,  problem->D->rowptr,  problem->D->rowptr+1, yi, &c2, gtmp);         
  memset(lam, 0, An * sizeof(double));                                  // lam := 0;
  trans = 'T';
  c1 = 1.0;
  c2 = 1.0;
  matdescra = "GXXCXX";
  mkl_dcsrmv(&trans, &Dm, &Dn, &c1, matdescra, problem->D->val,         // lam += D'*gtmp;
	     problem->D->colind, problem->D->rowptr, problem->D->rowptr+1, gtmp, &c2, lam);

  memcpy(gtmp, &(problem->yhat[(problem->Nt-1)*Cm]), Cm*sizeof(double));
  c1 = -1.0;
  c2 = 1.0;
  matdescra = "GXXCXX";
  mkl_dcsrmv(&notrans, &Cm, &Cn, &c1, matdescra,  problem->C->val,      // gtmp = ynhat - C*yi;
	      problem->C->colind,  problem->C->rowptr,  problem->C->rowptr+1, yi, &c2, gtmp);        

  c1 = (problem->a1)*dt;
  c2 = 1.0;
  matdescra = "GXXCXX";
  mkl_dcsrmv(&trans, &Cm, &Cn, &c1, matdescra, problem->C->val,         // lam += a1*dt*C'*gtmp;
	     problem->C->colind, problem->C->rowptr, problem->C->rowptr+1, gtmp, &c2, lam);

  opt = MKL_DSS_DEFAULTS + MKL_DSS_TRANSPOSE_SOLVE;
  nRHS = 1;
  info = dss_solve_real(solver->handle, opt, lam, nRHS, yi);         // yi := M'\lam;
  memcpy(lam, yi, An * sizeof(double));                              // lam := yi;

  // Make copy of lam for reconstruction of outbound traces if necessary
  memcpy(lamsave, yi, An*sizeof(double));
  
  u -= Bn;
  trans = 'T';
  c1 = -1.0;
  c2 = 1.0;
  matdescra = "GXXCXX";
  mkl_dcsrmv(&trans, &Bm, &Bn, &c1, matdescra, problem->B->val,       // u(:,end) -= B'*lam;
	     problem->B->colind, problem->B->rowptr, problem->B->rowptr+1, lam, &c2, u);        

  for (i = problem->Nt-1; i > 0; i--) {
    u -= Bn;

    memcpy(yi, &(problem->yhat[(i-1)*Cm]), Cm*sizeof(double));
    notrans = 'N';
    c1 = -1.0;
    c2 = 1.0;
    matdescra = "GXXCXX";
    mkl_dcsrmv(&notrans, &Cm, &Cn, &c1, matdescra, problem->C->val,    //yi := yhat(:,i)-C*yy(:,i)
	       problem->C->colind, problem->C->rowptr, problem->C->rowptr+1, &(yy[(i-1)*An]), &c2, yi);     

    trans = 'T';
    c1 = dt*problem->a1;
    c2 = 1.0;
    matdescra = "GXXCXX";
    mkl_dcsrmv(&trans, &Cm, &Cn, &c1, matdescra, problem->C->val,      //lam += a1*dt*C'*yi;
	       problem->C->colind, problem->C->rowptr, problem->C->rowptr+1, yi, &c2, lam);   

    opt = MKL_DSS_DEFAULTS + MKL_DSS_TRANSPOSE_SOLVE;
    nRHS = 1;
    info = dss_solve_real(solver->handle, opt, lam, nRHS, yi);         // yi := M'\lam;
    memcpy(lam, yi, An * sizeof(double));                              // lam := yi;

    trans = 'T';
    c1 = -1.0;
    c2 = 1.0;
    matdescra = "GXXCXX";
    mkl_dcsrmv(&trans, &Bm, &Bn, &c1, matdescra, problem->B->val,       // u(:,end) -= B'*lam;
	       problem->B->colind, problem->B->rowptr, problem->B->rowptr+1, lam, &c2, u);        
    
  }
  cblas_daxpy(An, -(problem->q), lam, 1, r, 1);                        // r := y0 - q*lam;
  cblas_daxpy(An, -1.0, problem->g, 1, r, 1);                            // r -= g;

  return (0);
}

Matrix *IpcA(double c, Matrix *A)
{
  int i, j, m, n, nnz;
  int diagFound;
  Matrix *M;

  M = (Matrix *)malloc(sizeof(Matrix));
  m = A->m;
  n = A->n;
  M->m = m;
  M->n = n;
  M->rowptr = (int *)malloc((m+1)*sizeof(int));
  memcpy(M->rowptr, A->rowptr, (m+1)*sizeof(int));
  nnz = M->rowptr[m];
  M->colind = (int *)malloc(nnz*sizeof(int));
  M->val = (double *)malloc(nnz*sizeof(double));
  memcpy(M->colind, A->colind, nnz*sizeof(int));

  for (i=0; i < n; i++) {
    diagFound = 0;
    for (j=A->rowptr[i]; j < A->rowptr[i+1]; j++) {
      M->val[j] = c*A->val[j];
      if (A->colind[j]==i) {
	M->val[j] += 1.0;
	diagFound = 1;
      }
    }
    if (!diagFound) {
      printf("Matrix A is missing a diagonal element.\n");
      exit(1);
    }
  }
  return (M);
}


Matrix *createIdentityMatrix(int size)
{
  Matrix *M;
  int i;

  M = (Matrix *)malloc(sizeof(Matrix));
  M->m = size;
  M->n = size;
  M->rowptr = (int *)malloc((size+1)*sizeof(int));
  M->colind = (int *)malloc(size*sizeof(int));
  M->val = (double *)malloc(size*sizeof(double));

  for (i=0; i < size; i++) {
    M->rowptr[i] = i;
    M->colind[i] = i;
    M->val[i] = 1.0;
  }
  M->rowptr[size] = size;
  return (M);
}

/* Stores outbound traces in problem->g and problem->g2. If isFirst!=0 or
 * isLast != 0, then the corresponding trace comes from the physical problem
 * and is not modified. 
 */
void computeOutboundTraces(int isFirst, int isLast, double pp, double qq, 
			   Problem *problem)
{
  double *y0, *yn;
  double *lam0, *lamn;
  double *y0hat, *ynhat;
  double *g, *g2;
  int An, Cm;
  double a1;
  double c1, c2, dt;
  Matrix *A, *C;
  char notrans, trans, *matdescra;
  int i;

  //  printf("Starting computeOutboundTraces\n");
  A = problem->A;
  C = problem->C;
  An = A->n;
  Cm = C->m;
  a1 = problem->a1;
  dt = (problem->beta)/(problem->Nt);

  //  printf("Locating vectors\n");
  y0 = problem->sol;   // First An entries of solution gives the initial condition
  yn = &(problem->yy[An*(problem->Nt-1)]);    // stores y at T_n
  lam0 = &(problem->yy[An*(problem->Nt+1)]);  // stores lambda at T_0
  lamn = &(problem->yy[An*(problem->Nt+2)]);  // stores lambda at T_n
  y0hat = problem->yhat;
  ynhat = &(problem->yhat[Cm*(problem->Nt-1)]);
  g = problem->g;
  g2 = problem->g2;

  /***********************************************************
   * If not the first subdomain, compute
   *      g = pp*y0 + lam0 
   ***********************************************************/

  //  i = cblas_idamax(An, y0, 1);
  //  printf("y0 for rank ?? = %10.5G\n",y0[i]);
  
  if (!isFirst) {
    /*
    // g = y0hat - C*y0
    memcpy(g, y0hat, Cm*sizeof(double));
    c1 = -1.0;
    c2 = 1.0;
    notrans = 'N';
    matdescra = "GXXCXX";
    mkl_dcsrmv(&notrans, &(C->m), &(C->n), &c1, matdescra, C->val,
	       C->colind, C->rowptr, C->rowptr+1, y0, &c2, g);

    // lam0 = lam0 + a1*dt*C'*g
    c1 = a1*dt;
    c2 = 1.0;
    trans = 'T';
    mkl_dcsrmv(&trans, &(C->m), &(C->n), &c1, matdescra, C->val,
	       C->colind, C->rowptr, C->rowptr+1, g, &c2, lam0);
    */
    // g = lam0 + pp*y0
    memcpy(g, lam0, An*sizeof(double));
    cblas_daxpy(An, pp, y0, 1, g, 1);
  }

  /***********************************************************
   * If not the last subdomain, compute
   *      g2 = yn - qq*lamn - qq*dt*A'*lamn +a1*qq*dt*C'*(ynhat - C*yn))
   ***********************************************************/
  if (!isLast) {
    // g2 = yn - qq*dt*A'*lamn
    memcpy(g2, yn, An*sizeof(double));

    c1 = -qq*dt;
    c2 = 1.0;
    trans = 'T';
    matdescra = "GXXCXX";
    mkl_dcsrmv(&trans, &(A->m), &(A->n), &c1, matdescra, A->val,
	       A->colind, A->rowptr, A->rowptr+1, lamn, &c2, g2);

    // g2 = g2 - qq*lamn
    cblas_daxpy(An, -qq, lamn, 1, g2, 1);

    // lamn = ynhat - C*yn
    memcpy(lamn, ynhat, Cm*sizeof(double));
    c1 = -1.0;
    c2 = 1.0;
    notrans = 'N';
    mkl_dcsrmv(&notrans, &(C->m), &(C->n), &c1, matdescra, C->val,
	       C->colind, C->rowptr, C->rowptr+1, yn, &c2, lamn);

    // g2 = g2 + a1*qq*dt*C'*lamn
    c1 = a1*qq*dt;
    c2 = 1.0;
    mkl_dcsrmv(&trans, &(C->m), &(C->n), &c1, matdescra, C->val,
	       C->colind, C->rowptr, C->rowptr+1, lamn, &c2, g2);

  }
}
