#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "mkl.h"
#include "control.h"

void factorize(Solver *solver);

Solver *newSolver(Matrix *M) {
  Solver *solver;
  solver = (Solver *)malloc(sizeof(Solver));
  solver->M = M;
  factorize(solver);
  return (solver);
}


void factorize(Solver *solver) {
  MKL_INT opt = MKL_DSS_DEFAULTS;
  MKL_INT zero_based = opt + MKL_DSS_ZERO_BASED_INDEXING;
  MKL_INT sym = MKL_DSS_NON_SYMMETRIC;
  int info;
  int n;

  n = solver->M->n;

  // printf(" DSS_CREATE \n\n ");
  info = dss_create(solver->handle, zero_based);
  if (info != MKL_DSS_SUCCESS) {
	printf("Solver returned error code %d\n", info);
	exit(1);
  }

  // printf(" Structure definition \n\n ");
  info = dss_define_structure(solver->handle, sym, solver->M->rowptr, n, n, solver->M->colind, solver->M->rowptr[n]);
  if (info != MKL_DSS_SUCCESS) {
    printf("Solver returned error code %d\n", info);
    exit(1);
  }
  
  // printf(" Reordering \n\n ");
  info = dss_reorder(solver->handle, opt, 0);
  if ( info != MKL_DSS_SUCCESS ) {
    printf("Solver returned error code %d\n", info);
    exit(1);
  }

  // printf(" Factorization \n\n ");
  info = dss_factor_real(solver->handle, opt, solver->M->val);
  if (info != MKL_DSS_SUCCESS) {
    printf("Solver returned error code %d\n", info);
    exit(1);
  }

}

void finalize_solver(Solver *solver)
{
  MKL_INT opt = MKL_DSS_DEFAULTS;
  dss_delete(solver->handle, opt);
}
