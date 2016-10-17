#ifndef INCL_CONTROL_H
#define INCL_CONTROL_H

typedef struct {
  int m, n;
  double *val;
  int *rowptr;
  int *colind;
} Matrix;

typedef struct {
  _MKL_DSS_HANDLE_t handle;
  Matrix *M;
} Solver;

typedef struct {
  int Nt;
  double beta;
  Matrix *A, *B, *C, *D;
  double *f, *g, *yhat, *g2;
  double a1, a2;
  double q;
  double *sol;
  double *yy;
  Solver *solver;
} Problem;

typedef struct {
  // Time horizon and time steps 
  double beta;
  int Nt;

  // Parameters for global problem
  double a1, a2;
  double q;

  // Optimized paramters for subproblems
  double pp, qq;
  int max_iter;
  double tol;

  // Parameters for GMRES
  int krylov;


} Param;

/* In Problem.c */
Problem *createProblem(double Nt, double beta, Matrix *A, Matrix *B,
		       Matrix *C, Matrix *D, double *f, double *g,
		       double *yhat, double *g2, double a1, double a2,
		       double q);

void freeProblem(Problem *problem);

int solveProblem(Problem *problem);

int residual(double *r, double *x, Problem *problem);

Matrix *IpcA(double c, Matrix *A);

Matrix *createIdentityMatrix(int size);

void computeOutboundTraces(int isFirst, int isLast, double pp, double qq, 
			   Problem *problem);

/* In Solver.c */
Solver *newSolver(Matrix *M);

void finalize_solver(Solver *solver);


/* In IO.c */

Matrix *readMatrix(char *filename);

Matrix *readMatrix_par(char *filename);

void checkAlloc(void *ptr, char *proc_name, int size);

void readRefSol(char *fname, double *tmpsol, int rank, int size, int Nt, int An, int Bn);

void readRefSol_par(char *fname, double *tmpsol, int rank, int size, int Nt, int An, int Bn);

void readParams(char *filename, Param *param, int numParams);

void readParams_par(char *filename, Param *param, int numParams);

void readDoubleVector(char *filename, double *dest, int n);

void readDoubleVector_par(char *filename, double *dest, int n, int start, int stop);

#endif
