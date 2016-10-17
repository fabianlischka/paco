#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "mkl.h"
#include "control.h"


void checkAlloc(void *ptr, char *proc_name, int size) {
  if (!ptr) {
    printf("%s: out of memory, size = %d\n", proc_name, size);
    exit(1);
  }
}

Matrix *readMatrix(char *filename) 
{
  Matrix *A;
  FILE *fid;
  int i, m, n, nnz;
  int *rowptr, *colind;
  double *val;

  fid = fopen(filename,"r");
  if (fid == NULL) {
    printf("File %s does not exist\n",filename);
    exit(1);
  }
  else {
    printf("Reading from %s\n",filename);
  }
  fscanf(fid, "%d",&m);
  fscanf(fid, "%d",&n);
  A = (Matrix *)malloc(sizeof(Matrix));
  A->m = m;
  A->n = n;
  rowptr = (int *)malloc((m+1)*sizeof(int));
  for (i=0; i <= m; i++) {
    fscanf(fid, "%d",&(rowptr[i]));
  }
  nnz = rowptr[m];
  colind = (int *)malloc(nnz*sizeof(int));
  val = (double *)malloc(nnz*sizeof(double));
  for (i=0; i < nnz; i++) {
    fscanf(fid, "%d",&(colind[i]));
  }
  for (i=0; i < nnz; i++) {
    fscanf(fid, "%le",&(val[i]));
  }
  A->rowptr = rowptr;
  A->colind = colind;
  A->val = val;
  fclose(fid);
  return (A);
}


Matrix *readMatrix_par(char *filename) 
{
  Matrix *A;
  MPI_File fid;
  int i, m, n, nnz;
  int *rowptr, *colind;
  double *val;
  int ierr;
  MPI_Offset filesize;
  char *chunk, *tok;

  ierr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);
  if (ierr) {
    printf("File %s does not exist\n",filename);
    exit(1);
  }
  MPI_File_get_size(fid, &filesize);
  chunk = (char *)malloc(filesize*sizeof(char));
  MPI_File_read_all(fid, chunk, filesize-1, MPI_CHAR, MPI_STATUS_IGNORE);
  chunk[filesize-1] = '\0';

  tok = strtok(chunk, " \n");
  sscanf(tok, "%d",&m);

  tok = strtok(NULL, " \n");
  sscanf(tok, "%d",&n);

  A = (Matrix *)malloc(sizeof(Matrix));
  A->m = m;
  A->n = n;
  rowptr = (int *)malloc((m+1)*sizeof(int));
  for (i=0; i <= m; i++) {
    tok = strtok(NULL, " \n");
    sscanf(tok, "%d",&(rowptr[i]));
  }
  nnz = rowptr[m];
  colind = (int *)malloc(nnz*sizeof(int));
  val = (double *)malloc(nnz*sizeof(double));
  for (i=0; i < nnz; i++) {
    tok = strtok(NULL, " \n");
    sscanf(tok, "%d",&(colind[i]));
  }
  for (i=0; i < nnz; i++) {
    tok = strtok(NULL, " \n");
    sscanf(tok, "%le",&(val[i]));
  }
  A->rowptr = rowptr;
  A->colind = colind;
  A->val = val;
  free(chunk);
  MPI_File_close(&fid);
  return (A);
}

void io_main(int argc, char **argv) {
  /* Unit test for reading in matrices */
  char *filename;
  Matrix *A;
  int i, j;
  
  filename = argv[1];
  A = readMatrix(filename);
  printf("Size of A is %d by %d\n",A->m, A->n);
  for (i=0; i < A->m; i++) {
    for (j=A->rowptr[i]; j < A->rowptr[i+1]; j++) {
      printf("(%2d,%2d)  %.16e\n",i, A->colind[j], A->val[j]);
    }
  }
}


void readDoubleVector(char *filename, double *dest, int n)
{
  FILE *fid;
  int i;

  fid = fopen(filename,"r");
  for (i=0; i < n; i++) {
    fscanf(fid, "%lg", &(dest[i]));
  }
  fclose(fid);
}

void readDoubleVector_par(char *filename, double *dest, int n, int start, int stop)
{
  MPI_File fid;
  int i, j, count;
  double tmp;
  int ierr;
  char *chunk, *tok;
  MPI_Offset filesize;

  ierr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);
  if (ierr) {
    printf("File %s does not exist\n",filename);
    exit(1);
  }
  else {
    printf("Reading from %s\n", filename);
  }
  MPI_File_get_size(fid, &filesize);
  chunk = (char *)malloc(filesize*sizeof(char));
  MPI_File_read_all(fid, chunk, filesize-1, MPI_CHAR, MPI_STATUS_IGNORE);
  chunk[filesize-1] = '\0';

  count = 0;
  
  tok = strtok(chunk, " \n");
  for (i=0; i < n*start; i++) {
    tok = strtok(NULL, " \n");
  }
  for (i=n*start; i < n*stop; i++){
    sscanf(tok, "%lg", &dest[count++]);
    tok = strtok(NULL, " \n");
  }
  free(chunk);
  MPI_File_close(&fid);
}

void readRefSol(char *filename, double *tmpsol, int rank, int size, int Nt, int An, int Bn)
{
  FILE *fid;
  int i, j, count;
  double tmp;
  int start, stop;

  fid = fopen(filename,"r");
  for (i=0; i < An; i++) {
    fscanf(fid,"%le", &tmp);
  }
  start = rank * Nt/size;
  for (i=0; i < start; i++) {
    fscanf(fid,"%le",&tmp);
  }
  stop = (rank+1) * Nt/size;
  count = An;
  for (i=start; i < stop; i++) {
    for (j=0; j < Bn; j++) {
      fscanf(fid,"%le",&(tmpsol[count++]));
    }
  }
  fclose(fid);
}

void readRefSol_par(char *filename, double *tmpsol, int rank, int size, int Nt, int An, int Bn)
{
  MPI_File fid;
  int i, j, count;
  double tmp;
  int start, stop;
  int ierr;
  char *chunk, *tok;
  //  char *filename = "/home/felix_kwok/tests/control/sol.txt";
  MPI_Offset filesize;

  ierr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);
  if (ierr) {
    printf("File %s does not exist\n",filename);
    exit(1);
  }
  MPI_File_get_size(fid, &filesize);
  chunk = (char *)malloc(filesize*sizeof(char));
  MPI_File_read_all(fid, chunk, filesize-1, MPI_CHAR, MPI_STATUS_IGNORE);
  chunk[filesize-1] = '\0';

  tok = strtok(chunk, " \n");
  for (i=1; i < An; i++) {
    tok = strtok(NULL, " \n");
    sscanf(tok,"%le", &tmp);
  }
  start = rank * Nt/size;
  for (i=0; i < start; i++) {
    for (j=0; j < Bn; j++) {
      tok = strtok(NULL, " \n");
      sscanf(tok,"%le",&tmp);
    }
  }
  stop = (rank+1) * Nt/size;
  count = An;
  for (i=start; i < stop; i++) {
    for (j=0; j < Bn; j++) {
      tok = strtok(NULL, " \n");
      sscanf(tok,"%le",&(tmpsol[count++]));
    }
  }
  free(chunk);
  MPI_File_close(&fid);
}

void readParams(char *filename, Param *param, int numParams)
{
  FILE *fid;
  int i, j, count;
  double tmp;
  int start, stop;
  int ierr;
  char *chunk, tok[257];
  MPI_Offset filesize;

  fid = fopen(filename, "r");
  if (fid == NULL) {
    printf("File %s does not exist\n", filename);
    exit(1);
  }
  else {
    printf("Reading from %s\n", filename);
  }
  count = 0;
  while (count < numParams) {
    fgets(tok, 256, fid);
    if (tok[0]=='/') {
      continue;
    }
    else {
      count++;
      switch(count) {
      case 1:
	sscanf(tok, "%lg", &(param->beta));
	break;
      case 2:
	sscanf(tok, "%d", &(param->Nt));
	break;
      case 3:
	sscanf(tok, "%lg", &(param->a1));
	break;
      case 4:
	sscanf(tok, "%lg", &(param->a2));
	break;
      case 5:
	sscanf(tok, "%lg", &(param->q));
	break;
      case 6:
	sscanf(tok, "%lg", &(param->pp));
	break;
      case 7:
	sscanf(tok, "%lg", &(param->qq));
	break;
      case 8:
	sscanf(tok, "%d", &(param->max_iter));
	break;
      case 9:
	sscanf(tok, "%lg", &(param->tol));
	break;
      case 10:
	sscanf(tok, "%d", &(param->krylov));
	break;
      default:
	printf("Error, should not be here\n");
	exit(1);
	break;
      }
    }
  }
  fclose(fid);
}



void readParams_par(char *filename, Param *param, int numParams)
{
  MPI_File fid;
  int i, j, count;
  double tmp;
  int start, stop;
  int ierr;
  char *chunk, *tok;
  MPI_Offset filesize;

  ierr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);
  if (ierr) {
    printf("File %s does not exist\n",filename);
    exit(1);
  }
  else {
    printf("Reading from %s\n", filename);
  }
  MPI_File_get_size(fid, &filesize);
  chunk = (char *)malloc(filesize*sizeof(char));
  MPI_File_read_all(fid, chunk, filesize-1, MPI_CHAR, MPI_STATUS_IGNORE);
  chunk[filesize-1] = '\0';

  count = 0;
  
  tok = strtok(chunk, "\n");
  while (count < numParams) {
    //    printf("%s",tok);
    if (tok[0]!='/') {
      //      printf("(Non-comment)\n");
      count++;
      switch(count) {
      case 1:
	sscanf(tok, "%lg", &(param->beta));
	//	printf("beta = %10.5G\n", param->beta);
	break;
      case 2:
	sscanf(tok, "%d", &(param->Nt));
	//	printf("Nt = %d\n", param->Nt);
	break;
      case 3:
	sscanf(tok, "%lg", &(param->a1));
	//	printf("a1 = %10.5G\n", param->a1);
	break;
      case 4:
	sscanf(tok, "%lg", &(param->a2));
	//	printf("a2 = %10.5G\n", param->a2);
	break;
      case 5:
	sscanf(tok, "%lg", &(param->q));
	//	printf("q = %10.5G\n", param->q);
	break;
      case 6:
	sscanf(tok, "%lg", &(param->pp));
	//	printf("pp = %10.5G\n", param->pp);
	break;
      case 7:
	sscanf(tok, "%lg", &(param->qq));
	//	printf("qq = %10.5G\n", param->qq);
	break;
      case 8:
	sscanf(tok, "%d", &(param->max_iter));
	//	printf("max_iter = %d\n", param->max_iter);
	break;
      case 9:
	sscanf(tok, "%lg", &(param->tol));
	//	printf("tol = %10.5G\n", param->tol);
	break;
      case 10:
	sscanf(tok, "%d", &(param->krylov));
	//	printf("krylov = %d\n", param->krylov);
	break;
      default:
	printf("Error, should not be here\n");
	exit(1);
	break;
      }
    }
    else {
      //      printf("(Comment)\n");
    }
    tok = strtok(NULL, "\n");
  }
  free(chunk);
  MPI_File_close(&fid);
}
