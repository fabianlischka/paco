#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "mkl.h"


int main(int argc, char **argv)
{
  // MPI-related variables
  int rank, size;

  MPI_File fid;
  int ierr;


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("Hello, this is node %d of %d.\n",rank,size);
  ierr = MPI_File_open(MPI_COMM_WORLD, "tests/control/A.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);
  if (ierr) {
    fprintf(stderr, "Rank %d: Couldn't open file A.txt\n", rank);
    MPI_Finalize();
    exit(2);
  }
  else {
    fprintf(stderr, "Rank %d: A.txt opened successfully\n", rank);
  }
  
  MPI_File_close(&fid);

  MPI_Finalize();
}
