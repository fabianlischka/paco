CC = mpicc
MPI_COMPILE_FLAGS = -O2
MPI_LINK_FLAGS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIBRARY_PATH = /u1/local/intel/mkl/lib/intel64
INCLUDE_PATH = /u1/local/intel/mkl/include
DEPS = control.h
OBJ = Problem.o Solver.o IO.o

%.o: %.c $(DEPS)
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_COMPILE_FLAGS) -c -o $@ $<

IO: IO.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ IO.o

control_main: $(OBJ) control_main.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) control_main.o

control_test: $(OBJ) control_test.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) control_test.o

control_main_full: $(OBJ) control_main_full.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) control_main_full.o


dd_main: $(OBJ) dd_main.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) dd_main.o

dd_main_full: $(OBJ) dd_main_full.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) dd_main_full.o


read_parallel: read_parallel.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ read_parallel.o

dd_gmres: $(OBJ) dd_gmres.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) dd_gmres.o

dd_gmres_full: $(OBJ) dd_gmres_full.o
	$(CC) -Wl,--no-as-needed -L$(LIBRARY_PATH) -I$(INCLUDE_PATH) $(MPI_LINK_FLAGS) -o $@ $(OBJ) dd_gmres_full.o
