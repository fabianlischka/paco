CC = mpicc
MPI_COMPILE_FLAGS = -O2 -Wl,--no-as-needed
MPI_LINK_FLAGS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIBRARY_PATHS = -L/u1/local/intel/mkl/lib/intel64 -L../lib
INCLUDE_PATHS = -I/u1/local/intel/mkl/include
DEPS = control.h
OBJ = Problem.o Solver.o IO.o
TARG = control_test control_main control_main_full \
       						  dd_main      dd_main_full      \
										dd_gmres     dd_gmres_full   read_parallel

all: $(TARG)

%.o: %.c $(DEPS)
	$(CC) $(LIBRARY_PATHS) $(INCLUDE_PATHS) $(MPI_COMPILE_FLAGS) -c -o $@ $<

$(TARG): % : %.o $(OBJ)
	$(CC) -Wl,--no-as-needed $(LIBRARY_PATHS) $(INCLUDE_PATHS) $(MPI_LINK_FLAGS) -o $@ $(OBJ) $<

.PHONY : clean
clean:
	-rm *.o $(TARG)
