SHELL = /bin/sh

FC=mpif90
LDLIBS= 
FFLAGS= -g -O3 -mcmodel=large  
PP=  param.o adaptint.o func.o

# object files
OBJS=  subr_support.o subr_readfile.o read_halofile.o  calpha.o c_alpha_mpi.o 

run: c_alpha


c_alpha:$(PP) $(OBJS)
	$(FC) $(FFLAGS) -openmp  $^ -o $@  $(LDLIBS)

$(OBJS): $(PP)

%.o: %.f90
	$(FC) $(FFLAGS)  -c $< $(LDLIBS)

clean:
	rm -f *.o *.mod *~ c_alpha
	find . -type f | xargs -n 5 touch
	rm -rf $(OBJS)
