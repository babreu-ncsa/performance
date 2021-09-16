FC=ifort
MKLPATH=${MKLROOT}/lib/intel64
MKLINCLUDE=${MKLROOT}/mkl/include

LDFLAGS=-mkl=parallel -L${MKLPATH} -I${MKLINCLUDE} -I${MKLINCLUDE}/intel64/lp64 -lmkl_lapack95_lp64 -Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_intel_thread.a ${MKLPATH}/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm

SOURCE=mkl-dgemm.f90

all:
	$(FC) $(SOURCE) $(LDFLAGS)