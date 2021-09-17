# this is configured for SDSC's Expanse
BLIS_PREFIX	= /cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/gcc-10.2.0/amdblis-2.2-jcmqdhq7cozl3yw3xkvcec4xsbt7o7kl
BLIS_INC	= $(BLIS_PREFIX)/include/blis
BLIS_LIB	= $(BLIS_PREFIX)/lib/libblis-mt.a

OTHER_LIBS	= -lm -lpthread -fopenmp

FC		= gfortran
CFLAGS		= -I$(BLIS_INC)
LINKER		= $(FC)

OBJS		= amd-dgemm.o

%.o: %.f90
	$(FC) $(CFLAGS) -c $< -o $@

all: $(OBJS)
	$(LINKER) $(OBJS) $(BLIS_LIB) $(OTHER_LIBS) -o amd-dgemm.x

clean:
	rm $(OBJS) *.x