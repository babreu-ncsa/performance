CUDACC=nvcc
CUDAFC=nvfortran
CFLAGS=-c
FFLAGS=-acc
CUDAFLAGS=-L${CUDA_HOME}/lib64 -lcublas -lcudart

CUDAFILE=mycublas.cu

cuda-dgemv: cuda-dgemv.f90
	${CUDACC} ${CFLAGS} ${CUDAFILE}
	${CUDAFC} ${FFLAGS} ${CUDAFLAGS} *.o $< -o $@

clean:
	rm *.o *.mod cuda-dgemv