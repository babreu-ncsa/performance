/***
* File: mycublas.cu
* Description: C code with call to cuBLAS Dgemm 
* Author: Bruno R. de Abreu  |  babreu at illinois dot edu
* National Center for Supercomputing Applications (NCSA)
*  
* Creation Date: Friday, 17th September 2021, 9:27:08 am
* Last Modified: Friday, 17th September 2021, 10:06:03 am
*  
* Copyright (c) 2021, Bruno R. de Abreu, National Center for Supercomputing Applications.
* All rights reserved.
* License: This program and the accompanying materials are made available to any individual
*          under the citation condition that follows: On the event that the software is
*          used to generate data that is used implicitly or explicitly for research
*          purposes, proper acknowledgment must be provided in the citations section of
*          publications. This includes both the author's name and the National Center
*          for Supercomputing Applications. If you are uncertain about how to do
*          so, please check this page: https://github.com/babreu-ncsa/cite-me.
*          This software cannot be used for commercial purposes in any way whatsoever.
*          Omitting this license when redistributing the code is strongly disencouraged.
*          The software is provided without warranty of any kind. In no event shall the
*          author or copyright holders be liable for any kind of claim in connection to
*          the software and its usage.
***/

#include "cublas_v2.h"

// Declared exter "C" to disable C++ name mangling
extern "C" void launchCUDGEMM(
    int m, int n, int k,
    const double alpha,
    const double *A, int lda,
    const double *B, int ldb,
    const double beta,
    double *C, int ldc,
    void *stream)
{
    cublasHandle_t handle;
    cublasCreate(&handle); // create
    cublasSetStream(handle, (cudaStream_t)stream);

    cublasOperation_t trans = CUBLAS_OP_N; // transpositon: only 'N'
    cublasDgemm(handle, trans, trans,
                m, n, k,
                &alpha,
                A, lda,
                B, ldb,
                &beta,
                C, ldc);

    cublasDestroy(handle); // destroy
}
