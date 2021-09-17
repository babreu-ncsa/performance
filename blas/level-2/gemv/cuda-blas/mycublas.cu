/***
* File: mycublas.cu
* Description: C code with call to cuBLAS Dgemm (Thanks Davide del Vento @ UCAR)
* Author: Bruno R. de Abreu  |  babreu at illinois dot edu
* National Center for Supercomputing Applications (NCSA)
*  
* Creation Date: Friday, 17th September 2021, 11:26:10 am
* Last Modified: Friday, 17th September 2021, 11:27:32 am
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

// Declared extern "C" to disable C++ name mangling
extern "C" void launchCUDGEMV(
    int m, int n,
    const double alpha,
    const double *A, int lda,
    const double *x, int incx,
    const double beta,
    double *y, int incy,
    void *stream)
{
    cublasHandle_t handle;
    cublasCreate(&handle); // this needs to match the destroy below
    cublasSetStream(handle, (cudaStream_t)stream);

    cublasOperation_t trans = CUBLAS_OP_N; // supporting only 'N' for now
    cublasDgemv(handle, trans,
                m, n,
                &alpha,
                A, lda,
                x, incx,
                &beta,
                y, incy);

    cublasDestroy(handle); // this needs to match the create above
                           // for performance, minimize the number of cublasCreate()/cublasDestroy()
}
