!!!!
!! File: cuda-dgemv.f90
!! Description: CUDA BLAS (cublas) DGEMV performance analysis (Thanks Davide del Vento @ UCAR)
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Friday, 17th September 2021, 11:24:19 am
!! Last Modified: Friday, 17th September 2021, 11:31:19 am
!!  
!! Copyright (c) 2021, Bruno R. de Abreu, National Center for Supercomputing Applications.
!! All rights reserved.
!! License: This program and the accompanying materials are made available to any individual
!!          under the citation condition that follows: On the event that the software is
!!          used to generate data that is used implicitly or explicitly for research
!!          purposes, proper acknowledgment must be provided in the citations section of
!!          publications. This includes both the author's name and the National Center
!!          for Supercomputing Applications. If you are uncertain about how to do
!!          so, please check this page: https://github.com/babreu-ncsa/cite-me.
!!          This software cannot be used for commercial purposes in any way whatsoever.
!!          Omitting this license when redistributing the code is strongly disencouraged.
!!          The software is provided without warranty of any kind. In no event shall the
!!          author or copyright holders be liable for any kind of claim in connection to
!!          the software and its usage.
!!!!

module mycublas
    interface
      subroutine launchCUDGEMV(m,n,alpha,A,lda,x,incx,beta,y,incy,stream) BIND (C, NAME='launchCUDGEMV')
        use ISO_C_BINDING
        implicit none
        integer (C_INT), value :: m
        integer (C_INT), value :: n
        double precision (C_DOUBLE), value :: alpha
        type (C_PTR), value :: A
        integer (C_INT), value :: lda
        type (C_PTR), value :: x
        integer (C_INT), value :: incx
        double precision (C_DOUBLE), value :: beta
        type (C_PTR), value :: y
        integer (C_INT), value :: incy
        integer (C_INT), value :: stream
      end subroutine
    end interface
  end module mycublas
  
  program   MAIN
  
  use ISO_C_BINDING
  use mycublas
  use openacc
  use iso_fortran_env
  
  implicit none
  
  double precision(C_DOUBLE) :: ALPHA, BETA
  integer          M, P, I, J, nloops
  parameter        (M=10000, P=10000)
  parameter        (nloops=100)
  double precision(C_DOUBLE) :: A(M,P), x(P), y(M)
  integer(C_INT) :: stream
  double precision :: startT, endT
  
  ALPHA = 1.0 
  BETA = 0.0
  
  call random_seed()
  call random_number(A)
  call random_number(x)
  y = 0.0
  
  call cpu_time(startT)
  ! Copy data to device at start of region and back to host and end of region
  !$acc data copy(A, x, y)
  
      ! Inside this region the device data pointer will be used
      !$acc host_data use_device(A, x, y)
          stream = acc_get_cuda_stream(acc_async_sync)
          do = 1, nloops
            call launchCUDGEMV(M,P,ALPHA,        &
                             C_LOC(A),M,       &
                             C_LOC(x),1,BETA,  &
                             C_LOC(y),1,stream)
          enddo
      !$acc end host_data
  
  !$acc end data
  call cpu_time(endT)
  
  PRINT *, "== Matrix-vector multiplication using CUBLAS =="
  PRINT 50, " == completed at ",(endT-startT)*1000/nloops," milliseconds =="
50  FORMAT(A,F12.5,A)
  PRINT *, ""
  
  end program
  
