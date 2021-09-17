!!!!
!! File: mkl-dgemv.f90
!! Description: Intel MKL's DGEMV performance analysis
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Friday, 17th September 2021, 10:40:58 am
!! Last Modified: Friday, 17th September 2021, 10:48:15 am
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

program mkl_dgemv
    use, intrinsic :: iso_fortran_env
    use :: mkl_service
    implicit none
    include "mkl_lapack.fi"
    integer, parameter :: dp = REAL64 ! double precision float
    integer, parameter :: i32 = INT32 ! 32-bit integer
    integer(i32), parameter :: ord1=10000_i32  ! leading dim of matrix m
    integer(i32), parameter :: ord2=10000_i32   ! dim of vector v
    integer(i32), parameter :: nloops=100       ! number of times to call DGEMV
    real(dp) :: startT, endT
    real(dp), dimension(:,:), allocatable :: m  ! (m*v = p)
    real(dp), dimension(:), allocatable :: v, p
    integer(i32) :: i, MAX_THREADS

    ! allocate
    allocate(m(ord1,ord2))
    allocate(v(ord2))
    allocate(p(ord1))

    ! fill in with random stuff
    call random_seed()
    call random_number(m)
    call random_number(v)
    p = 0.0_dp

    ! set correct number of threads
    MAX_THREADS = MKL_GET_MAX_THREADS()
    PRINT 20," Running Intel(R) MKL with ",MAX_THREADS," threads"
20   FORMAT(A,I2,A)
    PRINT *, ""
    CALL MKL_SET_NUM_THREADS(MAX_THREADS)
    ! call MKL (syntax below), this is just to warm things up
    ! dgemv('N', M, N, alpha, a, lda, x, incx, beta, y, incy)
    call dgemv('N', ord1, ord2, 1.0_dp, m, ord1, v, 1, 0.0_dp, p, 1)

    ! not time it for real
    startT = dsecnd()
    do i = 1, nloops
        call dgemv('N', ord1, ord2, 1.0_dp, m, ord1, v, 1, 0.0_dp, p, 1)
    enddo
    endT = dsecnd()

    PRINT *, "== Matrix-vector multiplication using Intel(R) MKL DGEMV =="
    PRINT 50, " == completed at ",(endT-startT)*1000/nloops," milliseconds =="
50  FORMAT(A,F12.5,A)
    PRINT *, ""

end program mkl_dgemv
