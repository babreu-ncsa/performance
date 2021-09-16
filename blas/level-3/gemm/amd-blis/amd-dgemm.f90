!!!!
!! File: amd-dgemm.f90
!! Description: AMD BLIS DGEMM performance analysis A*B = C
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Thursday, 16th September 2021, 12:10:24 pm
!! Last Modified: Thursday, 16th September 2021, 12:14:58 pm
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

program amd_dgemm
    use, intrinsic :: iso_fortran_env
    use :: omp_lib
    implicit none
    integer, parameter :: dp = REAL64 ! double precision float
    integer, parameter :: i32 = INT32 ! 32-bit integer
    integer(i32), parameter :: ord1=10000_i32  ! leading dim of matrix A
    integer(i32), parameter :: ord2=10000_i32   ! lower dim of matrix A
    integer(i32), parameter :: ord3=10000_i32   ! other dim of B
    integer(i32), parameter :: nloops=100_i32   ! number of DGEMM calls
    real(dp) :: startT, endT
    real(dp), dimension(:,:), allocatable :: m, v, p ! (m*v=p)
    integer(i32) :: i

    ! allocate
    allocate(m(ord1,ord2))
    allocate(v(ord2,ord3))
    allocate(p(ord1,ord3))

    ! fill in with random stuff
    call random_seed()
    call random_number(m)
    call random_number(v)
    p = 0.0_dp

    ! call AMD BLIS (syntax below, first call usually a query)
    ! dgemm('N', 'N', M, N, K, ALPHA, A, M, B, K, BETA, C, M)
    call dgemm('N', 'N', ord1, ord3, ord2, 1.0_dp, m, ord1, v, ord2, 0.0_dp, p, ord1)

    ! now time it
    startT = omp_get_wtime()
    do i = 1, nloops
            call dgemm('N', 'N', ord1, 1, ord2, 1.0_dp, m, ord1, v, ord2, 0.0_dp, p, ord1)
    enddo
    endT = omp_get_wtime()

    PRINT *, "== Matrix multiplication using AMD BLIS DGEMM =="
    PRINT 50, "== completed at ",(endT-startT)*1000/nloops," milliseconds =="
50   FORMAT(A,F12.5,A)
    PRINT *, ""

end program amd_dgemm

