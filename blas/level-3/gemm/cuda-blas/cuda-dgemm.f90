!!!!
!! File: cuda-dgemm.f90
!! Description: CUDA BLAS (cublas) DGEMM performance analysis A*B = C
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Friday, 17th September 2021, 7:22:49 am
!! Last Modified: Friday, 17th September 2021, 10:16:33 am
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
        subroutine launchCUDGEMM(M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, stream) BIND (C, NAME='launchCUDGEMM')
            use, intrinsic :: iso_c_binding
            implicit none
            integer (C_INT), value :: M
            integer (C_INT), value :: N
            integer (C_INT), value :: K
            double precision (C_DOUBLE), value :: ALPHA
            type (C_PTR), value :: A
            integer (C_INT), value :: LDA
            type (C_PTR), value :: B
            integer (C_INT), value :: LDB
            double precision (C_DOUBLE), value :: BETA
            type (C_PTR), value :: C
            integer (C_INT), value :: LDC
            integer (C_INT), value :: stream 
        end subroutine
    end interface
end module mycublas


program MAIN
    use, intrinsic :: iso_c_binding ! this translates data types from FORTRAN to C
                        ! https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html
    use openacc         ! OpenACC API
    use, intrinsic :: iso_fortran_env
    use mycublas
    implicit none
    double precision(C_DOUBLE) :: alpha, beta
    integer, parameter :: ord1=10000  ! leading dim of matrix A
    integer, parameter :: ord2=10000   ! lower dim of matrix A
    integer, parameter :: ord3=10000   ! other dim of B
    integer, parameter :: nloops=100   ! number of DGEMM calls
    double precision :: startT, endT
    double precision(C_DOUBLE),dimension(:,:),allocatable :: m, v, p ! (m*v=p)
    integer :: i
    integer(C_INT) :: stream

    ! allocate
    allocate(m(ord1,ord2))
    allocate(v(ord2,ord3))
    allocate(p(ord1,ord3))

    ! fill in with random stuff
    call random_seed()
    call random_number(m)
    call random_number(v)
    p = 0.0

    alpha = 1.0
    beta = 0.0

    call cpu_time(startT)
    ! copy data to device at beginning and back to host at the end
    !$acc data copy(m,v,p)
        !$acc host_data use_device(m,v,p)
            stream = acc_get_cuda_stream(acc_async_sync)
            do i = 1, nloops
                call launchCUDGEMM(     &
                    ord1, ord2, ord3,   &
                    alpha,              &
                    C_LOC(m), ord1,     &
                    C_LOC(v), ord2,     &
                    beta,               &
                    C_LOC(p), ord1,     &
                    stream)
            enddo
        !$acc end host_data
    !$acc end data
    call cpu_time(endT)

    PRINT *, "== Matrix multiplication using CUBLAS DGEMM =="
    PRINT 50, "== completed at ",(endT-startT)*1000/nloops," milliseconds =="
50   FORMAT(A,F12.5,A)
    PRINT *, ""

    deallocate(m)
    deallocate(v)
    deallocate(p)

end program
