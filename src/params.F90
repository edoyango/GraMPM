module params

    use iso_fortran_env

    implicit none
#ifdef SP
    integer, parameter:: f = real32
#else
    integer, parameter:: f = real64
#endif

end module params