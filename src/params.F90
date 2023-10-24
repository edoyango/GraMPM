module params

    use iso_fortran_env

    implicit none
#ifdef SP
    integer, parameter:: f = real32
#else
    integer, parameter:: f = real64
#endif

    ! geometric parameters
    integer, parameter:: dims = 3, maxn = 5000
    real(f), parameter:: dx = 0.002_f

    ! grid dimensions
    integer, parameter:: n_buffer_cells = 2
    real(f), parameter:: dcell = 2._f*dx

end module params