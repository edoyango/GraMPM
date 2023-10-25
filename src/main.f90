program mpm_main

    use params, only: f, dims, maxn, dx, n_buffer_cells, dcell
    use mpm_class_m, only: mpm

    implicit none
    type(mpm):: mympm
    real(f), parameter:: mingrid(dims) = 0._f, maxgrid(dims) = [0.3_f, 0.1_f, 0.02_f], rho_ini=2650._f, mass_ini=rho_ini*dx*dx*dx
    integer, parameter:: ngrid(dims) = int((maxgrid-mingrid)/dcell)+1
    ! integer:: i, j, k, n

    call mympm%p%init(maxn)
    call mympm%g%init(mingrid, maxgrid, dcell)

    call init_particles(mympm%p)

contains

    subroutine init_particles(p)

        use particles_m, only: particles

        implicit none
        type(particles), intent(inout):: p
        integer:: n, i, j, k
        integer, parameter:: npx=50, npy=10, npz=25

        p%ntotal = 0
        do i = 1, npx
            do j = 1, npy
                do k = 1, npz
                    call p%append_position( &
                        (i-0.5_f)*dx, &
                        (j-0.5_f)*dx, &
                        (k-0.5_f)*dx &
                    )
                end do
            end do
        end do

        p%rho(1:p%ntotal) = rho_ini
        p%mass(1:p%ntotal) = mass_ini

    end subroutine init_particles

end program mpm_main