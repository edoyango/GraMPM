program mpm_main

    use params, only: f
    use mpm_class_m, only: mpm

    implicit none
    type(mpm):: mympm
    integer, parameter:: dims=3, npx=50, npy=10, npz=25, ntotal=npx*npy*npz
    real(f), parameter:: dx=0.002_f, dcell=0.004_f, mingrid(dims) = [-dcell*2._f, -dcell*2._f, -dcell*2._f], &
        maxgrid(dims) = [0.3_f+2._f*dcell, 0.1_f, 0.02_f+2._f*dcell], rho_ini=2650._f, mass_ini=rho_ini*dx*dx*dx
    integer, parameter:: ngrid(dims) = int((maxgrid-mingrid)/dcell)+1
    integer:: i, j, k, n

    call mympm%p%init(ntotal)
    call mympm%g%init(ngrid(1), ngrid(2), ngrid(3))

    n = 0
    do i = 1, npx
        do j = 1, npy
            do k = 1, npz
                n = n + 1
                mympm%p%x = (i-0.5_f)*dx
                mympm%p%y = (j-0.5_f)*dx
                mympm%p%z = (k-0.5_f)*dx
            end do
        end do
    end do

    mympm%p%rho(1:ntotal) = rho_ini
    mympm%p%mass(1:ntotal) = mass_ini

end program mpm_main