program test_inits

    use iso_fortran_env, only: output_unit

    implicit none

    write(output_unit, "(A)") "Testing particle inits..."
    call test_particle_init()

    write(output_unit, "(A)") "Testing grid inits..."
    call test_grid_init()

contains

    character(len=4) function pass_or_fail(tf)

        logical, intent(in):: tf

        if (tf) then
            pass_or_fail = "PASS"
        else
            pass_or_fail = "FAIL"
        end if

    end function pass_or_fail

    subroutine test_particle_init()
        ! test that all the members of the particles type are intialized and relevant attributes have values set

        use params, only: dx, maxn, f
        use particles_m, only: particles

        type(particles):: p

        call p%init(maxn)

        write(output_unit, "(4x, A)") "p%ntotal == 0"
        write(output_unit, "(4x, A)") pass_or_fail(p%ntotal == 0)
        write(output_unit, "(4x, A)")  "p%maxn == maxn"
        write(output_unit, "(4x, A)")  pass_or_fail(p%maxn == maxn)
        write(output_unit, "(4x, A)")  "size(p%x) == maxn .and. all(p%x == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%x) == maxn .and. all(p%x == 0._f))
        write(output_unit, "(4x, A)")  "size(p%y) == maxn .and. all(p%y == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%y) == maxn .and. all(p%y == 0._f))
        write(output_unit, "(4x, A)")  "size(p%z) == maxn .and. all(p%z == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%z) == maxn .and. all(p%z == 0._f))
        write(output_unit, "(4x, A)")  "size(p%vx) == maxn .and. all(p%vx == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%vx) == maxn .and. all(p%vx == 0._f))
        write(output_unit, "(4x, A)")  "size(p%vy) == maxn .and. all(p%vy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%vy) == maxn .and. all(p%vy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%vz) == maxn .and. all(p%vz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%vz) == maxn .and. all(p%vz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%sigxx) == maxn .and. all(p%sigxx == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%sigxx) == maxn .and. all(p%sigxx == 0._f))
        write(output_unit, "(4x, A)")  "size(p%sigyy) == maxn .and. all(p%sigyy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%sigyy) == maxn .and. all(p%sigyy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%sigzz) == maxn .and. all(p%sigzz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%sigzz) == maxn .and. all(p%sigzz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%sigxy) == maxn .and. all(p%sigxy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%sigxy) == maxn .and. all(p%sigxy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%sigxz) == maxn .and. all(p%sigxz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%sigxz) == maxn .and. all(p%sigxz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%sigyz) == maxn .and. all(p%sigyz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%sigyz) == maxn .and. all(p%sigyz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%mass) == maxn .and. all(p%mass == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%mass) == maxn .and. all(p%mass == 0._f))
        write(output_unit, "(4x, A)")  "size(p%rho) == maxn .and. all(p%rho == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%rho) == maxn .and. all(p%rho == 0._f))
        write(output_unit, "(4x, A)")  "size(p%ax) == maxn .and. all(p%ax == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%ax) == maxn .and. all(p%ax == 0._f))
        write(output_unit, "(4x, A)")  "size(p%ay) == maxn .and. all(p%ay == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%ay) == maxn .and. all(p%ay == 0._f))
        write(output_unit, "(4x, A)")  "size(p%az) == maxn .and. all(p%az == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%az) == maxn .and. all(p%az == 0._f))
        write(output_unit, "(4x, A)")  "size(p%dx) == maxn .and. all(p%dx == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%dx) == maxn .and. all(p%dx == 0._f))
        write(output_unit, "(4x, A)")  "size(p%dy) == maxn .and. all(p%dy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%dy) == maxn .and. all(p%dy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%dz) == maxn .and. all(p%dz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%dz) == maxn .and. all(p%dz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%epsxx) == maxn .and. all(p%epsxx == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%epsxx) == maxn .and. all(p%epsxx == 0._f))
        write(output_unit, "(4x, A)")  "size(p%epsyy) == maxn .and. all(p%epsyy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%epsyy) == maxn .and. all(p%epsyy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%epszz) == maxn .and. all(p%epszz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%epszz) == maxn .and. all(p%epszz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%epsxy) == maxn .and. all(p%epsxy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%epsxy) == maxn .and. all(p%epsxy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%epsxz) == maxn .and. all(p%epsxz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%epsxz) == maxn .and. all(p%epsxz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%epsyz) == maxn .and. all(p%epsyz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%epsyz) == maxn .and. all(p%epsyz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%omegaxy) == maxn .and. all(p%omegaxy == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%omegaxy) == maxn .and. all(p%omegaxy == 0._f))
        write(output_unit, "(4x, A)")  "size(p%omegaxz) == maxn .and. all(p%omegaxz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%omegaxz) == maxn .and. all(p%omegaxz == 0._f))
        write(output_unit, "(4x, A)")  "size(p%omegayz) == maxn .and. all(p%omegayz == 0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(p%omegayz) == maxn .and. all(p%omegayz == 0._f))

    end subroutine test_particle_init

    subroutine test_grid_init() 

        use params, only: f, n_buffer_cells
        use grid_m, only: grid

        type(grid):: g
        real(f), parameter:: dcell = 0.1_f, mingrid(3) = [-0.2_f, 0.1_f, 1._f], &
            maxgrid(3) = [0.1_f, 1._f, 2.05_f]
        integer, parameter:: correct_ngrid(3) = [4, 10, 12] + 2*n_buffer_cells, &
            ncells = product(correct_ngrid)

        ! x-direction
        ! -0.4 -0.3 -0.2 -0.1  0.0  0.1  0.2  0.3
        ! |    |    |    |    |    |    |    |
        ! +----+----+----+----+----+----+----+
        ! 1    2    3    4    5    6    7    8

        ! y-direction
        ! -0.1  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0  1.1  1.2
        ! |    |    |    |    |    |    |    |    |    |    |    |    |    |
        ! +----+----+----+----+----+----+----+----+----+----+----+----+----+
        ! 1    2    3    4    5    6    7    8    9    10   11   12   13   14

        ! z-direction
        ! 0.8  0.9  1.0  1.1  1.2  1.3  1.4  1.5  1.6  1.7  1.8  1.9  2.0  2.1  2.2  2.3
        ! |    |    |    |    |    |    |    |    |    |    |    |    |    |    |    |
        ! +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        ! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16

        call g%init(mingrid, maxgrid, dcell)

        write(output_unit, "(4x, A)")  "g%dcell == dcell"
        write(output_unit, "(4x, A)")  pass_or_fail(g%dcell == dcell)
        write(output_unit, "(4x, A)")  "all(g%mingrid == mingrid-dcell*n_buffer_cells)"
        write(output_unit, "(4x, A)")  pass_or_fail(all(g%mingrid == mingrid-dcell*n_buffer_cells))
        write(output_unit, "(4x, A)")  "all(g%maxgrid == maxgrid+dcell*n_buffer_cells)"
        write(output_unit, "(4x, A)")  pass_or_fail(all(g%maxgrid == maxgrid+dcell*n_buffer_cells))
        write(output_unit, "(4x, A)")  "all(g%ngrid == correct_ngrid) "
        write(output_unit, "(4x, A)")  pass_or_fail(all(g%ngrid == correct_ngrid))
        write(output_unit, "(4x, A)")  "size(g%mass) == ncells .and. all(g%mass==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%mass) == ncells .and. all(g%mass==0._f))
        write(output_unit, "(4x, A)")  "size(g%momentumx) == ncells .and. all(g%momentumx==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%momentumx) == ncells .and. all(g%momentumx==0._f))
        write(output_unit, "(4x, A)")  "size(g%momentumy) == ncells .and. all(g%momentumy==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%momentumy) == ncells .and. all(g%momentumy==0._f))
        write(output_unit, "(4x, A)")  "size(g%momentumz) == ncells .and. all(g%momentumz==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%momentumz) == ncells .and. all(g%momentumz==0._f))
        write(output_unit, "(4x, A)")  "size(g%forcex) == ncells .and. all(g%forcex==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%forcex) == ncells .and. all(g%forcex==0._f))
        write(output_unit, "(4x, A)")  "size(g%forcey) == ncells .and. all(g%forcey==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%forcey) == ncells .and. all(g%forcey==0._f))
        write(output_unit, "(4x, A)")  "size(g%forcez) == ncells .and. all(g%forcez==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%forcez) == ncells .and. all(g%forcez==0._f))
        write(output_unit, "(4x, A)")  "size(g%vx) == ncells .and. all(g%vx==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%vx) == ncells .and. all(g%vx==0._f))
        write(output_unit, "(4x, A)")  "size(g%vy) == ncells .and. all(g%vy==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%vy) == ncells .and. all(g%vy==0._f))
        write(output_unit, "(4x, A)")  "size(g%vz) == ncells .and. all(g%vz==0._f)"
        write(output_unit, "(4x, A)")  pass_or_fail(size(g%vz) == ncells .and. all(g%vz==0._f))

    end subroutine test_grid_init

end program