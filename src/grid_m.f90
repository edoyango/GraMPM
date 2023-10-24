module grid_m

    use params, only: f, n_buffer_cells

    implicit none

    type:: grid
        integer:: ngrid(3)
        real(f):: mingrid(3), maxgrid(3), dcell
        real(f), allocatable, dimension(:, :, :):: mass, momentumx, momentumy, momentumz, forcex, forcey, forcez, vx, &
            vy, vz
    contains
        procedure:: init
    end type grid

contains

    subroutine init(self, mingrid, maxgrid, dcell)

        real(f), intent(in):: mingrid(3), maxgrid(3), dcell
        class(grid), intent(inout):: self

        self%dcell = dcell
        self%mingrid = mingrid - n_buffer_cells*dcell
        self%maxgrid = maxgrid + n_buffer_cells*dcell
        self%ngrid = ceiling((self%maxgrid(:)-self%mingrid(:))/dcell) + 1

        allocate(self%mass(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%momentumx(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%momentumy(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%momentumz(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%forcex(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%forcey(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%forcez(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%vx(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%vy(self%ngrid(1), self%ngrid(2), self%ngrid(3)), &
            self%vz(self%ngrid(1), self%ngrid(2), self%ngrid(3)), source=0._f)

    end subroutine init

end module grid_m