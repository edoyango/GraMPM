module grid_m

    use params, only: f

    implicit none

    type:: grid
        integer:: ngridx, ngridy, ngridz
        real(f), allocatable, dimension(:, :, :):: mass, momentum, force, v
    contains
        procedure:: init
    end type grid

contains

    subroutine init(self, ngridx, ngridy, ngridz)

        integer, intent(in):: ngridx, ngridy, ngridz
        class(grid), intent(inout):: self

        allocate(self%mass(ngridx, ngridy, ngridz), self%momentum(ngridx, ngridy, ngridz), &
            self%force(ngridx, ngridy, ngridz), self%v(ngridx, ngridy, ngridz), source=0._f)

    end subroutine init

end module grid_m