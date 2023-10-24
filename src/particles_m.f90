module particles_m

    use params, only: f

    implicit none

    type:: particles
        integer:: maxn, ntotal
        real(f), allocatable, dimension(:):: x, y, z, vx, vy, vz, sigxx, sigyy, sigzz, sigxy, sigxz, sigyz, mass, rho, &
            ax, ay, az, dx, dy, dz, epsxx, epsyy, epszz, epsxy, epsxz, epsyz, omegaxy, omegaxz, omegayz
    contains
        procedure:: init, append_position
    end type particles

contains

    subroutine init(self, maxn)

        integer, intent(in):: maxn
        class(particles), intent(inout):: self

        self%maxn = maxn
        self%ntotal = 0

        allocate(self%x(maxn), self%y(maxn), self%z(maxn), self%vx(maxn), self%vy(maxn), self%vz(maxn), &
            self%sigxx(maxn), self%sigyy(maxn), self%sigzz(maxn), self%sigxy(maxn), self%sigxz(maxn), &
            self%sigyz(maxn), self%mass(maxn), self%rho(maxn), self%ax(maxn), self%ay(maxn), self%az(maxn), &
            self%dx(maxn), self%dy(maxn), self%dz(maxn), self%epsxx(maxn), self%epsyy(maxn), self%epszz(maxn), &
            self%epsxy(maxn), self%epsxz(maxn), self%epsyz(maxn), self%omegaxy(maxn), self%omegaxz(maxn), &
            self%omegayz(maxn), source=0._f)

    end subroutine init

    subroutine append_position(self, x, y, z)

        class(particles), intent(inout):: self
        real(f), intent(in):: x, y, z

        self%ntotal = self%ntotal + 1
        self%x(self%ntotal) = x
        self%y(self%ntotal) = y
        self%z(self%ntotal) = z

    end subroutine append_position

end module particles_m