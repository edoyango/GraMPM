module particles_m

    use params, only: f

    implicit none

    type:: particles
        integer:: maxn
        real(f), allocatable, dimension(:):: x, y, z, vx, vy, vz, sigxx, sigyy, sigzz, sigxy, sigxz, sigyz, mass, rho, &
            ax, ay, az, dx, dy, dz, epsxx, epsyy, epszz, epsxy, epsxz, epsyz, omegaxy, omegaxz, omegayz
    contains
        procedure:: init
    end type particles

contains

    subroutine init(self, maxn)

        integer, intent(in):: maxn
        class(particles), intent(inout):: self

        self%maxn = maxn

        allocate(self%x(maxn), self%y(maxn), self%z(maxn), self%vx(maxn), self%vy(maxn), self%vz(maxn), &
            self%sigxx(maxn), self%sigyy(maxn), self%sigzz(maxn), self%sigxy(maxn), self%sigxz(maxn), &
            self%sigyz(maxn), self%mass(maxn), self%rho(maxn), self%ax(maxn), self%ay(maxn), self%az(maxn), &
            self%dx(maxn), self%dy(maxn), self%dz(maxn), self%epsxx(maxn), self%epsyy(maxn), self%epszz(maxn), &
            self%epsxy(maxn), self%epsxz(maxn), self%epsyz(maxn), self%omegaxy(maxn), self%omegaxz(maxn), &
            self%omegayz(maxn), source=0._f)

    end subroutine init

end module particles_m