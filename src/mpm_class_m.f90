module mpm_class_m

    use params, only: f
    use particles_m, only: particles
    use grid_m, only: grid

    type:: mpm
        type(particles):: p
        type(grid):: g
        integer, allocatable:: pgpairs(:, :)
    contains
        procedure:: find_particle_gridnode_pairs
    end type mpm

contains

    elemental integer function x2idx(x, mingridx, dcell)

        real(f), intent(in):: x, mingridx, dcell

        x2idx = int((x-mingridx)/dcell) + 1

    end function x2idx

    subroutine find_particle_gridnode_pairs(self)

        class(mpm), intent(inout):: self

        self%p%idx(1:self%p%ntotal) = x2idx(self%p%x(1:self%p%ntotal), self%g%mingrid(1), self%g%dcell)
        self%p%idy(1:self%p%ntotal) = x2idx(self%p%y(1:self%p%ntotal), self%g%mingrid(2), self%g%dcell)
        self%p%idz(1:self%p%ntotal) = x2idx(self%p%z(1:self%p%ntotal), self%g%mingrid(3), self%g%dcell)

    end subroutine find_particle_gridnode_pairs

end module mpm_class_m