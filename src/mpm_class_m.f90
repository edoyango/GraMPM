module mpm_class_m

    use params, only: f
    use particles_m, only: particles
    use grid_m, only: grid

    type:: mpm
        type(particles):: p
        type(grid):: g
    end type mpm

end module mpm_class_m