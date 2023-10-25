program test_p2g_idx

    use params, only: f
    use pass_or_fail_m, only: pass_or_fail
    use mpm_class_m, only: mpm, x2idx

    implicit none
    type(mpm):: mympm
    real(f), parameter:: dcell = 0.2_f, mingrid(3) = [-1._f, 0.05_f, 1._f], maxgrid(3) = [1._f, 0.3_f, 2.5_f]

    call mympm%p%init(5)

    mympm%p%x(:) = [-0.51_f, -0.58_f, -0.26_f, -0.46_f, -0.64_f]
    mympm%p%y(:) = [0.11, 0.21,	0.17, 0.07,	0.18]
    mympm%p%z(:) = [1.20, 1.07, 1.01, 1.66, 1.29]

    call test_x2idx(mympm)

contains

    subroutine test_x2idx(mpm_obj)
        
        type(mpm):: mpm_obj

        ! x-direction
        ! -1.0  -0.8  -0.6  -0.4  -0.2  0.0   0.2   0.4
        ! |     |     |     |     |     |     |     |
        ! x--1--x--2--x--3--x--4--x--5--x--6--x--7--x

        ! y-direction
        ! 0.05  0.25  0.45
        ! |     |     |
        ! x--1--x--2--x

        ! z-direction
        ! 1.0   1.2   1.4   1.6   1.8   2.0   2.2   2.4   2.6
        ! |     |     |     |     |     |     |     |     |
        ! x--1--x--2--x--3--x--4--x--5--x--6--x--7--x--8--x
        mympm%p%idx = x2idx(mympm%p%x, mingrid(1), dcell)
        mympm%p%idy = x2idx(mympm%p%y, mingrid(2), dcell)
        mympm%p%idz = x2idx(mympm%p%z, mingrid(3), dcell)
        call pass_or_fail("Correct idx", all(mympm%p%idx == [3, 3, 4, 3, 2]))
        call pass_or_fail("Correct idy", all(mympm%p%idy == [1, 1, 1, 1, 1]))
        call pass_or_fail("Correct idz", all(mympm%p%idz == [2, 1, 1, 4, 2]))

    end subroutine test_x2idx

end program test_p2g_idx