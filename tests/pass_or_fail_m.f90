module pass_or_fail_m

    use iso_fortran_env, only: output_unit

contains

    subroutine pass_or_fail(teststr, tf)

        character(len=*), intent(in):: teststr
        logical, intent(in):: tf
        character(len=4):: outcome

        if (tf) then
            outcome = "PASS"
        else
            outcome = "FAIL"
        end if

        write(output_unit, "(A, 4x, A)") outcome, trim(teststr)

    end subroutine pass_or_fail

end module pass_or_fail_m