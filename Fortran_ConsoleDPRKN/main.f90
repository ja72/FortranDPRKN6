program main
    use mod_common
    use mod_ode_system
    use mod_harmonic_system
    implicit none

    call run_harmonic_test()

contains

    subroutine run_harmonic_test()
        real(wp) :: h, tf
        real(wp), allocatable :: y_exact(:)
        integer :: target_cycles, resolution
        integer, parameter :: system_size = 1
        type(harmonic_system), target :: sys
        type(ode_solver(system_size)), pointer :: solver

        print *, "==========================================================="
        print *, "Testing mod_ode with Harmonic Oscillator Example - Fortran"
        print *, ""

        sys = harmonic_system(omega=1.0_wp, y_init=[1.0_wp], dy_init=[0.0_wp])
        solver => init_ode_solver(sys)

        target_cycles = 10
        resolution = 6000

        tf = target_cycles * (2 * pi)/sys%omega
        h = tf/(resolution*target_cycles)

        call solver%run(h=h, tf=tf)

        allocate(y_exact(sys%n))
        call sys%exact_solution(tf, y_exact)

        print *, "--- Integration Statistics ---"
        print '(1x,A,I6)', "Successful Steps: ", solver%steps
        print '(1x,A,I6)', "Rejected Steps:   ", solver%rejected_steps
        print '(1x,A,F12.9)', "Final Position:   ", solver%y(1)
        print '(1x,A,F12.9)', "Analytical Exact: ", y_exact
        print '(1x,A,F12.9)', "Absolute Error: ", abs(solver%y(1) - y_exact)
        print '(1x,A,F0.3,"%")', "Relative Error: ", 100*abs((solver%y(1) - y_exact) / y_exact)
    end subroutine run_harmonic_test

end program main