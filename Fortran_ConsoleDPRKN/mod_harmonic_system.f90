module mod_harmonic_system
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use mod_ode_system
    implicit none
    private
    
    public :: harmonic_system
    
    type, extends(ode_system) :: harmonic_system
        real(wp) :: omega
    contains
        procedure :: calc_ddy => harmonic_calc_ddy
        procedure :: ode_update => harmonic_update
    end type
    
    interface harmonic_system
        module procedure :: init_harmonic_system
    end interface

    contains
    
    pure function init_harmonic_system(omega, y_init, dy_init) result(sys)
    real(wp), intent(in) :: omega
    real(wp), intent(in), dimension(:) :: y_init, dy_init
    type(harmonic_system) :: sys

        sys%n = size(y_init)
        allocate(sys%y_init(sys%n))
        allocate(sys%dy_init(sys%n))
        sys%y_init = y_init
        sys%dy_init = dy_init
        sys%omega = omega
        sys%exact_solution => harmonic_exact_solution
    end function init_harmonic_system

    pure subroutine harmonic_calc_ddy(self, t, y, dy, ddy)
        class(harmonic_system), intent(in) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in), dimension(self%n) :: y, dy
        real(wp), intent(out), dimension(self%n) :: ddy
        ddy(1) = -self%omega**2 * y(1)
    end subroutine harmonic_calc_ddy

    pure subroutine harmonic_exact_solution(self, t, y_exact)
        class(ode_system), intent(in) :: self
        real(wp), intent(in) :: t
        real(wp), intent(out), dimension(self%n) :: y_exact
        
        select type (self)
        type is (harmonic_system)
            y_exact(1) = self%y_init(1) * cos(self%omega * t) + (self%dy_init(1) / self%omega) * sin(self%omega * t)
        class default 
            error stop "Error: harmonic_exact_solution called with an object that is not of type harmonic_system."
        end select
    end subroutine harmonic_exact_solution
    
    subroutine harmonic_update(self, solver, h)
    ! This subroutine can be used to update any parameters of the harmonic system if needed.
    ! For now, it does nothing as the harmonic oscillator parameters are constant.        
    class(harmonic_system), intent(in) :: self
    real(wp), intent(in) :: h
    class(ode_solver(self%n)), intent(in) :: solver
    
    integer, parameter :: poll_steps = 40
    
        if( solver%steps == 0) then
            write(*,997) "Step", "Time", "Position", "Speed", "Step Size"
        end if
    
        if( mod(solver%steps, poll_steps) == 0)  then
            write(*,998) solver%steps, solver%t, solver%y(1), solver%dy(1), h
        end if
        
997     format (1x,A6,1x,A6,1x,A10,1x,A10,1x,A10)
998     format (1x,I6,1x,F6.2,1x,F10.6,1x,F10.6,1x,F10.6)
    
    end subroutine harmonic_update

end module mod_harmonic_system