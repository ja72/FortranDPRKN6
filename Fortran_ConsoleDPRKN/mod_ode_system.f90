module mod_ode_system
    use mod_common
    implicit none
    private
    
    public :: ode_system, ode_solver, init_ode_solver
    public :: i_calc_ddy, i_exact_solution, i_ode_update
    
    type, abstract :: ode_system
        integer :: n
        real(wp), allocatable :: y_init(:)
        real(wp), allocatable :: dy_init(:)
        procedure(i_exact_solution), pointer :: exact_solution => NULL()
    contains
        procedure(i_ode_update), deferred :: ode_update
        procedure(i_calc_ddy), deferred :: calc_ddy
    end type
    
    type :: ode_solver(n)
        integer , len:: n
        real(wp) :: t
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: dy(:)
        integer :: steps, rejected_steps
        real(wp) :: reltol
        real(wp) :: abstol
        class(ode_system), pointer :: ode_sys => NULL()
    contains 
        procedure :: reset => reset_solver
        procedure :: step => step_solver
        procedure :: run => run_solver
    end type

    interface
        pure subroutine i_calc_ddy(self, t, y, dy, ddy)
        import :: wp, ode_system
        class(ode_system), intent(in) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in), dimension(self%n) :: y, dy
        real(wp), intent(out), dimension(self%n) :: ddy
        end subroutine i_calc_ddy
        
        pure subroutine i_exact_solution(self, t, y_exact)
        import :: wp, ode_system
        class(ode_system), intent(in) :: self
        real(wp), intent(in) :: t
        real(wp), intent(out), dimension(self%n) :: y_exact
        end subroutine i_exact_solution
        
        subroutine i_ode_update(self, solver, h)
        import :: wp, ode_system, ode_solver
        class(ode_system), intent(in) :: self
        class(ode_solver(self%n)), intent(in) :: solver
        real(wp), intent(in) :: h
        end subroutine i_ode_update
    end interface
    
    interface ode_solver
        module procedure :: init_ode_solver
    end interface

contains
    
    function init_ode_solver(ode_sys) result(solver)
    !integer, intent(in) :: n
    class(ode_system), intent(in), pointer :: ode_sys
    type(ode_solver(ode_sys%n)), pointer :: solver
    integer :: n
        n = ode_sys%n
        
        if( .not. associated(ode_sys)) then
            error stop "Error: ode_system object is not associated."
        end if
        allocate(solver)
        solver%ode_sys => ode_sys
        solver%reltol = 1.0e-8_wp
        solver%abstol = 1.0e-8_wp
        allocate(solver%y(n))
        allocate(solver%dy(n))
        call reset_solver(solver)
    end function init_ode_solver
    
    pure subroutine reset_solver(self)
    class(ode_solver(*)), intent(inout) :: self
    integer :: n
        n = self%ode_sys%n
        self%t = 0.0_wp
        self%steps = 0
        self%rejected_steps = 0
        self%y(1:n) = self%ode_sys%y_init(1:n)
        self%dy(1:n) = self%ode_sys%dy_init(1:n)
    end subroutine reset_solver
    
    subroutine run_solver(self, h, tf)
    class(ode_solver(*)), intent(inout) :: self
    real(wp), intent(inout) :: h
    real(wp), intent(in) :: tf
    class(ode_system), pointer :: ode_sys
    
    integer :: n
    real(wp) :: h_next, scale
    real(wp) :: y_next(self%n), dy_next(self%n)
    real(wp) :: err_y(self%n), err_dy(self%n), total_err
    
    if( .not. associated(self%ode_sys) ) then
        error stop "Error: ode_solver object has no associated ode_system."
    end if
    
    ode_sys => self%ode_sys

    n = self%n
    self%steps = 0
    self%rejected_steps = 0

    call ode_sys%ode_update(self, h)

    do while (self%t < tf)
        ! Prevent overshooting the bounds of tf
        if (self%t + h > tf) h = tf - self%t
    
        ! Take trial step
        !call dprkn6_step(harmonic_oscillator, n, self%t, h, y, dy, y_next, dy_next, err_y, err_dy)
        call self%step(h, y_next, dy_next, err_y, err_dy)
    
        ! Calculate scalar error norm (Max normalized metric)
        total_err = maxval(abs(err_y) / (self%abstol + self%reltol * abs(y_next)))
        total_err = max(total_err, maxval(abs(err_dy) / (self%abstol + self%reltol * abs(dy_next))))
    
        ! Determine step acceptance based on local error safety limits
        if (total_err <= 1.0_wp) then
            ! Step accepted! Advance states
            self%t = self%t + h
            self%y = y_next
            self%dy = dy_next
            self%steps = self%steps + 1
        
            ! Optional: print milestones
            call ode_sys%ode_update(self, h)
            
            ! Compute next adaptive step size (safety factored)
            scale = 0.9_wp * (1.0_wp / total_err)**(1.0_wp / 6.0_wp)
            h_next = h * max(0.2_wp, min(5.0_wp, scale))
            h = h_next
        else
            ! Step Rejected! Shrink h and re-evaluate without changing t
            self%rejected_steps = self%rejected_steps + 1
            scale = 0.9_wp * (1.0_wp / total_err)**(1.0_wp / 5.0_wp)
            h = h * max(0.1_wp, min(0.9_wp, scale))
        end if
    end do    
    end subroutine
    
    pure subroutine step_solver(self, h, y_next, dy_next, err_y, err_dy)
    class(ode_solver(*)), intent(inout) :: self
    real(wp), intent(in) :: h
    real(wp), intent(out) :: y_next(self%n), dy_next(self%n)
    real(wp), intent(out) :: err_y(self%n), err_dy(self%n)
    
        call dprkn6_step(f, self%n, self%t, h, self%y, self%dy, y_next, dy_next, err_y, err_dy)
    
    contains
    
    pure subroutine f(t, n, y, dy, ddy)
    real(wp), intent(in) :: t
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n), dy(n)
    real(wp), intent(out) :: ddy(n)
        call self%ode_sys%calc_ddy(t,y,dy,ddy)            
    end subroutine f
        
    pure subroutine dprkn6_step(f, n, t, h, y, dy, y_next, dy_next, err_y, err_dy)
    ! Dormand & Prince 1987 RKN(6,4,6) Coefficients 
        integer, parameter :: STAGES = 6
    
        real(wp), parameter :: c(STAGES) = [ &
            0.0_wp, 0.0870311222683072_wp, 0.2227182283087289_wp, &
            0.5103445105261314_wp, 0.8118021966205844_wp, 1.0_wp]
        
        real(wp), parameter :: a(STAGES, STAGES) = reshape([ &
            0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp, &
            0.0037872061902041_wp,0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp, &
            0.0094056263590518_wp,0.0154131557026779_wp,0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp, &
            0.0381643445203362_wp,0.0016641849646879_wp,0.0903332468351608_wp,0.0_wp,               0.0_wp,               0.0_wp, &
            0.1118128522770744_wp,-0.1264875619727402_wp,0.2798939912071854_wp,0.0643264629472304_wp,0.0_wp,               0.0_wp, &
            0.1172837260538965_wp,0.0461825129618174_wp,0.1783510580975610_wp,0.1189311654388318_wp,0.0392515374478933_wp,0.0_wp  &
        ], [STAGES, STAGES])

        real(wp), parameter :: a_prime(STAGES, STAGES) = reshape([ &
            0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp, &
            0.0870311222683072_wp,0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp, &
            0.0526019342795898_wp,0.1701162940291391_wp,0.0_wp,               0.0_wp,               0.0_wp,               0.0_wp, &
            0.1039860667357774_wp,-0.1471018599424707_wp,0.5534603037328247_wp,0.0_wp,               0.0_wp,               0.0_wp, &
            0.2736183602568603_wp,-0.5731997327318057_wp,0.8524451044437175_wp,0.2589384646518123_wp,0.0_wp,               0.0_wp, &
            0.1147544074211158_wp,0.0891040854495914_wp,0.2286704179354032_wp,0.4627443831776992_wp,0.1047267060161904_wp,0.0_wp  &
        ], [STAGES, STAGES])

        ! 6th-order integration vectors
        real(wp), parameter :: b(STAGES)       = [0.1172837260538965_wp, 0.0461825129618174_wp, 0.1783510580975610_wp, 0.1189311654388318_wp, 0.0392515374478933_wp, 0.0_wp]
        real(wp), parameter :: b_prime(STAGES) = [0.1147544074211158_wp, 0.0891040854495914_wp, 0.2286704179354032_wp, 0.4627443831776992_wp, 0.1047267060161904_wp, 0.0_wp]
    
        ! 4th-order embedding weights for adaptive error estimation (E = O(h^5) position error, E' = O(h^5) velocity error)
        real(wp), parameter :: e(STAGES)       = [0.1221711200155609_wp, 0.0249767280803513_wp, 0.1764619729909204_wp, 0.1384074360341775_wp, 0.0326462725514081_wp, 0.0053364703275818_wp]
        real(wp), parameter :: e_prime(STAGES) = [0.1039860667357774_wp, 0.1118671569424855_wp, 0.2104595292415174_wp, 0.4578132009893992_wp, 0.0833131707567084_wp, 0.0325608753341121_wp]
        interface
            pure subroutine f(t, n, y, dy, ddy)
                import :: wp
                integer, intent(in) :: n
                real(wp), intent(in) :: t
                real(wp), intent(in) :: y(n), dy(n)
                real(wp), intent(out) :: ddy(n)
            end subroutine f
        end interface
        
        integer, intent(in)          :: n
        real(wp), intent(in)         :: t, h
        real(wp), intent(in)         :: y(n), dy(n)
        real(wp), intent(out)        :: y_next(n), dy_next(n)
        real(wp), intent(out)        :: err_y(n), err_dy(n)
        
        real(wp) :: k(n, STAGES)
        real(wp) :: y_stage(n), dy_stage(n)
        real(wp) :: t_stage
        integer  :: s, j
        
        ! Stage 1 evaluation
        call f(t, n, y, dy, k(:, 1))
        
        ! Stages 2 through 6
        do s = 2, STAGES
            t_stage = t + c(s) * h
            y_stage = y + c(s) * h * dy
            dy_stage = dy
            
            do j = 1, s - 1
                y_stage  = y_stage  + (h**2) * a(j, s) * k(:, j)
                dy_stage = dy_stage + h * a_prime(j, s) * k(:, j)
            end do
            
            call f(t_stage, n, y_stage, dy_stage, k(:, s))
        end do
        
        ! Compute Final 6th-order accurate update
        y_next  = y  + h * dy
        dy_next = dy
        do s = 1, STAGES
            y_next  = y_next  + (h**2) * b(s) * k(:, s)
            dy_next = dy_next + h * b_prime(s) * k(:, s)
        end do
        
        ! Compute local truncation error using embedded weights
        err_y  = (b - e)       * (h**2)
        err_dy = (b_prime - e_prime) * h
        ! Project error vectors relative to stages evaluated
        ! (Actual raw vector diff is computed against the aggregated stage results)
        err_y  = 0.0_wp
        err_dy = 0.0_wp
        do s = 1, STAGES
            err_y  = err_y  + (h**2) * (b(s) - e(s)) * k(:, s)
            err_dy = err_dy + h * (b_prime(s) - e_prime(s)) * k(:, s)
        end do
    end subroutine dprkn6_step
    
    end subroutine step_solver
    
end module