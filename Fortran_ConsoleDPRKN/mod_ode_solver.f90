module mod_ode_solver
    use mod_common
    implicit none
    private

    public :: calc_ddy, ode_solver
    
    type :: ode_solver
        integer :: n
        integer :: steps, rejected_steps
        real(wp) :: t
        real(wp), allocatable :: y(:), dy(:)
        real(wp), allocatable :: y_init(:), dy_init(:)
        real(wp) :: reltol, abstol
        procedure(calc_ddy), pointer :: f => null()
    contains
        procedure :: reset => ode_reset
        procedure :: step => ode_step
        procedure :: run => ode_simulate
    end type ode_solver
    
    interface ode_solver
        module procedure :: ode_init
    end interface
    
    interface 
        pure subroutine calc_ddy(obj, t, y, dy, ddy)
        import
        class(ode_solver), intent(in) :: obj
        real(wp), intent(in) :: t
        real(wp), intent(in), dimension(obj%n) :: y, dy
        real(wp), intent(out), dimension(obj%n) :: ddy
        end subroutine calc_ddy
    end interface
    
    contains
    
    pure function ode_init(f_ptr, n) result(obj)
    type(ode_solver) :: obj
    procedure(calc_ddy), intent(in), pointer :: f_ptr
    integer, intent(in) :: n
        obj%f => f_ptr
        obj%reltol = 1.0e-8_wp
        obj%abstol = 1.0e-8_wp
        obj%n = n
    
        call ode_reset(obj, t=0.0_wp)
    end function ode_init
    
    pure subroutine ode_reset(obj, t, y_init, dy_init)
    class(ode_solver), intent(inout) :: obj
    real(wp), intent(in), optional :: t
    real(wp), intent(in), optional :: y_init(obj%n), dy_init(obj%n)
    
        if( present(t) ) then
            obj%t = t
        end if
    
        if( present(y_init) ) then
            obj%y_init = y_init
        end if
        if( present(dy_init) ) then
            obj%dy_init = dy_init
        end if
    
        obj%y  = obj%y_init
        obj%dy = obj%dy_init
        obJ%steps = 0
        obj%rejected_steps = 0
    end subroutine ode_reset
    
    pure subroutine ode_step(obj, h, y_next, dy_next, err_y, err_dy)
    !use mod_dprkn6
    class(ode_solver), intent(inout) :: obj
    real(wp), intent(in) :: h
    real(wp), intent(out) :: y_next(obj%n), dy_next(obj%n)
    real(wp), intent(out) :: err_y(obj%n), err_dy(obj%n)
    
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
    

    call dprkn6_step(f, obj%n, obj%t, h, obj%y, obj%dy, y_next, dy_next, err_y, err_dy)
    
    contains
        pure subroutine f(t, n, y, dy, ddy)
        real(wp), intent(in) :: t
        integer, intent(in) :: n
        real(wp), intent(in) :: y(n), dy(n)
        real(wp), intent(out) :: ddy(n)
            call obj%f(t,y,dy,ddy)            
        end subroutine f
        
        pure subroutine dprkn6_step(f, n, t, h, y, dy, y_next, dy_next, err_y, err_dy)
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
        
    end subroutine ode_step
    
    subroutine ode_simulate(obj, h, tf, update)
    class(ode_solver), intent(inout) :: obj
    real(wp), intent(inout) :: h
    real(wp), intent(in) :: tf
    
    interface
        subroutine ode_update(obj, h)
        import
        class(ode_solver), intent(inout) :: obj
        real(wp), intent(in) :: h
        end subroutine
    end interface
    procedure(ode_update), optional :: update
    
    integer :: n
    real(wp) :: h_next, scale
    real(wp) :: y_next(obj%n), dy_next(obj%n)
    real(wp) :: err_y(obj%n), err_dy(obj%n), total_err
    logical :: has_update
    
    if( present(update) ) then
        has_update = .true.
    else
        has_update = .false.
    end if
    
    n = obj%n
    obj%steps = 0
    obj%rejected_steps = 0
    
    ! Optional: print milestones
    if (has_update)  then 
        call update(obj, h)
    end if
    
    do while (obj%t < tf)
        ! Prevent overshooting the bounds of tf
        if (obj%t + h > tf) h = tf - obj%t
        
        ! Take trial step
        !call dprkn6_step(harmonic_oscillator, n, obj%t, h, y, dy, y_next, dy_next, err_y, err_dy)
        call obj%step(h, y_next, dy_next, err_y, err_dy)
        
        ! Calculate scalar error norm (Max normalized metric)
        total_err = maxval(abs(err_y) / (obj%abstol + obj%reltol * abs(y_next)))
        total_err = max(total_err, maxval(abs(err_dy) / (obj%abstol + obj%reltol * abs(dy_next))))
        
        ! Determine step acceptance based on local error safety limits
        if (total_err <= 1.0_wp) then
            ! Step accepted! Advance states
            obj%t = obj%t + h
            obj%y = y_next
            obj%dy = dy_next
            obj%steps = obj%steps + 1
            
            ! Optional: print milestones
            if (has_update)  then 
                call update(obj, h)
            end if
            ! Compute next adaptive step size (safety factored)
            scale = 0.9_wp * (1.0_wp / total_err)**(1.0_wp / 6.0_wp)
            h_next = h * max(0.2_wp, min(5.0_wp, scale))
            h = h_next
        else
            ! Step Rejected! Shrink h and re-evaluate without changing t
            obj%rejected_steps = obj%rejected_steps + 1
            scale = 0.9_wp * (1.0_wp / total_err)**(1.0_wp / 5.0_wp)
            h = h * max(0.1_wp, min(0.9_wp, scale))
        end if
    end do
    

    end subroutine ode_simulate
    
end module