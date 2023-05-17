! Methods for solving ordinary differential equations
module solveode

  implicit none

contains

  !! MAIN FUNCTION FOR SOLVING ODE USING RUNGE-KUTTA
  subroutine rk4(ode, N, ti, tf, yi, y_array, t_array)
    ! solve a differential equation using RK4
    ! https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    ! ode takes form
    !     dy/dt = f(y,t)
    ! and results are returned using
    ! t_array for the times,
    ! y_array for the function values

    ! arguments
    integer, intent(in) :: N ! number of steps
    real, intent(in) :: ti, tf ! initial and final times respectively
    real, intent(in), dimension(:) :: yi ! initial value of y, 1D array

    ! ode: function returning f(y,t), takes arguments y and t representing the
    !     current y value (real array) and the current time value (real scalar)
    interface
      subroutine ode(y, t, dydt, errored)
        real, intent(in), dimension(:) :: y
        real, intent(out), dimension(size(y)):: dydt
        real, intent(in) :: t
        logical, intent(out) :: errored
      end subroutine ode
    end interface

    ! output
    real, dimension(N+1, size(yi)), intent(out) :: y_array ! array of values for y
    !                                         first index: time
    !                                         second index: space
    real, dimension(N+1) :: t_array ! array of time values corresponding with
    !                             first index in y_array

    ! local variables
    real :: t_n ! current_time
    real, dimension(size(yi)) :: y_n ! current y_value
    integer :: j ! current time index for y_array
    real :: h ! step size
    logical :: error

    ! evaluate step size
    h = (tf-ti)/N

    ! set initial values
    y_n = yi
    t_n = ti
    y_array(1,:) = yi
    t_array(1) = ti

    ! evolve y and t using rk4_step
    do j=1,N
      ! print*,y_n
      call rk4_step(ode, y_n, t_n, h, error)
      y_array(j+1,:) = y_n
      t_array(j+1) = t_n
    end do
  end subroutine rk4

  subroutine rk4_step(ode, y_n, t_n, h, error)
    ! evaluate step at time t_n for rk4 to solve dy/dt = f(y,t)
    ! https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

    ! arguments:
    real, intent(inout), dimension(:)  :: y_n  ! y at step n, 1d array
    real, intent(inout)                :: t_n  ! time at step n
    real, intent(in)                :: h    ! step size
    logical, intent(out) :: error

    ! ode: subroutine returning f(y,t), takes arguments y and t representing the
    !     current y value (real array) and the current time value (real scalar),
    !     errored = TRUE if ode failed to complete
    interface
      subroutine ode(y, t, dydt, errored)
        real, intent(in), dimension(:) :: y
        real, intent(out), dimension(size(y)):: dydt
        real, intent(in) :: t
        logical, intent(out) :: errored
      end subroutine ode
    end interface


    ! local variables
    real, dimension(size(y_n)):: k1, k2, k3, k4
    logical :: err1, err2, err3, err4

    ! determine slopes
    call ode(y_n, t_n, k1, err1)
    call ode(y_n + h*k1/2, t_n + h/2, k2, err2)
    call ode(y_n + h*k2/2, t_n + h/2, k3, err3)
    call ode(y_n + h*k3, t_n + h, k4, err4)

    ! evaluate step
    t_n = t_n + h
    y_n = y_n + (1/6.)*h*(k1 + 2*k2 + 2*k3 + k4)
    
    ! return error=TRUE if any slope failed to evaluate
    error = err1 .or. err2 .or. err3 .or. err4
  end subroutine rk4_step
end module solveode
