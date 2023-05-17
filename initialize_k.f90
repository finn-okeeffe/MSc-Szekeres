! contains methods to calculate k(r_i) numerically

module initialize_k
  use szmodel
  implicit none

  real, allocatable :: M_array(:) ! M(r_i) values on the radial grid

contains
  subroutine init_k_vec_newt(r, num_steps, max_iters)
    ! find k values by using Newton's method on the vector k(r_i)

    ! arguments:
    real, dimension(:), intent(in) :: r ! values to find sz_k at (not including 0)
    integer, intent(in) :: num_steps ! number of points to use in integrand
    integer, intent(in) :: max_iters ! maximum number of iterations to do

    ! local variables
    integer :: i
    real, dimension(size(r)) :: k_guess, err, err_h ! Current guess, residual at r_i, residual at r_i+h
    real, parameter :: h = 1e-10 ! step size for taking numerical derivative

    M_array = sz_M(r) ! pre-compute M values on grid
    k_guess = 0. ! initial guess of k(r_i)

    ! Determine k(r_i) iteratively using Newton's method,
    ! breaking if the desired accuracy is reached or if
    ! we hit the max number of allowed iterations
    do i = 1, max_iters
      ! error of current guess
      err = k_error_vec(k_guess, r, num_steps)

      ! break if desired accuracy is reached
      if (maxval(abs(err)) < 1e-10) exit

      ! use backward difference to calculate derr/dk(r_i) to avoid positive k values
      !   since k has a max value of (9 M(r)^2 Lambda)^(1/3) if we want a
      !   bang time
      err_h = k_error_vec(k_guess-h, r, num_steps)

      ! update k(r_i) guess
      k_guess = k_guess - err*h/(err - err_h)
    end do

    ! set variable in szmodel.f90, manually setting k(0)=0 to avoid errors
    k_array = [0., k_guess]
    ! determine k_0 for large r FLRW approximation
    k_0 = k_array(size(r))/(max_r**2)

    ! log a warning if a large number of iterations is required
    if ( i > 50 ) then
      print*,"Warning: Large number of iterations required in init_k_vec_newt"
      print*,"sum error:",sum(err)
      print*,"iteration:",i
    end if
  end subroutine init_k_vec_newt

  subroutine init_k_vec_bin(r, num_steps, max_iters)
    ! find k values using binary search on the vector k(r_i), as per Bolejko's code

    ! arguments:
    real, dimension(:), intent(in) :: r ! values to find sz_k at (not including 0)
    integer, intent(in) :: num_steps ! number of points to use in integrand
    integer, intent(in) :: max_iters ! maximum number of iterations to do

    ! local variables
    integer :: i, j
    real, dimension(size(r)) :: k_guess, k_error, k_upper, k_lower ! current guess of k(r_i), current error, and upper and lower bounds of k(r_i)

    M_array = sz_M(r) ! set M(r_i) values

    ! initial values
    k_guess = 0.
    k_upper = 0.
    k_lower = -0.25

    ! Perform binary search
    ! iterating until desired accuracy reached, or a maximum number of iteration reached
    do i = 1, max_iters
      ! evaluate error
      k_error = k_error_vec(k_guess, r, num_steps)

      ! update k_upper or k_lower
      do j = 1, size(r)
        if ( k_error(j)<0 ) then
          k_upper(j) = k_guess(j)
        else
          k_lower(j) = k_guess(j)
        end if
      end do

      ! update guess
      k_guess = 0.5*(k_lower+k_upper)

      ! if k_upper and k_lower close enough, end iteration
      if ( maxval(k_upper - k_lower) < 1e-10) then
        exit
      end if

      ! progress logs, since this method can take a while
      if (mod(i,10)==0) print*,"i:",i
    end do
    ! set k_array variable in szmodel.f90
    k_array = [0., k_guess]
  end subroutine init_k_vec_bin

  function k_error_vec(k_guess, r, num_steps)
    ! calculate residual for a given guess of the vector k(r_i)
    ! performs integration using sum of rectangle areas

    !! arguments
    real, dimension(:), intent(in) :: k_guess, r ! guess of k(r_i) to evaluate, and r_i values
    integer, intent(in) :: num_steps ! number of rectangles to use in integration

    !! return value
    real, dimension(size(k_guess)) :: k_error_vec ! residuals of k(r_i) to return

    !! local variables
    real, dimension(size(k_guess), num_steps) :: integrand ! values of integrand, used in integration
    ! real, dimension(size(k_guess)) :: R_vals
    integer :: i

    !! calculate integrand
    forall (i=1:num_steps)
      integrand(:,i) = sqrt((r*(i-0.5)/num_steps)/(2*M_array - &
      &               k_guess*(r*(i-0.5)/num_steps) + &
      &               (1./3.)*Lambda*(r*(i-0.5)/num_steps)**3))
    end forall

    ! do i = 1, num_steps
    !   R_vals = r*(i-0.5)/num_steps
    !   integrand(:,i) = sqrt(R_vals/(2*M_array - k_guess*R_vals + (1./3.)*Lambda*R_vals**3))
    ! end do

    ! print*,"integrand at r=50, R=10",integrand(2, 80)

    !! perform integration
    k_error_vec = t_0 - sum(integrand, dim=2)*r/num_steps
  end function k_error_vec


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! scalar versions with OpenMP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function k_error_scalar(k_guess, r, M_array_index, num_steps)
    ! calculate residual of k(r_i) at a single r value
    ! performs integration using sum of rectangle areas

    ! arguments
    real, intent(in) :: k_guess, r ! k(r) guess, r value of guess
    integer, intent(in) :: M_array_index, num_steps ! index to use in M_array, number of rectangles to use in integration

    ! return value
    real :: k_error_scalar ! residual of k(r) at r

    ! local variables
    real :: R_tilde_values(num_steps), integrand(num_steps)
    integer :: i
    real, parameter :: ONE_THIRD = 1./3.


    R_tilde_values = [(r*(i-0.5)/num_steps, i=1,num_steps)] ! values of integrating variable

    ! evaluate integrand
    integrand = sqrt(R_tilde_values/(2*M_array(M_array_index) - &
    &               k_guess*R_tilde_values + &
    &               ONE_THIRD*Lambda*R_tilde_values**3))

    ! numerically integrate
    k_error_scalar = t_0 - sum(integrand)*r/num_steps
  end function k_error_scalar

  subroutine init_k_newt_scalar(r, num_steps, max_iters)
    ! find k values by using Newton's method, using a loop rather than a vector.
    ! this method uses openMP and is faster than init_k_vec_newt on uni servers

    ! arguments:
    real, dimension(:), intent(in) :: r ! values to find sz_k at (not including 0)
    integer, intent(in) :: num_steps ! number of points to use in integrand
    integer, intent(in) :: max_iters ! maximum number of iterations to do

    ! local variables
    integer :: i, k_index
    real :: k_guess, err, err_h, r_val ! Current guess, residual at r_i, residual at r_i+h, current r value
    real, parameter :: h = 1e-10 ! step size for taking derivative

    ! calculate M_array
    M_array = sz_M(r)

    ! set k(0) = 0
    k_array(1) = 0.

    !! loop over r grid using OpenMP, determining k(r_i) at each r_i value
    !$OMP PARALLEL DO PRIVATE(k_index, k_guess, r_val, i, err, err_h) SHARED(k_array)
    do k_index = 1, size(r)
      ! initial guess
      k_guess = 0.

      ! set r value
      r_val = r(k_index)

      ! perform Newton's method to find k(r_i)
      do i = 1, max_iters
        err = k_error_scalar(k_guess, r_val, k_index, num_steps)
        if (abs(err) < 1e-10) exit

        ! use backward difference to calculate derr/dk to avoid positive k values
        !   since k has a max value of (9 M(r)^2 Lambda)^(1/3) if we want a
        !   bang time
        err_h = k_error_scalar(k_guess-h, r_val, k_index, num_steps)

        k_guess = k_guess - err*h/(err - err_h)
      end do

      ! print a warning if a large number of iterations is required for the current r_i value
      if ( i > 50 ) then
        print*,"Warning: Large number of iterations required in init_k_newt_scalar"
        print*,"error:",err
        print*,"iteration:",i
      end if

      ! set k_array variable in szmodel.f90
      k_array(k_index+1) = k_guess
    end do
    !$OMP END PARALLEL DO

    ! set k_0 for large r FLRW approximation
    k_0 = k_array(size(r))/(max_r**2)
  end subroutine init_k_newt_scalar
end module initialize_k
