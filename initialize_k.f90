module initialize_k
  use szmodel
  implicit none

  real, allocatable :: M_array(:)

contains
  subroutine init_k_vec_newt(r, num_steps, max_iters)
    ! find k values by using Newton's method on a vectorised k(r_i)
    ! arguments:
    real, dimension(:), intent(in) :: r ! values to find sz_k at (not including 0)
    integer, intent(in) :: num_steps ! number of points to use in integrand
    integer, intent(in) :: max_iters ! maximum number of iterations to do

    ! local variables
    integer :: i
    real, dimension(size(r)) :: k_guess, err, err_h
    real, parameter :: h = 1e-10 ! step size for taking derivative

    M_array = sz_M(r)
    k_guess = 0.
    do i = 1, max_iters
      err = k_error_vec(k_guess, r, num_steps)
      if (maxval(abs(err)) < 1e-10) exit

      ! use backward difference to calculate derr/dk to avoid positive k values
      !   since k has a max value of (9 M(r)^2 Lambda)^(1/3) if we want a
      !   bang time
      err_h = k_error_vec(k_guess-h, r, num_steps)

      k_guess = k_guess - err*h/(err - err_h)
    end do
    k_array = [0., k_guess]
    k_0 = k_array(size(r))/(max_r**2)

    if ( i > 50 ) then
      print*,"sum error:",sum(err)
      print*,"iteration:",i
    end if
  end subroutine init_k_vec_newt

  subroutine init_k_vec_bin(r, num_steps, max_iters)
    ! find k values using binary search on a vectorised k(r_i)
    ! arguments:
    real, dimension(:), intent(in) :: r ! values to find sz_k at (not including 0)
    integer, intent(in) :: num_steps ! number of points to use in integrand
    integer, intent(in) :: max_iters ! maximum number of iterations to do

    ! local variables
    integer :: i, j
    real, dimension(size(r)) :: k_guess, k_error, k_upper, k_lower

    M_array = sz_M(r)

    k_guess = 0.
    k_upper = 0.
    k_lower = -0.25
    do i = 1, max_iters
      ! evaluate error
      k_error = k_error_vec(k_guess, r, num_steps)
      do j = 1, size(r)
        if ( k_error(j)<0 ) then
          k_upper(j) = k_guess(j)
        else
          k_lower(j) = k_guess(j)
        end if
      end do

      k_guess = 0.5*(k_lower+k_upper)

      ! if k_upper and k_lower close enough, end iteration
      if ( maxval(k_upper - k_lower) < 1e-10) then
        exit
      end if

      if (mod(i,10)==0) print*,"i:",i
      ! print*,"k_guess:",k_guess
    end do
    k_array = [0., k_guess]
  end subroutine init_k_vec_bin

  function k_error_vec(k_guess, r, num_steps)
    real, dimension(:), intent(in) :: k_guess, r
    integer, intent(in) :: num_steps
    real, dimension(size(k_guess)) :: k_error_vec

    real, dimension(size(k_guess), num_steps) :: integrand
    ! real, dimension(size(k_guess)) :: R_vals
    integer :: i
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
    k_error_vec = t_0 - sum(integrand, dim=2)*r/num_steps
  end function k_error_vec


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! scalar versions with OpenMP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function k_error_scalar(k_guess, r, M_array_index, num_steps)
    real, intent(in) :: k_guess, r
    integer, intent(in) :: M_array_index, num_steps
    real :: k_error_scalar


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
    ! find k values by using Newton's method, this method uses openMP and is
    !   faster than init_k_vec_newt on systems with multiple cores
    ! arguments:
    real, dimension(:), intent(in) :: r ! values to find sz_k at (not including 0)
    integer, intent(in) :: num_steps ! number of points to use in integrand
    integer, intent(in) :: max_iters ! maximum number of iterations to do

    ! local variables
    integer :: i, k_index
    real :: k_guess, err, err_h, r_val
    real, parameter :: h = 1e-10 ! step size for taking derivative

    M_array = sz_M(r)
    k_array(1) = 0.
    !$OMP PARALLEL DO PRIVATE(k_index, k_guess, r_val, i, err, err_h) SHARED(k_array)
    do k_index = 1, size(r)
      k_guess = 0.
      r_val = r(k_index)
      do i = 1, max_iters
        err = k_error_scalar(k_guess, r_val, k_index, num_steps)
        if (abs(err) < 1e-10) exit

        ! use backward difference to calculate derr/dk to avoid positive k values
        !   since k has a max value of (9 M(r)^2 Lambda)^(1/3) if we want a
        !   bang time
        err_h = k_error_scalar(k_guess-h, r_val, k_index, num_steps)

        k_guess = k_guess - err*h/(err - err_h)
      end do

      if ( i > 50 ) then
        print*,"error:",err
        print*,"iteration:",i
      end if

      k_array(k_index+1) = k_guess
    end do
    !$OMP END PARALLEL DO
    k_0 = k_array(size(r))/(max_r**2)
  end subroutine init_k_newt_scalar

end module initialize_k
