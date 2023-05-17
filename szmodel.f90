module szmodel
  ! functions specifying the Szekeres model and related methods

  use solveode
  use healpix_modules
  implicit none

  ! model parameters
  ! units: c=1, with distance Mpc.h^{-1}

  ! background FLRW units
  real, parameter :: H_0 = 100./(3.e5)
  real, parameter :: Omega_M = 0.315
  real, parameter :: Omega_Lambda = 1-Omega_M
  real, parameter :: Lambda = Omega_Lambda*3*H_0**2
  real, parameter :: q0 = -0.5275
  real, parameter :: j0 = 1.
  real, parameter :: t_0 = (1/(3*sqrt(Omega_Lambda)*H_0))* &
  &                       log((1+sqrt(Omega_Lambda))/(1-sqrt(Omega_Lambda)))

  ! Szekeres model parameters
  real, parameter :: M_0 = 0.5 * Omega_M * H_0**2
  real :: delta_0
  real :: alpha
  real :: r_0
  real :: Delta_r

  ! Piecewise Haantjes transformation scaling and shift parameters
  real :: scaling_factor, shift

  ! arrays to hold R(t,r) and k(r)
  real, allocatable :: k_array(:) ! holds values of k(r_i) up to max_r
  real, allocatable :: r_array(:) ! holds values of r_i
  real :: max_r ! max value in r_array
  real :: k_0 ! k(r) -> k_0 r**2 as r -> infinity
  real, allocatable :: arealR_array(:,:) ! values of R(t_i, r_i),
  !                   format: arealR_array(t_index, r_index)
  real, allocatable :: t_array(:) ! values of t_i
  real, dimension(2) :: t_domain ! max and min vals of t_array

contains

  subroutine set_params(params)
    ! set parameters to values in given vector
    real, intent(in), dimension(3) :: params
    delta_0 = params(1)
    alpha = params(2)
    r_0 = params(3)
    Delta_r = 0.1*r_0
  end subroutine set_params

  ! Mass functions
  pure function sz_M(r)
    real, intent(in), dimension(:) :: r
    real, dimension(size(r)) :: sz_M
    sz_M = M_0 * r**3 * (1+sz_delta_M(r))         ! BNW
  end function sz_M

  pure function sz_delta_M(r)
    real, intent(in), dimension(:) :: r
    real, dimension(size(r)) :: sz_delta_M
    sz_delta_M = 0.5*delta_0*(1-tanh((r-r_0)/(2*Delta_r))) ! BNW
  end function sz_delta_M

  pure function sz_dMdr(r,M)
    real, intent(in) :: r ! radial coordinate
    real, intent(in) :: M ! sz_M(r), pass as argument to minimise operations
    real :: sz_dMdr
    sz_dMdr = 3*M/r - (delta_0*M_0*r**3)/&                ! BNW
    &    (4*Delta_r*(cosh((r-r_0)/(2*Delta_r))**2))
  end function sz_dMdr



  ! dipole functions and derivatives
  subroutine get_dipoles_and_derivatives(r, S, P, Q, dS, dP, dQ, ddS, ddP, ddQ)
    ! returns the dipole functions (S,P,Q) and their first and second derivative w.r.t. r

    !! Arguments
    real, intent(in) :: r
    real, intent(out) ::  S, P, Q, dS, dP, dQ, ddS, ddP, ddQ

    !! Local Variables
    real :: r_to_the_alpha_minus_two

    ! shorthand
    r_to_the_alpha_minus_two = r**(alpha-2)

    ! S(r) and its derivatives
    S = r_to_the_alpha_minus_two*r*r
    dS = alpha*r_to_the_alpha_minus_two*r
    ddS = alpha*(alpha-1)*r_to_the_alpha_minus_two

    ! P(r) and its derivatives
    P = 0
    dP = 0
    ddP = 0

    ! Q(r) and its derivatives
    Q = 0
    dQ = 0
    ddQ = 0
  end subroutine get_dipoles_and_derivatives

  function sz_S(r)
    ! S(r) dipole function
    real, intent(in) :: r
    real :: sz_S
    
    sz_S = r**alpha
  end function sz_S

  function sz_dS(r)
    ! S'(r) dipole function first derivative
    real, intent(in) :: r
    real :: sz_dS

    sz_dS = alpha*r**(alpha-1)
  end function sz_dS

  function sz_ddS(r)
    ! second derivate of S(r)
    real, intent(in) :: r
    real :: sz_ddS

    sz_ddS = alpha*(alpha-1)*r**(alpha-2)
  end function sz_ddS

  function sz_P(r)
    ! P(r) dipole function
    real, intent(in) :: r
    real :: sz_P

    sz_P = 0.
  end function sz_P

  function sz_dP(r)
    ! P'(r) dipole function first derivative
    real, intent(in) :: r
    real :: sz_dP

    sz_dP = 0.
  end function sz_dP

  function sz_ddP(r)
    ! second derivative of P(r)
    real, intent(in) :: r
    real :: sz_ddP

    sz_ddP = 0.
  end function sz_ddP

  function sz_Q(r)
    ! Q(r) dipole function
    real, intent(in) :: r
    real :: sz_Q
    sz_Q = 0*r
  end function sz_Q

  function sz_dQ(r)
    ! Q'(r) dipole function first derivative
    real, intent(in) :: r
    real :: sz_dQ
    sz_dQ = 0*r
  end function sz_dQ

  function sz_ddQ(r)
    ! second derivative of Q(r)
    real, intent(in) :: r
    real :: sz_ddQ
    sz_ddQ = 0*r
  end function sz_ddQ

  ! anisotropy functions E and E'/E in projective coordinates
  function sz_E(r, p, q)
    real, intent(in) :: r, p, q
    real :: sz_E

    sz_E = ((p-sz_P(r))**2 + (q-sz_Q(r))**2 + sz_S(r)**2)/(2*sz_S(r))
  end function sz_E

  function sz_EdashOverE(r,p,q)
    real, intent(in) :: r,p,q
    real :: sz_EdashOverE

    real :: numer, denom

    numer = sz_dP(r)*(p-sz_P(r)) + sz_dQ(r)*(q-sz_Q(r)) - sz_dS(r)*sz_S(r)
    denom = (p-sz_P(r))**2 + (q-sz_Q(r))**2 + sz_S(r)**2

    sz_EdashOverE = -2*numer/denom - sz_dS(r)/sz_S(r)
  end function sz_EdashOverE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! R(t,r) interpolation and it's initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sz_arealR(t, r)
    !! use bilinear interpolation to find R(t,r) for any r and t value

    !! Arguments
    real, intent(in) :: t ! t value to find R at
    real, intent(in) :: r ! r value to find R at

    !! return value
    real :: sz_arealR ! value of R(t, r)

    !! local variables
    real :: r_index, r_dec_part, t_index, t_dec_part
    real :: w11, w21, w12, w22


    ! !! DEBUG
    ! if ( r<0.1 ) then
    !   print*,"sz_arealR: r=",r
    ! end if

    !! find indices in grid
    r_index = 1 + ((size(arealR_array,dim=2)-1) * r/max_r)
    r_dec_part = r_index - floor(r_index)
    t_index = 1 + ((size(arealR_array,dim=1)-1) * (t-t_domain(1))/(t_domain(2) - t_domain(1)))
    t_dec_part = t_index - floor(t_index)


    !! calculate interpolatation weights for bilinear interpolation
    w11 = (1-r_dec_part)*(1-t_dec_part)
    w21 = (1-r_dec_part)*t_dec_part
    w12 = r_dec_part*(1-t_dec_part)
    w22 = r_dec_part*t_dec_part

    !! interpolate
    sz_arealR = w11*arealR_array(floor(t_index), floor(r_index)) + &
    &           w21*arealR_array(ceiling(t_index), floor(r_index)) + &
    &           w12*arealR_array(floor(t_index), ceiling(r_index)) + &
    &           w22*arealR_array(ceiling(t_index), ceiling(r_index))
  end function sz_arealR


  pure function dRdt(arealR,r)
    !! Calculates derivative of R(t,r) with respect to t

    !! Arguments
    real, dimension(:), intent(in) :: arealR ! 1D array holding R(t,r_i) values
    real, dimension(:), intent(in) :: r ! r values where R(t,r_i) are taken.

    !! output
    real, dimension(size(r)) :: dRdt

    !! evaluate dRdt
    dRdt = sqrt(2*sz_M(r)/arealR - sz_k(r) + (Lambda*arealR**2)/3)
  end function dRdt

  function dRdr(t,r,h)
    ! evaluate R'(t,r)

    !! Arguments
    real, intent(in) :: t, r ! coords
    real, intent(in) :: h ! step size

    !! Return value
    real :: dRdr
    dRdr = (sz_arealR(t,r+h) - sz_arealR(t,r-h))/(2*h)
  end function dRdr

  subroutine find_arealR_array_omp(ti, tf, N, arealR0)
    ! solve dRdt to find R(t,r) on a grid using openMP

    !! Arguments
    real, intent(in) :: ti, tf ! initial and final t values
    integer, intent(in) :: N ! number of time steps
    real, dimension(size(r_array)-1), intent(in) :: arealR0 ! initial values of arealR on the grid r at ti excluding r=0

    !! Local variables
    integer :: r_index, t_index
    real, dimension(1) :: arealR, r, k1, k2, k3, k4
    real :: t, h
    logical :: error
    real, parameter :: ONE_SIXTH = 1./6.

    ! evaluate step size
    h = (tf-ti)/N

    !! set easy R(t,r) values
    ! set values at r=0
    arealR_array(:,1) = 0.
    ! set values at present time (t=t_0)
    arealR_array(1,2:) = arealR0

    !! set easy t values
    t_domain = [ti, tf]
    t_array(1) = ti

    ! solve ode for remaining variables
    ! ODE doesn't depend on any derivatives w.r.t. r, so paralelllize and loop
    !   first over r utilising OpenMP
    !$OMP PARALLEL DO PRIVATE(r_index,arealR,r,t,t_index,k1,k2,k3,k4) SHARED(arealR_array, t_array)
    do r_index = 1, size(arealR0)
      ! initial values
      arealR = arealR0(r_index)
      t = ti

      ! radial coordinate
      r = r_array(r_index+1)

      !! evolve one time step using rk4
      do t_index = 1, N
        ! determine slopes
        k1 = dRdt(arealR, r)
        k2 = dRdt(arealR + h*k1/2, r)
        k3 = dRdt(arealR + h*k2/2, r)
        k4 = dRdt(arealR + h*k3, r)
        ! evaluate step
        t = t + h
        arealR = arealR + ONE_SIXTH*h*(k1 + 2*k2 + 2*k3 + k4)

        ! set values into arrays
        arealR_array(t_index+1,r_index+1) = arealR(1)
        t_array(t_index+1) = t
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine find_arealR_array_omp

  subroutine find_arealR_array(ti, tf, N, arealR0)
    ! solve dRdt to find R(t,r) on a grid, using methods from solveode.f90

    !! Arguments
    real, intent(in) :: ti, tf ! initial and final t values
    integer, intent(in) :: N ! number of time steps
    real, dimension(size(r_array)), intent(in) :: arealR0 ! initial values of arealR on the grid r at ti

    ! set t_domain with range where R(t,r) is found on
    t_domain = [ti, tf]

    ! call RK4 from solveode.f90 to solve for R(t,r)
    call rk4(dRdt_wrapper, N, ti, tf, arealR0, arealR_array, t_array)

    ! set R(t,0) = 0 to remove NaNs
    arealR_array(:,1) = 0.

  contains
    pure subroutine dRdt_wrapper(y, t, darealRdt, errored)
      ! wrapper for dRdt to use it with rk4
      real, intent(in), dimension(:) :: y
      real, intent(out), dimension(size(y)):: darealRdt
      real, intent(in) :: t
      logical, intent(out) :: errored

      darealRdt = dRdt(y,r_array)+t*0
      errored = .FALSE.
    end subroutine dRdt_wrapper
  end subroutine find_arealR_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! k(r) and it's derivatives, initialization found in initialize_k.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function sz_k(r)
    ! interpolate k_array linearly to find k(r)
    ! arguments
    real, intent(in), dimension(:) :: r ! values r_i to find k(r_i) at

    ! return value
    real, dimension(size(r)) :: sz_k ! values of k(r_i)

    ! local variables
    real, dimension(size(r)) :: indices, dec_part

    ! find indices
    indices = 1 + ((size(k_array)-1) * r/max_r)
    dec_part = indices - floor(indices)
    ! interpolate
    sz_k = dec_part*k_array(ceiling(indices))+(1-dec_part)*k_array(floor(indices))
  end function sz_k

  pure function dkdr(r,h)
    ! evaluate k'(t,r) at r using sz_k and a finite difference
    real, intent(in) :: r ! coords
    real :: dkdr
    real,intent(in) :: h ! step size, must be larger than grid pixel size
    real, dimension(1) :: one_element_array

    ! evaluate finite difference
    one_element_array = (sz_k([r+h]) - sz_k([r-h]))/(2*h)
    dkdr = one_element_array(1)
  end function dkdr


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Other Useful functions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sz_F(t,r,p,q)
    ! à la B&S 2020, shorthand for g_pp and g_qq
    real, intent(in) :: t,r,p,q ! coords
    real :: sz_F

    sz_F = (sz_arealR(t,r)/sz_E(r,p,q))**2
  end function sz_F

  function sz_H(t,r,p,q)
    ! à la B&S 2020, shorthand for g_rr
    real, intent(in) :: t, r, p, q ! coords
    real :: sz_H
    real, dimension(1) :: k
    k = sz_k([r])
    sz_H = ((dRdr(t,r,r_array(2)) - sz_arealR(t,r)*sz_EdashOverE(r,p,q))**2)/(1-k(1))
  end function sz_H

  function inner_product(vector1, vector2, coords)
    ! finds inner product of two vectors, i.e. g(vector1, vector2)
    ! vector1 and vector2 must be given in projective coordinates
    real, dimension(4), intent(in) :: vector1, vector2, coords
    real :: inner_product

    inner_product = -vector1(1)*vector2(1) + &
    &               sz_H(coords(1),coords(2),coords(3),coords(4))*vector1(2)*vector2(2) +&
    &               sz_F(coords(1),coords(2),coords(3),coords(4))* &
    &               (vector1(3)*vector2(3) + vector1(4)*vector2(4))
  end function inner_product
end module szmodel
