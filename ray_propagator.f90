module ray_propagator


  ! methods and types for ray tracing
  use szmodel
  implicit none

  real, parameter :: max_precision = 100
  real, parameter :: min_theta = pi/2
  integer, parameter :: test_steps = 20001
  real, parameter :: r_stop = 400 ! point at which geometry ~ FLRW

  type ray
    ! represents a ray, holds relevant data
    real(8), dimension(2) :: sky_angles ! position on sky of initial ray
    !                               [theta, phi]
    real, dimension(10) :: state_vector ! holds current state
    ! state_vector = [coords, 4 momentum, d_A, dd_A/ds]
    real :: affine_parameter ! current affine parameter value
    real :: d_L
    logical :: inverted = .FALSE. ! whether photon stats are inverted from
    !                             the original coordinates
    logical :: to_be_inverted = .FALSE.


    integer :: times_inverted

    logical :: errored ! if true, ray must be discarded. Set to true when:
    !                         - r<0 at any point in ray tracing process.
    character(len = 20) :: err_message

    real, allocatable :: stv_history(:,:)
    real, allocatable :: aff_param_hist(:)
    logical, allocatable :: inv_hist(:)
  end type ray


contains
  subroutine propagate_ray(photon, d_L, lam_step_size, test)
    ! evolves a ray back in time from the observer to distance d_L from the observer
    type(ray), intent(inout) :: photon ! ray to propagate
    real, intent(in) :: d_L ! luminosity distance to max out at
                            ! if negative, will run until reaching r_stop
    real, intent(in) :: lam_step_size ! size of step in affine parameter
    logical, intent(in) :: test

    logical :: errored
    integer :: i
    real :: step, S, P, Q, norm

    !! TESTING
    real, dimension(10) :: old_sv

    i = 2
    old_sv = photon%state_vector
    if ( test ) then
      photon%stv_history(1,:) = photon%state_vector
      photon%aff_param_hist(1) = 0
    end if

      ! normal code
      ! print*,"test diff_state_vector", diff_state_vector(photon%state_vector, 0.)
      ! solve using rk4 in Szekeres up until photon reaches desired d_L or r_max
    do
      ! rk4 update step
      step = lam_step_size
      ! if ( photon%state_vector(2) < 1 ) then
      !   step = step/100
      ! end if
      call rk4_step(diff_state_vector, photon%state_vector, &
      &             photon%affine_parameter, step, errored)
      ! set luminosity distance using reciprocity theorem
      photon%d_L = photon%state_vector(9)*(photon%state_vector(5)**2)

      ! S = sz_S(photon%state_vector(2))
      ! P = sz_P(photon%state_vector(2))
      ! Q = sz_Q(photon%state_vector(2))
      ! if ( photon%inverted ) then
      !   norm = P**2 + Q**2 + S**2
      !   S = S/norm
      !   P = P/norm
      !   Q = Q/norm
      ! end if



      old_sv = photon%state_vector

      !! check if p and q too large, if so perform inversion
      ! if ( max(abs(photon%state_vector(3)),abs(photon%state_vector(4))) > max_precision) then
      !   call invert_ray(photon)
      ! end if
      ! invert when (p,q) gets below theta=min_theta rather than at a fixed (p,q),
      !   that seems to avoid rapid changes in (p,q) around blow up,
      !   keeping the error low
      ! if (sum((photon%state_vector(3:4)-[P,Q])**2) > (S/tan(min_theta/2))**2) then
      !   call invert_ray(photon)
      ! end if
      ! if (sum(photon%state_vector(3:4)**2) > 1/tan(min_theta/2)**2) then
      !   call invert_ray(photon)
      ! end if
      if ( photon%to_be_inverted ) then
        call invert_ray(photon)
        photon%to_be_inverted = .FALSE.
      end if


      ! record geodesic history if test=TRUE
      if ( test .and. i .le. test_steps) then
        photon%stv_history(i,:) = photon%state_vector
        photon%inv_hist(i) = photon%inverted
        photon%aff_param_hist(i) = photon%affine_parameter
      end if

      ! if errored, stop propagation
      if ( errored ) then
        photon%errored = .TRUE.
        exit
      end if

      ! exit if passed target d_L (tracing COMPOSITE)
      if (d_L > 0 .and. ((lam_step_size>0 .and. photon%d_L < d_L) .or. &
      &       (lam_step_size<0 .and. photon%d_L > d_L))) then
        exit
      end if

      ! exit when tracing CMB if reached FLRW region
      if ( d_L < 0 .and.  photon%state_vector(2) > r_stop) then
        exit
      end if

      i = i+1
    end do


  contains
    subroutine diff_state_vector(s_vec,aff_param, dsdt, error)
      ! derivative of the state vector
      ! arguments and return value
      real, intent(in), dimension(:) :: s_vec ! state vector of photon
      real, intent(in) :: aff_param ! affine parameter value
      real, intent(out), dimension(size(s_vec)):: dsdt ! derivative for rk4_step
      logical, intent(out) :: error

      ! local variables
      real, dimension(4) :: diff_4p ! derivative of the 4-momentum
      real :: diff_d_A ! derivative of angular distance w.r.t. affine parameter

      ! variables for coding in geodesic equations
      real :: kpq2, kr2, F, H
      real, dimension(4) :: dF, dH
      real :: arealR, darealRdt, darealRdr, darealRdrdt, Rdashdash
      real :: r_step, Rplush, Rminush
      real :: S, P, Q, dS, dP, dQ, ddS, ddP, ddQ
      real :: E, EdashOverE, k, one_over_one_minus_k, dk, a1, pdiff, qdiff, a2
      real :: EEdr, EEdp, EEdq, EEpSS
      real :: M, dM
      real, dimension(1) :: one_element_array

      real :: norm, dnorm, ddnorm
      real :: iS, iP, iQ, idS, idP, idQ, iddS, iddP, iddQ

      r_step=r_array(2)

      ! if NaN's appear
      if (any(.not. s_vec == s_vec)) then
              print*,"NaN in s_vec:",s_vec
              print*,"sky angles:",photon%sky_angles
      end if


      !! check for errored
      if ( s_vec(2) < r_step ) then
        ! print*,"too small r photon",s_vec(2)
        dsdt = 0.
        ! photon%errored = .TRUE.
        photon%err_message = "radius too small"
        error = .TRUE.
      else
        !! calculate derivative
        ! values needed

        kpq2 = s_vec(7)**2 + s_vec(8)**2                           ! kp**2 + kq**2
        kr2 = s_vec(6)**2                                                  ! kr**2

        if ( s_vec(2) < r_stop ) then ! within inhomogeneity
          arealR = sz_arealR(s_vec(1),s_vec(2))                             ! R(t,r)

          Rplush = sz_arealR(s_vec(1),s_vec(2)+r_step)
          Rminush = sz_arealR(s_vec(1),s_vec(2)-r_step)
          darealRdr = (Rplush-Rminush)/(2*r_step)                          ! R'(t,r)
          Rdashdash = (Rplush - 2*arealR + Rminush)/(r_step**2)           ! R''(t,r)

          one_element_array = sz_k([s_vec(2)])                              ! k(r)
          k = one_element_array(1)
          dk = dkdr(s_vec(2), r_step)                                       !k'(r)
        else ! in FLRW region
          arealR = s_vec(2)*sz_arealR(s_vec(1),max_r)/max_r
          darealRdr = sz_arealR(s_vec(1),max_r)/max_r
          Rdashdash = 0.
          k = k_0*s_vec(2)**2
          dk = 2*k_0*s_vec(2)
        end if
        one_over_one_minus_k = 1/(1-k)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! DIPOLE FUNCTIONS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_dipoles_and_derivatives(s_vec(2), S,P,Q,dS,dP,dQ,ddS,ddP,ddQ)

        ! adjust dipole function and derivatives if photon is inverted
        if ( photon%inverted ) then
          norm = P**2 + Q**2 + S**2
          dnorm = 2*(P*dP + Q*dQ + S*dS)
          ddnorm = 2*(dP**2 + P*ddP + dQ**2 + Q*ddQ + dS**2 + S*ddS)

          iP = P/norm
          iQ = Q/norm
          iS = S/norm

          idP = (dP*norm - P*dnorm)/(norm**2)
          idQ = (dQ*norm - Q*dnorm)/(norm**2)
          idS = (dS*norm - S*dnorm)/(norm**2)

          iddP = (ddP*norm - 2*dP*dnorm - P*ddnorm + 2*P*dnorm*dnorm/norm)/(norm**2)
          iddQ = (ddQ*norm - 2*dQ*dnorm - Q*ddnorm + 2*Q*dnorm*dnorm/norm)/(norm**2)
          iddS = (ddS*norm - 2*dS*dnorm - S*ddnorm + 2*S*dnorm*dnorm/norm)/(norm**2)

          P = iP
          Q = iQ
          S = iS
          dP = idP
          dQ = idQ
          dS = idS
          ddP = iddP
          ddQ = iddQ
          ddS = iddS
        end if


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! FUNCTIONS APPEARING IN CONNECTION
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pdiff = (s_vec(3)-P)/S                                          !(p-P)/S
        qdiff = (s_vec(4)-Q)/S                                          !(q-Q)/S

        E = 0.5*S*(pdiff**2 + qdiff**2 + 1)                           ! E(r,p,q)
        EdashOverE = -(dP*pdiff+dQ*qdiff-dS)/E - dS/S                     ! E'/E

        a2 = darealRdr - arealR*EdashOverE                             !R'-RE'/E

        F = (arealR/E)**2                                                     !F
        H = (a2**2)*one_over_one_minus_k                                      !H

        a1 = 2*arealR/(E**2)                                            ! 2R/E^2

        one_element_array = sz_M([s_vec(2)])                               !M(r)
        M = one_element_array(1)

        dM = sz_dMdr(s_vec(2),M)                                         ! M'(r)

        darealRdt= sqrt(2*M/arealR - k + (Lambda*arealR**2)/3)            ! dR/dt
        darealRdrdt = (dM/arealR - M*darealRdr/(arealR**2) - 0.5*dk + & !d^2R/drdt
        &               Lambda*arealR*darealRdr/3)/darealRdt

        EEpSS = EdashOverE + dS/S


        ! Derivatives of E'/E
        EEdr = -(ddP*pdiff+ddQ*qdiff-ddS - (dQ**2+dP**2+dS**2)/S)/E &
        &       -EEpSS**2 - ddS/S + (dS/S)**2
        EEdp = -dP/(E*S) - pdiff*EEpSS/E
        EEdq = -dQ/(E*S) - qdiff*EEpSS/E

        !! DERIVATIVES OF F
        dF(1) = a1*darealRdt                                               ! dF/dt
        dF(2) = a1*a2                                                      ! dF/dr
        dF(3) = -a1**2*E*pdiff/2                                           ! dF/dp
        dF(4) = -a1**2*E*qdiff/2                                           ! dF/dp

        !! DERIVATIVES OF H

        dH(1) = 2*a2*(darealRdrdt - darealRdt*EdashOverE)*one_over_one_minus_k !dH/dt
        dH(2) = (2*a2*(Rdashdash-darealRdr*EdashOverE-arealR*EEdr) &       ! dH/dr
        &         + H*dk)*one_over_one_minus_k
        dH(3) = -2*a2*arealR*EEdp*one_over_one_minus_k                     ! dH/dp
        dH(4) = -2*a2*arealR*EEdq*one_over_one_minus_k                     ! dH/dq


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! DERIVATIVES OF STATE VECTOR DERIVATIVE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        diff_4p = 0.
        diff_4p(1) = -0.5*(dH(1)*kr2 + dF(1)*kpq2)
        diff_4p(2) = (0.5*(dH(2)*kr2 + dF(2)*kpq2) - s_vec(6)*sum(dH*s_vec(5:8)))/H
        diff_4p(3) = (0.5*(dH(3)*kr2 + dF(3)*kpq2) - s_vec(7)*sum(dF*s_vec(5:8)))/F
        diff_4p(4) = (0.5*(dH(4)*kr2 + dF(4)*kpq2) - s_vec(8)*sum(dF*s_vec(5:8)))/F

        
        diff_d_A = -0.5*s_vec(9)*(dM-3*M*EdashOverE)*s_vec(5)**2 /(a2*arealR**2)
        ! set derivative vector
        dsdt = [s_vec(5:8), diff_4p, s_vec(10),diff_d_A]
        error = .FALSE.

        if (any(.not. dsdt == dsdt)) then
                print*,"NaN in dsdt:",dsdt
        end if

      
      ! check if photon needs to be inverted
      if (sum((s_vec(3:4)-[P,Q])**2) > (S/tan(min_theta/2))**2) then
        photon%to_be_inverted = .TRUE.
      end if



        ! if ( diff_4p(2) > 1.e2 ) then
        !   print*,"Large dk^r:", diff_4p(2)
        !   print*,"  coords:", s_vec(1:4)
        !   print*,"  k:", s_vec(5:8)
        !   print*,"  (F,H):",F,H
        !   print*,"  dF:",dF
        !   print*,"  dH:",dH
        !   print*,"  (S,P,Q):",S,P,Q
        !   print*,"  (dS,dP,dQ):",dS,dP,dQ
        !   print*,"  (ddS,ddP,ddQ):",ddS,ddP,ddQ
        ! end if
        ! if (photon%inverted) then
        !   print*," inverted"
        ! else
        !   print*,"  uninverted"
        ! end if

      end if



    end subroutine diff_state_vector
  end subroutine propagate_ray


  function new_ray(sky_angles, coords, record_history)
    ! initializes a new ray at coords pointing at direction sky_angles in the sky
    real(8), dimension(2), intent(in) :: sky_angles ! direction in sky that ray comes from
    !                                 [theta, phi]
    real, dimension(4), intent(in) :: coords ! coordinates of the observer
    logical, intent(in) :: record_history
    type(ray) :: new_ray ! ray object to return

    ! local variables
    real, dimension(4) :: four_momentum
    real, dimension(3) :: local_3p ! spatial components of momentum in the local frame

    !! initialise the easy stuff
    new_ray%sky_angles = sky_angles
    new_ray%affine_parameter = 0.
    new_ray%d_L = 0.
    new_ray%times_inverted = 0
    new_ray%errored = .FALSE.
    new_ray%to_be_inverted = .FALSE.

    if ( record_history ) then
      print*,"Allocating History arrays with",test_steps,"steps"
      allocate(new_ray%stv_history(test_steps, 10))
      allocate(new_ray%aff_param_hist(test_steps))
      allocate(new_ray%inv_hist(test_steps))
      print*,"rp --> stv_history shape:",shape(new_ray%stv_history)
      print*,"rp --> aff_param_hist shape:",shape(new_ray%aff_param_hist)
      print*,"rp --> inv_hist shape:",shape(new_ray%inv_hist)
      print*,"rp --> inv_hist allocated:",allocated(new_ray%inv_hist)
    end if

    !! initialise 4-momentum
    ! set up 3 momentum in local frame
    local_3p = -[cos(sky_angles(2))*sin(sky_angles(1)), & ! x
    &           sin(sky_angles(2))*sin(sky_angles(1)), & ! y
    &           cos(sky_angles(1))]                      ! z
    ! transform 4-momentum in local fram [1, local_3p] to 4-momentum in proj coords
    ! basis forms chosen so that:
    !   - theta=0 is in the dr direction
    !   - theta=pi/2, phi=0 is in the dp direction
    !   - theta=pi/2, phi=pi/2 is in the dq direction
    four_momentum = [1., & ! t
    &               local_3p(3)/sqrt(sz_H(coords(1),coords(2),coords(3),coords(4))), &  ! r
    &               local_3p(1)/sqrt(sz_F(coords(1),coords(2),coords(3),coords(4))), &  ! p
    &               local_3p(2)/sqrt(sz_F(coords(1),coords(2),coords(3),coords(4)))]    ! q

    !! set state vector
    new_ray%state_vector = [coords, four_momentum, 0.,-1.]
  end function new_ray

  subroutine invert_ray(photon)
    ! perform inversion on a ray
    type(ray), intent(inout) :: photon
    real :: pq_norm
    real :: p, q

    ! print*,"INVERTING"
    ! transform coordinates
    p = photon%state_vector(3)
    q = photon%state_vector(4)
    pq_norm = photon%state_vector(3)**2 + photon%state_vector(4)**2
    photon%state_vector(3:4) = photon%state_vector(3:4)/pq_norm

    ! update inverted property
    photon%inverted = .not. photon%inverted
    photon%times_inverted = photon%times_inverted + 1

    ! transform four momentum
    photon%state_vector(7:8) = photon%state_vector(7:8)/pq_norm - &
    &        2*[p,q]*(p*photon%state_vector(7) + &
    &        q*photon%state_vector(8))/(pq_norm**2)

    ! dipole function inversion is handled by setting photon%inverted == true
  end subroutine invert_ray
end module ray_propagator
