program main
  use model_statistics

  implicit none
  real :: param_error, angular_coords(3), coords(3), params(3), axis_theta
  integer :: i


    !! TESTING VARS
    real :: ray_theta, ray_phi
    type(ray) :: test_ray

  call pre_model_setup(.TRUE.)


  ! parameters: delta_0, alpha, r_0
  ! [-0.97174837  0.94948474 34.97985712 34.96396796  1.99449534] nlp: ~ 29.90 for COMP dip and CMB only
  ! [-0.58544466  0.85544953 57.57741181 55.77230181  2.63266338] nlp: ~ 28.50 for COMP dip and CMB only
  ! [-0.61739258  0.95       62.37448863 53.64650494  1.67415252] modified nlp: ~ 8.527 (COMP dip nlp/num_z + CMB dip nlp + CMB quad nlp)
 ! ([-0.80, 0.9, 55.0, 54.0, 0.7*np.pi]
  ! [-0.53404086  0.95       64.88764856 48.27292444  1.87559392] modified nlp: ~7.766
  params = [-0.53404086, 0.95, 64.88764856]
  call set_params(params)

  ! angular_coords: r, theta, phi
  angular_coords =  [48.27292444, 1.87559392, 0.0]


  ! convert from (r,theta,phi) to (r,p,q)
  coords = [angular_coords(1), & ! r
  &         sz_P(angular_coords(1)) + sz_S(angular_coords(1))*cos(angular_coords(3))*&
  &                  cos(angular_coords(2)/2)/sin(angular_coords(2)/2), & !p
  &         sz_Q(angular_coords(1)) + sz_S(angular_coords(1))*sin(angular_coords(3))*&
  &                  cos(angular_coords(2)/2)/sin(angular_coords(2)/2)] !q


  print*,"Projective Coords:",coords
  print*,"Parameters:", params

  print*,"executing saving iteration:"
  param_error = parameter_error(params, coords, .TRUE., .TRUE.)
  print*,"error:",param_error

  print*,"Doing Mocks"
  call save_composite_mocks_multipoles(params, coords, .FALSE.)
  print*,"Done!"

  ! test ray
  ! print*,"Solving Dynamics"
  ! call init_k_newt_scalar(r, K_ERROR_INTEGRATION_NUM, MAX_K_ITERS)
  ! call find_arealR_array_omp(t_0, 0., num_steps, r)
  ! print*,"Creating Ray"
  ! call pix2ang_ring(N_SIDE, 674, ray_theta, ray_phi)
  ! test_ray = new_ray([ray_theta, ray_phi], [t_0, coords], .FALSE.)
  ! allocate(test_ray%stv_history(test_steps, 10))
  ! allocate(test_ray%aff_param_hist(test_steps))
  ! allocate(test_ray%inv_hist(test_steps))
  ! print*,"stv_history shape:",shape(test_ray%stv_history)
  ! print*,"aff_param_hist shape:",shape(test_ray%aff_param_hist)
  ! print*,"inv_hist shape:",shape(test_ray%inv_hist)
  ! print*,"inv_hist allocated:",allocated(test_ray%inv_hist)
  ! print*,"initial state vector:",test_ray%state_vector
  ! print*,"Propagating Ray"
  ! call propagate_ray(test_ray, -1., lambda_step, .TRUE.)
  ! print*,"Saving Ray"
  ! open(1, file="ray_history.dat", status='replace')
  ! do i = 1, test_steps
  !   write(1,*) test_ray%aff_param_hist(i), test_ray%stv_history(i,:), test_ray%inv_hist(i)
  ! end do
  ! close(1)


contains
  function parameter_error(parameters, coordinates, initialise, save)
    ! evaluates error function for a chosen model
    real, intent(in) :: parameters(3), coordinates(3) ! parameters and coords
    real :: parameter_error
    ! params: (delta_0, alpha, r_0, D)
    ! coords: (r, p, q)
    logical, intent(in) :: save ! whether to save results
    logical, intent(in) :: initialise ! whether to recalculate k(r) and R(t,r)

    !! ERROR WEIGHTS
    real, parameter :: CMB_D1_WEIGHT = 1
    real, parameter :: HUB_C1_WEIGHT = 1
    real, parameter :: HUB_C2_WEIGHT = 1


    !! local variables
    real :: obs_coords(4), aniso_map(0:N_PIX-1), z_vals(NUM_DATA_POINTS), &
    &     d_vals(NUM_DATA_POINTS), sz_comp_C0(NUM_Z_POINTS), &
    &     sz_comp_C1(NUM_Z_POINTS), sz_comp_C2(NUM_Z_POINTS), &
    &     CMB_dipole, CMB_quad
    real(8) :: theta, phi
    complex :: CMB_a1m(-1:1)
    integer :: j
    integer :: istart, ifinish ! for timing
    logical :: CMB_errored(0:N_PIX-1)

    !! Set Szekeres model parameters
    call set_params(parameters)

    !! Observer Coordinates
    obs_coords = [t_0, coordinates(1), coordinates(2), coordinates(3)]

    !! solve dynamics if needed
    if ( initialise ) then
      ! initialise sz_k
      call system_clock(istart)
      call init_k_newt_scalar(r, K_ERROR_INTEGRATION_NUM, MAX_K_ITERS)
      call system_clock(ifinish)
      print*, "k(r) Initialisation time:",(ifinish-istart)/1000
      ! evolve R
      ! call find_arealR_array(t_0, t_0/2, num_steps, [0., r])
      call system_clock(istart)
      call find_arealR_array_omp(t_0, 0., num_steps, r)
      call system_clock(ifinish)
      print*, "R(t,r) Initialisation time:",(ifinish-istart)/1000

      ! ! save k(r) and R(t,r) if needed
      ! if ( save ) then
      !   print*,"Saving k(r_j) and R(t_i,r_j)"
      !   open(1, file="k_arr.dat", status="replace")
      !   open(2, file="R_arr.dat", status="replace")
      !   write(2,*) t_array
      !   do j = 1, N+1
      !     write(1,*) r_array(j), k_array(j)
      !     write(2,*) arealR_array(:,j)
      !   end do
      !   close(1)
      !   close(2)
      !   print*,"done"
      ! end if

    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        CMB
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! trace CMB
    call system_clock(istart)
    call healpix_CMB_aniso(obs_coords,lambda_step,aniso_map, CMB_errored)
    call system_clock(ifinish)
    print*, "CMB Raytracing time:",(ifinish-istart)/1000

    !! Compute dipole
    call system_clock(istart)
    if ( any(CMB_errored) ) then
      forall (i=0:N_PIX-1, CMB_errored(i)) aniso_map(i) = 0.
      call CMB_healpix_dipole_and_quadrupole(aniso_map, CMB_dipole, CMB_quad)
    else
      call CMB_healpix_dipole_and_quadrupole(aniso_map, CMB_dipole, CMB_quad)
    end if
    print*,"CMB DIPOLE:",CMB_dipole,"K"
    print*,"CMB QUADRUPOLE:",CMB_quad,"K^2"
    call system_clock(ifinish)
    print*, "CMB dipole calculation time:",(ifinish-istart)/1000

    ! save aniso array
    if ( save ) then
      print*, "Saving aniso_map"
      open(1, file="CMB_aniso_map.dat", status='replace')
      do j = 0, N_PIX-1
        call pix2ang_ring(N_SIDE, j, theta, phi)
        write(1,*) theta, phi, aniso_map(j)
      end do
      close(1)
      print*,"finished"
    end if



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        COMPOSITE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Ray trace composite
    print*,"tracing galaxies"
    call system_clock(istart)
    call raytrace_galaxies(obs_coords, lambda_step, comp_rays)
    call system_clock(ifinish)
    print*, "COMPOSITE raytracing time:",(ifinish-istart)/1000


    call system_clock(istart)
    z_vals = [(comp_rays(i)%state_vector(5)-1, i=1,NUM_DATA_POINTS)]
    d_vals = [(comp_rays(i)%d_L, i=1,NUM_DATA_POINTS)]
    call Hubble_Multipoles_true_d(z_vals, d_vals, sz_comp_C0, sz_comp_C1, sz_comp_C2)
    call system_clock(ifinish)
    print*, "Hubble power spectrum time:",(ifinish-istart)/1000

    if ( save ) then
      print*,"Saving mock composite z and d's"
      open(1, file="sz_z_d.dat", status="replace")
      do j = 1, NUM_DATA_POINTS
        write(1,*) z_vals(j), d_vals(j)
      end do
      close(1)
    end if

    if ( save ) then
      print*, "Saving mock composite dipole"
      open(1, file="sz_comp_C1_C2.dat", status='replace')
      do j = 1, NUM_Z_POINTS
        write(1,*) z_points(j), sz_comp_C0(j), sz_comp_C1(j), sz_comp_C2(j)
      end do
      close(1)
      print*,"finished"
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        EVALUATE ERROR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parameter_error = HUB_C1_WEIGHT*sum(abs((sz_comp_C1 - mean_comp_c1)/mean_comp_c1))/NUM_Z_POINTS+ &
    &     HUB_C2_WEIGHT*sum(abs((sz_comp_c2 - mean_comp_c2)/mean_comp_c2))/NUM_Z_POINTS + &
    &     CMB_D1_WEIGHT*abs((CMB_dipole/MEAN_CMB_TEMP - TRUE_CMB_D1)/TRUE_CMB_D1)
    print*,"--RELATIVE ERROR, NOT RELATED TO DISTRIBUTION--"
    print*,"HUBBLE DIPOLE ERROR:", (sz_comp_C1 - mean_comp_c1)/mean_comp_c1
    print*,"HUBBLE QUADRUPOLE MAX ERROR:", (sz_comp_c2 - mean_comp_c2)/mean_comp_c2
  end function parameter_error



  subroutine save_composite_mocks_multipoles(parameters, coordinates, initialise)
    real, intent(in) :: parameters(3), coordinates(3)
    logical, intent(in) :: initialise

    !! local variables
    real :: obs_coords(4), z_vals(NUM_DATA_POINTS), &
    &     d_vals(NUM_DATA_POINTS), sz_comp_C0(NUM_CATALOGUES, NUM_Z_POINTS), &
    &     sz_comp_C1(NUM_CATALOGUES, NUM_Z_POINTS), sz_comp_C2(NUM_CATALOGUES, NUM_Z_POINTS)
    integer :: istart, ifinish ! for timing
    integer :: j


    !! Observer Coordinates
    obs_coords = [t_0, coordinates(1), coordinates(2), coordinates(3)]

    !! solve dynamics if needed
    if ( initialise ) then
      call set_params(parameters)

      ! initialise sz_k
      call system_clock(istart)
      call init_k_newt_scalar(r, K_ERROR_INTEGRATION_NUM, MAX_K_ITERS)
      call system_clock(ifinish)
      print*, "k(r) Initialisation time:",(ifinish-istart)/1000

      ! evolve R
      call system_clock(istart)
      call find_arealR_array_omp(t_0, 0., num_steps, r)
      call system_clock(ifinish)
      print*, "R(t,r) Initialisation time:",(ifinish-istart)/1000
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        COMPOSITE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Ray trace composite
    print*,"tracing galaxies"
    call system_clock(istart)
    call raytrace_galaxies(obs_coords, lambda_step, comp_rays)
    call system_clock(ifinish)
    print*, "COMPOSITE raytracing time:",(ifinish-istart)/1000


    call system_clock(istart)
    z_vals = [(comp_rays(i)%state_vector(5)-1, i=1,NUM_DATA_POINTS)]
    d_vals = [(comp_rays(i)%d_L, i=1,NUM_DATA_POINTS)]
    call get_Hubble_Multipoles(z_vals, sz_comp_C0, sz_comp_C1, sz_comp_C2)
    call system_clock(ifinish)
    print*, "Hubble power spectrum time:",(ifinish-istart)/1000

    print*, "Saving COMPOSITE mocks multiples"
    open(1, file="sz_comp_Cl_mocks.dat", status='replace')
    do j = 1, NUM_Z_POINTS
      write(1,*) z_points(j), sz_comp_C0(:,j), sz_comp_C1(:,j), sz_comp_C2(:,j)
    end do
    close(1)
    print*,"finished"
    
  end subroutine save_composite_mocks_multipoles
end program main
