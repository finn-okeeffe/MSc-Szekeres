subroutine initialise_fortran_subroutines()
  ! run all that initialisation stuff that is in main.f90 (copied to model_statistics)
  use model_statistics

  implicit none
  call pre_model_setup(.TRUE.)
end subroutine initialise_fortran_subroutines

subroutine solve_model_dynamics(params)
  use model_statistics
  implicit none
  real(8), intent(in) :: params(3)

  !! Set Szekeres model parameters
  call set_params(params)

  ! initialise sz_k
  call init_k_newt_scalar(r, K_ERROR_INTEGRATION_NUM, MAX_K_ITERS)
  ! solve for R(t,r)
  call find_arealR_array_omp(t_0, 0., num_steps, r)
end subroutine solve_model_dynamics

subroutine cmb_dip_and_quad(angular_coords, CMB_dipole, CMB_quadrupole)
  use model_statistics
  implicit none

  real(8), intent(in) :: angular_coords(3) ! (r, theta, phi)
  real(8), intent(out) ::  CMB_dipole, CMB_quadrupole ! in units K, and K^2 respectively

  real :: aniso_map(0:N_PIX-1), coords(3)
  complex :: CMB_a1m(-1:1)
  logical :: CMB_errored(0:N_PIX-1)
  integer :: i

  !! convert to projective coords
  coords = [angular_coords(1), & ! r
  &         sz_P(angular_coords(1)) + sz_S(angular_coords(1))*cos(angular_coords(3))*&
  &                  cos(angular_coords(2)/2)/sin(angular_coords(2)/2), & !p
  &         sz_Q(angular_coords(1)) + sz_S(angular_coords(1))*sin(angular_coords(3))*&
  &                  cos(angular_coords(2)/2)/sin(angular_coords(2)/2)] !q

  !! Raytracing
  call healpix_CMB_aniso([t_0, coords],lambda_step,aniso_map, CMB_errored)

  !! Compute Dipole amplitude and Quadrupole power
  if ( any(CMB_errored) ) then
    forall (i=0:N_PIX-1, CMB_errored(i)) aniso_map(i) = 0.
    call CMB_healpix_dipole_and_quadrupole(aniso_map, CMB_dipole, CMB_quadrupole)
  else
    call CMB_healpix_dipole_and_quadrupole(aniso_map, CMB_dipole, CMB_quadrupole)
  end if
end subroutine cmb_dip_and_quad

subroutine hubble_multipoles(angular_coords, dipoles, quadrupoles)
  use model_statistics
  implicit none
  real(8), intent(in), dimension(3) :: angular_coords ! (r, theta, phi)
  real(8), intent(out), dimension(12) :: dipoles, quadrupoles

  real :: coords(3), z_vals(NUM_DATA_POINTS), d_vals(NUM_DATA_POINTS), &
        sz_comp_C0(NUM_Z_POINTS), sz_comp_C1(NUM_Z_POINTS), sz_comp_C2(NUM_Z_POINTS)
  integer :: i

  ! convert angular coords to projective coords
  coords = [angular_coords(1), & ! r
  &         sz_P(angular_coords(1)) + sz_S(angular_coords(1))*cos(angular_coords(3))*&
  &                  cos(angular_coords(2)/2)/sin(angular_coords(2)/2), & !p
  &         sz_Q(angular_coords(1)) + sz_S(angular_coords(1))*sin(angular_coords(3))*&
  &                  cos(angular_coords(2)/2)/sin(angular_coords(2)/2)] !q

  ! perform ray tracing
  call raytrace_galaxies([t_0, coords], lambda_step, comp_rays)

  ! calculate hubble flow dipole
  z_vals = [(comp_rays(i)%state_vector(5)-1, i=1,NUM_DATA_POINTS)]
  d_vals = [(comp_rays(i)%d_L, i=1,NUM_DATA_POINTS)]
  call Hubble_Multipoles_true_d(z_vals, d_vals, sz_comp_C0, sz_comp_C1, sz_comp_C2)

  ! return dipole and quadrupole power at z values
  dipoles = sz_comp_C1
  quadrupoles = sz_comp_C2
end subroutine hubble_multipoles
