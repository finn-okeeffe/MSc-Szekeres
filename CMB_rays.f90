module CMB_rays
  use ray_propagator
  use omp_lib
  use healpix_modules
  use healpix_types
  implicit none

  ! Parameters for HEALpix discretisation
  integer, parameter :: N_SIDE = 2**3 ! number of halvings
  integer, parameter :: N_PIX = 12*N_SIDE**2 ! number of pixels

  ! Measured CMB
  real, parameter :: MEAN_CMB_TEMP = 2.725 !K
  real, parameter :: TRUE_CMB_D1 = (5.77e-3)/MEAN_CMB_TEMP ! CMB DeltaT/T dipole amplitude
  real, parameter :: SIGMA_CMB_D1 = (0.36e-3)/MEAN_CMB_TEMP ! Uncertainty in above

contains

  subroutine healpix_CMB_aniso(obs_coords, lambda_step, aniso_map, errored)
    ! calculates CMB anisotropy healpix map in ring order

    !! arguments
    real, intent(in) :: lambda_step, obs_coords(4)
    real, dimension(0:N_PIX-1), intent(out) :: aniso_map
    logical, dimension(0:N_PIX-1) :: errored

    !! local variables:
    real(8) :: theta, phi
    real :: avg
    integer :: pix_num
    type(ray) :: private_ray

    real :: coords(4), fourmom(4)

    errored = .FALSE.
    !! Create and propagate rays
    open(2, file="CMB_state_vectors.dat", status="replace")
    !$OMP PARALLEL PRIVATE(private_ray, phi, theta, pix_num) SHARED(aniso_map, errored)
      !$OMP DO
      do pix_num = 0, N_PIX-1
        !! initialise_ray
        ! get angle
        call pix2ang_ring(N_SIDE, pix_num, theta, phi)
        ! initialise ray using new_ray in ray_propagator
        private_ray = new_ray([theta, phi], obs_coords, .FALSE.)

        !! PROPAGATE RAY
        call propagate_ray(private_ray, -1., lambda_step, .FALSE.)

        !! calculate factor accounting for remaining redshift
        if ( .not. private_ray%errored ) then
          aniso_map(pix_num) = max_r / (private_ray%state_vector(5) * &
          & sz_arealR(private_ray%state_vector(1),max_r)) ! 1/((1+z)*a(t))

          if ( private_ray%inverted ) then
            call invert_ray(private_ray)
          end if

          coords = private_ray%state_vector(1:4)
          fourmom = private_ray%state_vector(5:8)
          
          write(2,*) pix_num, theta, phi, inner_product(fourmom,fourmom,coords), private_ray%state_vector
        else
          print*,"Errored, pix num:", pix_num,"radius:",private_ray%state_vector(2)
          print*,"error message:", private_ray%err_message
        end if
        errored(pix_num) = private_ray%errored
      end do
      !$OMP END DO
    !$OMP END PARALLEL
    close(2)


    !! calculate anisotropy
    avg = sum(aniso_map, mask=.not.errored)/count(.not.errored)
    ! avg = sum(aniso_map)/N_PIX
    aniso_map = (aniso_map - avg)/avg ! = Delta T / T
  end subroutine healpix_CMB_aniso


  subroutine CMB_healpix_dipole_and_quadrupole(aniso_map, dipole, quadrupole)
    real(dp), intent(in) :: aniso_map(0:N_PIX-1) ! Delta T / T
    real, intent(out) ::  dipole, quadrupole
    integer(i4b), parameter :: lmax = 2
    complex(dpc) :: alm(1, 0:lmax, 0:lmax)
    real(dp) :: cl(0:lmax,1)

    call map2alm(N_SIDE, lmax, lmax, aniso_map, alm)
    call alm2cl(lmax, lmax, alm, cl)

    quadrupole = cl(2,1)*6*(MEAN_CMB_TEMP**2)/(2*pi) !K^2
    dipole = 1.5*sqrt(cl(1,1)/pi)*MEAN_CMB_TEMP ! K
    
  end subroutine CMB_healpix_dipole_and_quadrupole
end module CMB_rays
