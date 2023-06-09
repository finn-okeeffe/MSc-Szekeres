!! Contains methods and variables to simulate a COMPOSITE catalogue (utilising ray_propagator.f90) and analyse COMPOSITE catalogue
module trace_composite
  use ray_propagator
  use healpix_modules
  use omp_lib
  implicit none

  ! location of COMPOSITE data
  character(len=25), parameter :: COMPOSITE_FILENAME = "data/Composite_LG.txt"
  ! number of data points
  integer, parameter :: NUM_DATA_POINTS = 4534
  ! location of the CMB in galactic latitude and longitude
  real, parameter, dimension(2) :: COMPOSITE_CMB_DIPOLE = [-48.26+90, 236.99]*pi/180
  ! variables to hold COMPOSITE data
  real, dimension(NUM_DATA_POINTS) :: cz_arr, d_arr, sigma_v_arr, sigma_d_arr
  real(8), dimension(NUM_DATA_POINTS) :: theta_arr, phi_arr

  !! Composite_LG.txt layout:
  ! cz (km/sec)     d (Mpc/h)     v_pec (km/s)    sigma_v     l (degrees)     b(degrees)

  !! notes on data:
  ! v_pec = v - 100d

  !! Smoothing Parameters for calculating angular power using BNW method
  real, parameter :: SIGMA_Z = 0.01
  real, parameter :: SIGMA_THETA = 25*pi/180

  !! Dipole Variation Parameters and Variables
  integer :: m
  integer, parameter :: NUM_Z_POINTS = 12 ! number of z points to sample at when calculating the angular power
  real, parameter, dimension(NUM_Z_POINTS) :: z_points = [(m*0.045/(NUM_Z_POINTS-1), m=0,(NUM_Z_POINTS-1))] ! the specific z points to sample, ranging from z=0 to z=0.045
  integer, parameter :: COMP_N_SIDE = 2**3 ! HEALPix grid parameter for calculating COMPOSITE angular power, only use 2**0 for fast prototyping
  integer, parameter :: COMP_N_PIX = 12*COMP_N_SIDE**2 ! number of pixels in HEALPix map
  integer, parameter :: NUM_CATALOGUES = 10000 ! number of random mock catalogues to use when evaluating uncertainties a la BNW
  real, allocatable :: dL_cats(:,:) ! holds random catalogues of d_L corresponding to the COMPOSITE catalogue and its uncertainties

  !! Holds weights for fast calculation of H dipole
  real, dimension(0:COMP_N_PIX-1,NUM_DATA_POINTS) :: wangs

contains

  subroutine load_composite_data()
    ! reads composite data from COMPOSITE_FILENAME, with NUM_DATA_POINTS data points.

    real, dimension(NUM_DATA_POINTS) :: vpec_arr, l_arr, b_arr ! COMPOSITE data that needs to be transformed

    !! other local variables
    integer :: i, pix_num ! looping integer
    real(8) :: sky_pix_vec(3), sky_comp_vec(3), delta_ang ! angle and vector for computing w_ang

    !! read observed COMPOSITE data into arrays
    open(1, file=COMPOSITE_FILENAME, status='old')
    do i = 1, NUM_DATA_POINTS
      read(1,*) cz_arr(i), d_arr(i), vpec_arr(i), sigma_v_arr(i), l_arr(i), b_arr(i)
    end do
    close(1)

    ! convert sky angles from degrees to radians
    phi_arr = l_arr * pi/180
    theta_arr = (-b_arr + 90) * pi/180

    ! convert cz from km/s to unitless (c=1)
    cz_arr = cz_arr / (3e5)

    ! get uncertainty in distance measure from uncertainty in peculiar velocity
    sigma_d_arr = sigma_v_arr/100

    !! Setup random dL catalogues, for faster calculation of H_0
    allocate(dL_cats(NUM_CATALOGUES, NUM_DATA_POINTS))
    ! get random d_L arrays
    dL_cats(1,:) = d_arr
    do i = 2, NUM_CATALOGUES
      dL_cats(i,:) = rand_dL_array()
    end do

    !! precalculate angle weightings for faster calculation of H_0
    do pix_num = 0, COMP_N_PIX-1
      call pix2vec_ring(COMP_N_SIDE, pix_num, sky_pix_vec)
      do i = 1, NUM_DATA_POINTS
        call ang2vec(theta_arr(i),phi_arr(i), sky_comp_vec)
        call angdist(sky_pix_vec, sky_comp_vec, delta_ang)
        wangs(pix_num,i) = exp(-0.5*(delta_ang/SIGMA_THETA)**2)
      end do
    end do
  end subroutine load_composite_data

  subroutine raytrace_galaxies(obs_coords, lambda_step, comp_rays)
    !! raytrace composite data using ray_propagator.f90

    !! arguments
    real, dimension(4), intent(in) :: obs_coords ! location of observer, initial position of rays
    real, intent(in) :: lambda_step ! step size to use when raytracing
    type(ray), intent(out), dimension(NUM_DATA_POINTS) :: comp_rays ! holds fully evolved ray types to return

    ! local variables
    type(ray) :: private_ray
    integer :: i

    !! Create and propagate rays, using OpenMP since the result of one ray does not depend on another ray
    !$OMP PARALLEL PRIVATE(private_ray, i) SHARED(comp_rays)
      !$OMP DO
      do i = 1, NUM_DATA_POINTS
        !! initialise_ray
        private_ray = new_ray([theta_arr(i), phi_arr(i)], obs_coords, .FALSE.)

        !! PROPAGATE RAY
        call propagate_ray(private_ray, d_arr(i), lambda_step, .FALSE.)

        !! save to ray array
        comp_rays(i) = private_ray
      end do
      !$OMP END DO
    !$OMP END PARALLEL
  end subroutine raytrace_galaxies


  function rand_dL_array()
    ! gets a sample of d_L from the COMPOSITE d_L normal distributions
    real, dimension(NUM_DATA_POINTS) :: rand_dL_array ! return value: holds a single random realisation of each d_L in COMPOSITE
    real, dimension(2) :: r ! random number
    integer :: i ! loop index

    ! loop through each galaxy in COMPOSITE
    do i=1,NUM_DATA_POINTS
      call random_number(r) ! get a random number between 0 and 1, uniform distribution
      rand_dl_array(i) = max(sigma_d_arr(i)*sqrt(-2*log(r(1)))*cos(2*pi*r(2))+d_arr(i), 0.1) ! convert to normal distribution for a galaxy in COMPOSITE
      !     realisation of d_L is forced to be atleast 0.1 by `max' to avoid a nonpositive distance.
    end do
  end function rand_dL_array

  function fast_smoothed_hubble(zindex,hpindex,z_vals,d_vals)
    !! calculates the smoothed Hubble function as per BNW(2016) page 9

    !! arguments
    integer, intent(in) :: zindex, hpindex ! indices of z_points and HEALPix grid respectively, corresponding to arguments of H_0(l,b,z)
    real, intent(in), dimension(NUM_DATA_POINTS) :: z_vals, d_vals ! data from catalogue, be it observed or simulated

    !! return value
    real :: fast_smoothed_hubble

    !! local variables
    real, dimension(NUM_DATA_POINTS) :: wz, wd, H, zeta ! factors of H_0
    real :: z ! redshift

    !! get z value, and calculate redshift weight
    z = z_points(zindex)
    wz = exp(-0.5*((z-z_vals)/SIGMA_Z)**2)

    !! calculate remaining factors
    zeta = z_vals + 0.5*(1-q0)*z_vals**2 - (1-q0-3*q0**2+j0)*z_vals**3
    wd = zeta*d_vals/(sigma_d_arr**2)
    H = zeta/d_vals

    !! Evaluate H(theta,phi,z) using wangs precalculated in load_composite_data
    fast_smoothed_hubble = sum(H*wd*wz*wangs(hpindex,:))/sum(wd*wz*wangs(hpindex,:))
  end function fast_smoothed_hubble


  subroutine get_Hubble_Multipoles(z_vals, comp_C0, comp_C1, comp_C2)
    !! Calculate the mean and variance of C1 of smoothed Hubble over z for a given catalogue
    !! arguments
    real, intent(in), dimension(NUM_DATA_POINTS) :: z_vals ! redshift of catalogue
    real, intent(out), dimension(NUM_CATALOGUES, NUM_Z_POINTS) :: comp_C0, comp_C1, comp_C2 ! angular power to return

    ! HEALPix parameters
    integer, parameter :: L_MAX = 2
    integer, parameter :: M_MAX = L_MAX

    ! local variables
    real :: current_dL(NUM_DATA_POINTS) ! a chosen d_L catalogue
    complex(8) :: alm(1:1,0:L_MAX,0:M_MAX) ! holding returned values of a^l_m
    real(8) :: cl(0:L_MAX,1) ! holding returned values of C_l
    real :: H_avg
    real, dimension(0:COMP_N_PIX-1) :: curr_H
    real(8), dimension(0:COMP_N_PIX-1) :: aniso_map

    integer :: i,j,k

    !! calculate smooth H at representative points for original COMPOSITE
    ! loop through catalogues using OpenMP
    print*,"Starting loop..."
    ! calculate angular power for each catalogue
    !$OMP PARALLEL DO PRIVATE(i,j,k,current_dL, alm, cl, curr_H) SHARED(comp_C1, comp_C0, comp_C2) COLLAPSE(2)
    do i = 1, NUM_CATALOGUES
      current_dL = dL_cats(i,:) ! select a catalogue
      ! calculate angular power at each z for the current catalogue
      do j = 1, NUM_Z_POINTS
        curr_H = 0.
        ! loop through points on shell
        ! calculate H(theta,phi,z) for each point on shell
        do k = 1, 12*COMP_N_SIDE**2
          curr_H(k-1) = fast_smoothed_hubble(j,k-1,z_vals,current_dL)
        end do
        ! get average
        H_avg = sum(curr_H)/COMP_N_PIX
        ! get anisotropy Delta H/H_0
        aniso_map = (curr_H - H_avg)/H_avg
        ! get monopole, dipole, and quadrupole power
        call map2alm(COMP_N_SIDE, L_MAX, M_MAX, aniso_map, alm)
        call alm2cl(L_MAX, M_MAX, alm, cl)
        comp_C0(i,j) = cl(0,1) ! should be 0 down to float precision, since aniso_map has zero average by definition
        comp_C1(i,j) = cl(1,1)
        comp_C2(i,j) = cl(2,1)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine get_Hubble_Multipoles


  subroutine Hubble_Multipoles_true_d(z_vals, d_vals, comp_C0, comp_C1, comp_C2)
    ! Calculate the mean and variance of C1 of smoothed Hubble over z for an exactly known d_L (e.g., simulated COMPOSITE catalogue)
    !! Arguments
    real, intent(in), dimension(NUM_DATA_POINTS) :: z_vals, d_vals ! z and d_L of catalogue
    real, intent(out), dimension(NUM_Z_POINTS) :: comp_C0, comp_C1, comp_C2 ! angular power to return

    ! HEALPix parameters
    integer, parameter :: L_MAX = 2
    integer, parameter :: M_MAX = L_MAX

    ! local variables
    complex(8) :: alm(1:1,0:L_MAX,0:M_MAX) ! holding returned values of a^l_m
    real(8) :: cl(0:L_MAX,1) ! holding returned values of C_l
    real :: H_avg
    real, dimension(0:COMP_N_PIX-1) :: curr_H
    real(8), dimension(0:COMP_N_PIX-1) :: aniso_map

    integer :: j,k

    !! calculate smooth H at representative points for original COMPOSITE
    ! loop through catelogs
    ! loop through a number of sample z points
    !$OMP PARALLEL DO PRIVATE(j,k,alm,cl,curr_H,H_avg,aniso_map) SHARED(comp_C0,comp_C1,comp_C2)
    do j = 1, NUM_Z_POINTS
      ! loop through points on shell
      curr_H = 0.
      ! calculate H(theta,phi,z) for each point on shell
      do k = 1, 12*COMP_N_SIDE**2
        curr_H(k-1) = fast_smoothed_hubble(j,k-1,z_vals,d_vals)
      end do
      ! get average
      H_avg = sum(curr_H)/COMP_N_PIX
      ! get anisotropy Delta H/H_0
      aniso_map = (curr_H - H_avg)/H_avg
      ! get angular power
      call map2alm(COMP_N_SIDE, L_MAX, M_MAX, aniso_map, alm)
      call alm2cl(L_MAX, M_MAX, alm, cl)
      comp_C0(j) = cl(0,1)
      comp_C1(j) = cl(1,1)
      comp_C2(j) = cl(2,1)
    end do
    !$OMP END PARALLEL DO
  end subroutine Hubble_Multipoles_true_d
end module trace_composite
