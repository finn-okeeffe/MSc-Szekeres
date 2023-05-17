! model_statistics is an outdated name. Most of the functions in this file got moved to more
! natural places, now it just holds the function to load data and initialise the code
! before we know the parameters, and holds some parameters for the grid, integration, and
! ray tracing.

module model_statistics

  use initialize_k
  use CMB_rays
  use omp_lib
  use trace_composite

  implicit none

  !! INTEGRATION AND RAYTRACING PARAMETERS
  integer, parameter :: N=5000 ! number of points in r grid (not including 0)
  ! ! 1000 is too low, significantly reduces accuracy
  integer, parameter :: num_steps = 10000 ! number of points to use in integration of R(t,r)
  ! doubling to 20000 doesn't increase accuracy much
  real, parameter :: lambda_step = -0.01
  integer, parameter :: K_ERROR_INTEGRATION_NUM = 5000 ! number of points to use when
  !           doing integration for finding error in k(r)
  integer, parameter :: MAX_K_ITERS = 100 ! number to max iterations when finding k(r)

  !! GRID AND COMPOSITE VARIABLES
  real, allocatable :: r(:) ! holds r grid values
  real, dimension(NUM_CATALOGUES, NUM_Z_POINTS) :: comp_C0, comp_C1, comp_C2
  ! logical :: errored(0:N_PIX-1)
  type(ray), allocatable, dimension(:) :: comp_rays
  real, dimension(NUM_Z_POINTS) :: mean_comp_c1, mean_comp_c2

  integer, private :: i

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Pre-Model Loading
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  subroutine pre_model_setup(use_saved_multipoles)
    ! setup grid, load COMPOSITE data, and calculate C_l for the observed COMPOSITE
    ! catalogue if the comp_Cl.dat file is missing (which stores these values)

    ! arguments
    logical, intent(in) :: use_saved_multipoles
    ! TRUE: load C_l values for observed COMPOSITE from comp_Cl.dat
    ! FALSE: recalculate C_l values for observed COMPOSITE, and save to comp_Cl.dat

    real :: z

    ! allocate arrays
    allocate(r(N))
    allocate(k_array(N+1))
    allocate(r_array(N+1))
    allocate(arealR_array(num_steps+1, N+1))
    allocate(t_array(num_steps+1))
    allocate(M_array(N))
    allocate(comp_rays(NUM_DATA_POINTS))

    ! initialise grid values
    r = [(500.*i/N, i=1,N)]
    r_array = [0., r]
    max_r = maxval(r)

    !! load COMPOSITE observational data (redshift, d_L, (l,b), uncertainties, etc)
    print*,"Loading Composite Data..."
    call load_composite_data()

    !! get C_l values for COMPOSITE
    print*,"Acquiring Dipole Analysis on Composite..."
    ! load from file
    if ( use_saved_multipoles ) then
      open(1, file="comp_Cl.dat", status="old")
      do i = 1, NUM_Z_POINTS
        read(1,*) z, comp_C0(:,i), comp_C1(:,i), comp_C2(:,i)
        print*,"z:",z
      end do
    
    ! calculate from observed quantities
    else
      print*,"Calculating Multipoles"
      call get_Hubble_Multipoles(cz_arr, comp_C0, comp_C1, comp_C2)

      print*, "Saving composite dipole"

      ! save to comp_Cl.dat
      open(1, file="comp_Cl.dat", status='replace')
      do i = 1, NUM_Z_POINTS
        write(1,*) z_points(i), comp_C0(:,i), comp_C1(:,i), comp_C2(:,i)
      end do
      close(1)
      print*,"finished"
    end if

    ! set average C_l variables in trace_composite.f90
    forall (i=1:NUM_Z_POINTS) mean_comp_c1(i) = sum(comp_C1(:,i))/NUM_CATALOGUES
    forall (i=1:NUM_Z_POINTS) mean_comp_c2(i) = sum(comp_C2(:,i))/NUM_CATALOGUES
  end subroutine pre_model_setup
end module model_statistics
