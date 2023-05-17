module spherical_harmonics
  ! contains methods to evaluate needed spherical harmonic functions
  use szmodel
  implicit none

  complex, private, parameter :: i = (0, 1) ! imaginary number i
contains

  !! DIPOLE FUNCTIONS
  ! using same normalisation as scipy
  function Y_1m1(theta,phi)
    ! m=-1, l=1
    real :: theta(:), phi(:)
    complex, dimension(size(phi)) :: Y_1m1

    Y_1m1 = 0.5*sqrt(3/(2*pi))*sin(theta)*exp(-i*phi)
  end function Y_1m1

  function Y_10(theta,phi)
    ! m=0, l=1
    real :: theta(:), phi(:)
    complex, dimension(size(phi)) :: Y_10

    Y_10 = 0.5*sqrt(3/pi)*cos(theta)
  end function Y_10

  function Y_11(theta,phi)
    ! m=1, l=1
    real :: theta(:), phi(:)
    complex, dimension(size(phi)) :: Y_11

    Y_11 = -0.5*sqrt(3/(2*pi))*sin(theta)*exp(i*phi)
  end function Y_11
end module spherical_harmonics
