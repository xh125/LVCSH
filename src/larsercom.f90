module lasercom
  use kinds,only : dp
  implicit none
  
  real(kind=dp)   :: w_laser     ! the laser energy in eV.
  real(kind=dp)   :: fwhm, fwhm_T            
  ! The laser shape is assumed to be a Gaussian f(t)=exp(-t^2/2T^2)
  !T is related to the full width at half-maximum as fwhm = 2sqrt(2ln2)T
  real(kind=dp)   :: Efield_x,Efield_y,Efield_z
  ! The laser Electric field strength along the xyz direction
  contains
  
end module lasercom