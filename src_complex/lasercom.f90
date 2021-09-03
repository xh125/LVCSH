module lasercom
  use kinds,only : dp
  implicit none
  
  logical         :: llaser         ! If .TRUE. a laser electric field  is applied.
  real(kind=dp)   :: efield         ! Amplitude of the laser electric field (in Ry a.u.;1 a.u. = 36.3609*10^10 V/m)
  real(kind=dp)   :: efield_cart(3) ! laser electric field (in Ry a.u.=36.3609*10^10 V/m) in cartesian axis.
  real(kind=dp)   :: w_laser        ! the center energy of laser in eV.
  real(kind=dp)   :: fwhm, fwhm_2T2   !         
  ! The laser shape is assumed to be a Gaussian f(t)=exp(-t^2/2T^2)
  !T is related to the full width at half-maximum as fwhm = 2sqrt(2ln2)T
  ! F(w)=exp(-(w-w_laser)^2*T^2/2)

  real(kind=dp),allocatable :: W_cvk(:,:,:) !光激发下的跃迁几率大小
  contains
  
  real function f_t(t,fwhm_2T2)
    implicit none
    real(kind=dp), intent(in) :: t
    real(kind=dp), intent(in) :: fwhm_2T2
    ! fwhm_2T2 = fwhm**2.0/4.0*log(2.0)
    f_t = exp(-1.0*t**2/fwhm_2T2)
    return
  end function f_t
  
  
  !ref : 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
  real function f_w(w,fwhm_2T2)
    implicit none
    real(kind=dp),intent(in) :: w
    real(kind=dp),intent(in) :: fwhm_2T2
    real(kind=dp) :: wfwhm
    wfwhm = (w-w_laser)**2 *fwhm_2T2
    wfwhm = -0.25*wfwhm
    f_w = EXP(wfwhm)
    !f_w = exp(-1.0*(w-w_laser)**2 * fwhm_2T2/4.0)
    return
  end function
  
  
end module lasercom