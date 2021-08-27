module getwcvk
  use kinds,only    : dp,dpc
  use lasercom,only : W_cvk,efield_cart
  implicit none
  
  contains
  !ref : 1 S. Butscher et al., Physical Review B 72 (2005) 
  !ref : 2 <固体物理> (9-29)(9-31)
  subroutine get_Wcvk(ihband_min,ieband_max,fwhm,w_center)
    !得到光激发下垂直跃迁的跃迁几率
    use elph2,only  : vmef,nkf  !vmef(3,nbndsub,nbndsub,nkf)
    use readepw,only : etf,icbm
    implicit none
    integer , intent(in) :: ihband_min,ieband_max
    real(kind=dp),intent(in) :: fwhm
    real(kind=dp),intent(in) :: w_center
    
    real(kind=dp) :: fwhm_2T2
    
    ! W_cvk(cband,vband,ik)=|<E dot vmef(3,cband,vband,ik)>|^2 *f_w(w)
    !!光激发下的跃迁几率大小
    real(kind=dp) :: fcw
    complex(kind=dpc) :: Evmef
    real(kind=dp) :: E_mnk   !垂直激发能量
    integer :: ik,ikk,ibnd,jbnd,ipol
    integer :: ierr
    integer :: ivbm
    
    ivbm = icbm-1
    
    allocate(W_cvk(icbm:ieband_max,ihband_min:ivbm,nkf),stat=ierr)
    if(ierr /=0) call errore('getmcvk','Error allocating W_cvk',1)
    
    !ref : 1 S. Butscher et al., Physical Review B 72 (2005) 
    !ref : 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
    fwhm_2T2 = fwhm**2.0/4.0*log(2.0)
    W_cvk = 0.0
    do ik=1,nkf
      do ibnd=ihband_min,ivbm
        do jbnd=icbm,ieband_max
          Evmef = 0.0
          E_mnk = etf(jbnd,ik)-etf(ibnd,ik)
          do ipol=1,3
            Evmef =Evmef+ (efield_cart(ipol)*(vmef(ipol,jbnd,ibnd,ik)))
          enddo
          fcw = f_w(E_mnk,w_center,fwhm_2T2)
          W_cvk(jbnd,ibnd,ik) = REAL(Evmef*CONJG(Evmef))*fcw
        enddo
      enddo
    enddo  
    
  end subroutine get_Wcvk
  
  !ref : 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
  real function f_w(w,w_laser,fwhm_2T2)
    implicit none
    real(kind=dp),intent(in) :: w,w_laser
    real(kind=dp),intent(in) :: fwhm_2T2
    real(kind=dp) :: wfwhm
    wfwhm = (w-w_laser)**2 *fwhm_2T2
    wfwhm = -0.25*wfwhm
    f_w = EXP(wfwhm)
    !f_w = exp(-1.0*(w-w_laser)**2 * fwhm_2T2/4.0)
    return
  end function

  
end module getwcvk