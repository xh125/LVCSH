module getmcvk
  use lasercom
  implicit none
  
  contains
  !ref : 1 S. Butscher et al., Physical Review B 72 (2005) 
  !ref : 2 <固体物理> (9-31)
  subroutine get_Mcvk()
    !得到光激发下垂直跃迁的跃迁几率
    use elph2,only  : vmef,ibndmin,ibndmax,nbndfst,nkf  !vmef(3,nbndsub,nbndsub,nkf)
    use readepw,only : E_nk
    implicit none
    
    ! W_cvk(cband,vband,ik)=|<E dot vmef(3,cband,vband,ik)>|^2 *f_w(w)
    real(kind=dp) :: fcw
    complex :: Evmef
    real(kind=dp) :: W_mnk   !垂直激发能量
    integer :: ik,ivband,icband,ipol
    integer :: ierr
    if(.not. allocated(W_cvk)) then
      allocate(W_cvk(nbndfst,nbndfst,nkf),stat=ierr)
      if(ierr /=0) call errore('getmcvk','Error allocating W_cvk',1)
    endif
    
    !ref : 1 S. Butscher et al., Physical Review B 72 (2005) 
    !ref : 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
    fwhm_2T2 = fwhm**2.0/4.0*log(2.0)
    W_cvk = 0.0
    do ik=1,nkf
      do ivband=1,nbndfst-1
        do icband=ivband+1,nbndfst
          Evmef = 0.0
          W_mnk = E_nk(icband,ik)-E_nk(ivband,ik)
          do ipol=1,3
            Evmef =Evmef+ (efield_cart(ipol)*(vmef(ipol,icband+ibndmin-1,ivband+ibndmin-1,ik)))
          enddo
          fcw = f_w(W_mnk)
          W_cvk(icband,ivband,ik) = Evmef*CONJG(Evmef)*fcw
        enddo
      enddo
    enddo  
    
  end subroutine get_Mcvk  
end module getmcvk