module dynamics
  use kinds,only : dp,dpc
  use io,only :
  use parameters,only : temp
  !use elph2,only          : nktotf,nqtotf
  !use modes,only          : nmodes
  use epwcom,only         : kqmap
  
  implicit none
  
  contains

  subroutine set_gamma(nmodes,nq,gamma,ld_fric,wf,ld_gamma)
    use constants,only : eps10
    implicit none
    integer , intent(in) :: nmodes,nq
    real(kind=dp),intent(in) :: gamma
    real(kind=dp),intent(in) :: ld_fric
    real(kind=dp),intent(in) :: wf(nmodes,nq)
    real(kind=dp),intent(out):: ld_gamma(nmodes,nq) 
    
    if(abs(gamma) > eps10) then 
      !gamma /= 0.0
      ld_gamma = gamma
    else
      ld_gamma = ld_fric * wf
    endif
    
    
  end subroutine
  
  subroutine get_dEa_dQ(nmodes,nq,nband,nk,P_nk,gmnvkq,isurface,dEa_dQ)
    implicit none
    integer,intent(in) :: nmodes,nq,nband,nk
    integer,intent(in) :: isurface
    real(kind=dp),intent(in)  :: P_nk(nband,nk,nband*nk)
    real(kind=dp),intent(in)  :: gmnvkq(nband,nband,nmodes,nk,nq)
    real(kind=dp),intent(out) :: dEa_dQ(nmodes,nq)
    
    integer :: iq,imode,ik,ikq,iband1,iband2
    
    dEa_dQ = 0.0
    do iq=1,nq
      do imode=1,nmodes
        do ik=1,nk
          ikq = kqmap(ik,iq)
          do iband1=1,nband
            do iband2=1,nband
              dEa_dQ(imode,iq) = dEa_dQ(imode,iq) + &
              gmnvkq(iband1,iband2,imode,ik,iq)*P_nk(iband1,ik,isurface)*P_nk(iband2,ikq,isurface)                 
            enddo
          enddo
        enddo
      enddo
    enddo
    
  end subroutine get_dEa_dQ
  
  !=======================================================================!
  != rk4 method to obtain coordinate and velocitie after a time interval =!
  !=======================================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods               =!
  !=======================================================================!
  subroutine rk4_nuclei(nmodes,nq,dEa_dQ,ld_gamma,wf,xx,vv,tt)
    implicit none
    integer , intent(in) :: nmodes,nq
    real(kind=dp),intent(in)   :: ld_gamma(nmodes,nq)
    real(kind=dp),intent(in)   :: dEa_dQ(nmodes,nq)
    real(kind=dp),intent(in)   :: wf(nmodes,nq)
    real(kind=dp),intent(inout):: xx(nmodes,nq),vv(nmodes,nq)
    real(kind=dp),intent(in)   :: tt
    real(kind=dp):: tt2,tt6
    real(kind=dp):: xx0(nmodes,nq),dx1(nmodes,nq),dx2(nmodes,nq),dx3(nmodes,nq),dx4(nmodes,nq)
    real(kind=dp):: vv0(nmodes,nq),dv1(nmodes,nq),dv2(nmodes,nq),dv3(nmodes,nq),dv4(nmodes,nq)

    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx,vv,dx1,dv1)
    xx0=xx+tt2*dx1; vv0=vv+tt2*dv1
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx0,vv0,dx2,dv2)
    xx0=xx+tt2*dx2; vv0=vv+tt2*dv2
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx0,vv0,dx3,dv3)
    xx0=xx+tt*dx3; vv0=vv+tt*dv3
    call derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx0,vv0,dx4,dv4)
    xx=xx+tt6*(dx1+2.0d0*dx2+2.0d0*dx3+dx4)
    vv=vv+tt6*(dv1+2.0d0*dv2+2.0d0*dv3+dv4)
  endsubroutine
  
  !====================================================!
  != calculate derivative of coordinate and velocitie =!
  !====================================================!
  != ref: notebook page 630                           =!
  !====================================================!
  ! ref : Huangkun <固体物理> (3-10)
  ! ref : 1 D. M. F. M. Germana Paterlini, Chemical Physics 236 (1998) 243.
  ! ref : PPT-92
  ! ref : 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
  subroutine derivs_nuclei(nmodes,nq,dEa_dQ,wf,ld_gamma,xx,vv,dx,dv)
    implicit none
    integer,intent(in) :: nmodes,nq
    real(kind=dp),intent(in)  ::  dEa_dQ(nmodes,nq),wf(nmodes,nq)
    real(kind=dp),intent(in)  ::  ld_gamma(nmodes,nq)
    real(kind=dp),intent(in)  ::  xx(nmodes,nq),vv(nmodes,nq)
    real(kind=dp),intent(out) ::  dx(nmodes,nq),dv(nmodes,nq)
    
    integer :: iq,imode,ik,iband1,iband2,ikq
    
    dx = vv
    dv = -wf**2*xx - dEa_dQ - ld_gamma*vv
    
    !do iq=1,nq
    !  do imode=1,nmodes
    !    dx(imode,iq) = vv(imode,iq)
    !    dv(imode,iq) = -wf(imode,iq)**2*xx(imode,iq)-dEa_dQ(imode,iq)-ld_gamma(imode,iq)*vv(imode,iq)
    !  enddo
    !enddo
    
  endsubroutine derivs_nuclei
  

  
  subroutine derivs_electron_diabatic(nfre,HH,cc,dc)
    use f95_precision
    use blas95
    use constants,only : cmplx_i,cmplx_0
    implicit none
    integer,intent(in)           :: nfre
    real(kind=dp),intent(in)     :: HH(nfre,nfre)
    complex(kind=dpc),intent(in) :: cc(nfre)
    complex(kind=dpc),intent(inout):: dc(nfre)
    
    dc= cmplx_0
    
    !dc = MATMUL(HH,c)
    call gemv(HH,cc,dc)
  
    dc = dc *(-cmplx_i)

  endsubroutine derivs_electron_diabatic

  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  subroutine rk4_electron_diabatic(nfre,HH,cc,cc0,dc1,dc2,dc3,dc4,tt)
    implicit none
    integer,intent(in)              :: nfre
    complex(kind=dpc),intent(inout) :: cc(nfre)
    real(kind=dp),intent(in)        :: HH(nfre)
    complex(kind=dpc),intent(out)   :: cc0(nfre),dc1(nfre),&
                        dc2(nfre),dc3(nfre),dc4(nfre)
    real(kind=dp),intent(in)        :: tt
    real(kind=dp):: tt2,tt6
    real(kind=dp):: sum_cc2
    
    tt2=tt/2.0d0; tt6=tt/6.0d0
    
    call derivs_electron_diabatic(nfre,HH,cc,dc1)
    cc0=cc+tt2*dc1
    !sum_cc2 = REAL(SUM(cc0*CONJG(cc0)))
    !cc0= cc0/sqrt(sum_cc2)
    call derivs_electron_diabatic(nfre,HH,cc0,dc2)
    cc0=cc+tt2*dc2
    !sum_cc2 = REAL(SUM(cc0*CONJG(cc0)))
    !cc0= cc0/sqrt(sum_cc2)   
    call derivs_electron_diabatic(nfre,HH,cc0,dc3)
    cc0=cc+tt*dc3     
    !sum_cc2 = REAL(SUM(cc0*CONJG(cc0)))
    !cc0= cc0/sqrt(sum_cc2)
    call derivs_electron_diabatic(nfre,HH,cc0,dc4)
    
    cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)
    
    sum_cc2 = REAL(SUM(cc*CONJG(cc)))
    cc = cc/sqrt(sum_cc2)
    
  endsubroutine rk4_electron_diabatic
      
  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!


  subroutine get_dEa2_dQ2(nmodes,nq,nfre,nfre_sh,isurface,EE,dd,dEa2_dQ2)
    implicit none
    integer,intent(in) :: nmodes,nq,nfre,nfre_sh
    integer,intent(in) :: isurface
    real(kind=dp),intent(in)  :: EE(nfre),dd(nfre_sh,nfre_sh,nmodes,nq)
    real(kind=dp),intent(out) :: dEa2_dQ2(nmodes,nq)
    
    integer :: iq,imode,ifre
    
    dEa2_dQ2 = 0.0
    do iq=1,nq
      do imode=1,nmodes
        do ifre=1,nfre_sh
          if(ifre /= isurface) then
            dEa2_dQ2(imode,iq) = dEa2_dQ2(imode,iq) + (EE(ifre)-EE(isurface))*DD(ifre,isurface,imode,iq)*DD(isurface,ifre,imode,iq)
          endif
        enddo
      enddo
    enddo 
    
  end subroutine get_dEa2_dQ2
  
  
  !===============================================!
  != add bath effect to coordinate and velocitie =!
  !===============================================!
  != ref: notebook page 462 and 638              =!
  !===============================================!
  ! ref: 1 D. M. F. M. Germana Paterlini, Chemical Physics 236 (1998) 243.
  SUBROUTINE ADD_BATH_EFFECT(nmodes,nq,wf,ld_gamma,temp,dEa2_dQ2,TT,l_ph_quantum,XX,VV)
    use kinds,only : dp,dpc
    use randoms,only : GAUSSIAN_RANDOM_NUMBER_FAST
    use constants,only : KB=>K_B_Ryd,sqrt3,sqrt5,sqrt7
    implicit none
    
    integer , intent(in)      :: nq,nmodes
    real(kind=dp), intent(in) :: wf(nmodes,nq)
    real(kind=dp), intent(in) :: ld_gamma(nmodes,nq)
    real(kind=dp), intent(in) :: temp
    real(kind=dp), intent(in) :: dEa2_dQ2(nmodes,nq)
    real(kind=dp), intent(in) :: tt
    logical,intent(in)        :: l_ph_quantum
    real(kind=dp), intent(inout) :: XX(nmodes,nq),VV(nmodes,nq)
    
    integer :: imode,iq
    real(kind=dp) :: SIGMAR,R1,R2,R3,R4,Z1,Z2,Z3,Z4
    real(kind=dp) :: wwf2
    real(kind=dp) :: gamma,womiga,aver_E_T
    
    DO iq=1,nq
      do imode=1,nmodes
        womiga = wf(imode,iq)
        gamma = ld_gamma(imode,iq)
        if(womiga > 0.0) then
					if(l_ph_quantum) then
						aver_E_T = (bolziman(womiga,temp)+0.5)*womiga
						!aver_E_T = KB*TEMP
					else
						aver_E_T = KB*TEMP
					endif
				else
					aver_E_T = 0.0
				endif
        SIGMAR=DSQRT(2.0*gamma*aver_E_T*TT)
        wwf2 = wf(imode,iq)**2+dEa2_dQ2(imode,iq)
      
        R1=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        R2=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        R3=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        R4=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        Z1=R1    ! V=phP
        Z2=TT*(R1/2.0D0+R2/SQRT3/2.0D0)  !V*T
        Z3=TT**2*(R1/6.0D0+R2*SQRT3/12.0D0+R3/SQRT5/12.0D0) ! V*T**2
        Z4=TT**3*(R1/24.0D0+R2*SQRT3/40.0D0+R3/SQRT5/24.0D0+R4/SQRT7/120.0D0) ! V*T**3
        XX(imode,iq)=XX(imode,iq)+(Z2-GAMMA*Z3+(-wwf2+GAMMA**2)*Z4)
        VV(imode,iq)=VV(imode,iq)+(Z1-GAMMA*Z2+(-wwf2+GAMMA**2)*Z3+(2.0*gamma*wwf2-GAMMA**3)*Z4)
        
      enddo
    enddo
    
  ENDSUBROUTINE  
  
  !ref : 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !    : (3.24)  
  function bolziman(womiga,temp)
    use io,only :stdout
    use constants,only : K_B_Ryd
    implicit none
    real(kind=dp),intent(in)::womiga,temp
    real(kind=dp) :: bolziman
    if(womiga == 0) then
      write(stdout,*) "womiga == 0,phonon error"
      stop
    endif
    bolziman=1.0/(exp(womiga/(K_B_Ryd*temp))-1.0)
    !<nb>=1/(exp{hbar*w/kbT}-1)
  end function     
  
end module dynamics