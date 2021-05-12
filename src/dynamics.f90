module dynamics
  use kinds,only : dp,dpc
  use io,only :
  use parameters,only : gamma,temp
  use hamiltonian,only : nphfre,nefre,set_H_nk
  use elph2,only          : nk=>nktotf,nq=>nqtotf,wf,nband=>nbndfst,epcq
  use modes,only          : nmodes
  use epwcom,only         : kqmap
  use surfacehopping,only : iesurface,ihsurface
  
  implicit none
  logical :: lelec = .true.
  logical :: lhole = .false.
  
  contains
  
  !=======================================================================!
  != rk4 method to obtain coordinate and velocitie after a time interval =!
  !=======================================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods               =!
  !=======================================================================!
  subroutine rk4_nuclei(pp_nk,xx,vv,tt)
    implicit none
    real(kind=dp),intent(in)   :: pp_nk(nband,nk,nefre)
    real(kind=dp),intent(inout):: xx(nmodes,nq),vv(nmodes,nq)
    real(kind=dp),intent(in)   :: tt
    real(kind=dp):: tt2,tt6
    real(kind=dp):: xx0(nmodes,nq),dx1(nmodes,nq),dx2(nmodes,nq),dx3(nmodes,nq),dx4(nmodes,nq)
    real(kind=dp):: vv0(nmodes,nq),dv1(nmodes,nq),dv2(nmodes,nq),dv3(nmodes,nq),dv4(nmodes,nq)

    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_nuclei(pp_nk,xx,vv,dx1,dv1)
    xx0=xx+tt2*dx1; vv0=vv+tt2*dv1
    call derivs_nuclei(pp_nk,xx0,vv0,dx2,dv2)
    xx0=xx+tt2*dx2; vv0=vv+tt2*dv2
    call derivs_nuclei(pp_nk,xx0,vv0,dx3,dv3)
    xx0=xx+tt*dx3; vv0=vv+tt*dv3
    call derivs_nuclei(pp_nk,xx0,vv0,dx4,dv4)
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
  subroutine derivs_nuclei(pp_nk,xx,vv,dx,dv)
    implicit none
    real(kind=dp),intent(in)   :: pp_nk(nband,nk,nefre)
    real(kind=dp),intent(in)  ::  xx(nmodes,nq),vv(nmodes,nq)
    real(kind=dp),intent(out) ::  dx(nmodes,nq),dv(nmodes,nq)
    integer :: iq,imode,ik,iband1,iband2,ikq
    real(kind=dp)::dEa_dQ
    do iq=1,nq
      do imode=1,nmodes
        dx(imode,iq) = vv(imode,iq)
        dv(imode,iq) = (-wf(imode,iq)**2*xx(imode,iq)-gamma*vv(imode,iq))
        
        dEa_dQ = 0.0
        
        do ik=1,nk
          ikq = kqmap(ik,iq)
          do iband1=1,nband
            do iband2=1,nband
              ! 电子能量随简正模的变化
              dEa_dQ = dEa_dQ + &
              pp_nk(iband1,ik,iesurface)*pp_nk(iband2,ikq,iesurface)*epcq(iband1,iband2,ik,imode,iq)
              ! 空穴能量随简正模的变化
              dEa_dQ = dEa_dQ - &
              pp_nk(iband1,ik,ihsurface)*pp_nk(iband2,ikq,iesurface)*epcq(iband1,iband2,ik,imode,iq)                    
            enddo
          enddo
        enddo
        dEa_dQ = sqrt(2.0*wf(imode,iq)/nq) * dEa_dQ
        dv(imode,iq) = dv(imode,iq) - dEa_dQ
      enddo
    enddo
    
  endsubroutine derivs_nuclei
  
  subroutine derivs_electron_diabatic(xx,c_nk,dc_nk,llelec)
    use f95_precision
    use blas95
    use constants,only : cmplx_i,cmplx_0
    implicit none
    real(kind=dp),intent(in)     :: xx(nmodes,nq)
    complex(kind=dpc),intent(in) :: c_nk(nband,nk)
    complex(kind=dpc),intent(inout):: dc_nk(nband,nk)
    logical,intent(in) :: llelec
    
    complex(kind=dpc) :: c(nefre),dc(nefre)
    real(kind=dp) :: HH(nefre,nefre)
    call set_H_nk(xx)
    
    c = reshape(c_nk,(/nefre/))
    dc= cmplx_0
    
    !dc = MATMUL(HH,c)
    call gemv(HH,c,dc)
    
    if(llelec) then
      dc = dc *(-cmplx_i)
    else
      dc = dc *(cmplx_i)
    endif
    
    dc_nk = reshape(dc,(/nband,nk/))
    
  endsubroutine derivs_electron_diabatic

  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  subroutine rk4_electron_diabatic(xx,cc,tt,llelec)
    implicit none
    real(kind=dp),intent(in)        :: xx(nmodes,nq)
    complex(kind=dpc),intent(inout) :: cc(nband,nk)
    real(kind=dp),intent(in)        :: tt
    logical,intent(in)              :: llelec
    real(kind=dp):: tt2,tt6
    complex(kind=dpc):: cc0(nband,nk),dc1(nband,nk),&
                        dc2(nband,nk),dc3(nband,nk),dc4(nband,nk)
    real(kind=dp)    :: nn(nband,nk) ,sum_nn
    
    nn=CONJG(cc)*cc
    sum_nn = SUM(nn)
    
    tt2=tt/2.0d0; tt6=tt/6.0d0
    
    call derivs_electron_diabatic(xx,cc,dc1,llelec)
    cc0=cc+tt2*dc1
    nn=CONJG(cc0)*cc0
    sum_nn = SUM(nn)    
    call derivs_electron_diabatic(xx,cc0,dc2,llelec)
    cc0=cc+tt2*dc2
    nn=CONJG(cc0)*cc0
    sum_nn = SUM(nn)       
    call derivs_electron_diabatic(xx,cc0,dc3,llelec)
    cc0=cc+tt*dc3
    nn=CONJG(cc0)*cc0
    sum_nn = SUM(nn)       
    call derivs_electron_diabatic(xx,cc0,dc4,llelec)
    cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)
    nn=CONJG(cc)*cc
    sum_nn = SUM(nn)
    
    !cc = cc/dsqrt(sum_nn)
    !nn=CONJG(cc)*cc    
  endsubroutine rk4_electron_diabatic  
  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  
  !===============================================!
  != add bath effect to coordinate and velocitie =!
  !===============================================!
  != ref: notebook page 462 and 638              =!
  !===============================================!
  ! ref: 1 D. M. F. M. Germana Paterlini, Chemical Physics 236 (1998) 243.
  SUBROUTINE ADD_BATH_EFFECT(nfre,EE,E0,PP,DD,TT,XX,VV)
    use kinds,only : dp,dpc
    use randoms,only : GAUSSIAN_RANDOM_NUMBER_FAST
    use parameters,only : gamma,temp
    use constants,only : KB=>K_B_Ryd,sqrt3,sqrt5,sqrt7
    use surfacehopping,only : iesurface,ihsurface,SUM_ph_U,SUM_ph_T,SUM_ph_E,&
                              ph_T,ph_U
    implicit none
    
    integer , intent(in)      :: nfre
    real(kind=dp), intent(in) :: EE(nfre),E0(nfre)
    real(kind=dp), intent(in) :: PP(nfre,nfre)
    real(kind=dp), intent(in) :: DD(nfre,nfre,nmodes,nq)
    real(kind=dp), intent(in) :: tt
    real(kind=dp), intent(inout) :: XX(nmodes,nq),VV(nmodes,nq)
    
    integer :: imode,iq,ifre
    real(kind=dp) :: SIGMAR,R1,R2,R3,R4,Z1,Z2,Z3,Z4!,GAUSSIAN_RANDOM_NUMBER_FAST
    !EXTERNAL GAUSSIAN_RANDOM_NUMBER_FAST
    real(kind=dp) :: wwf2
    real(kind=dp) :: dEa2_dQ2
    
    SIGMAR=DSQRT(2.0*gamma*KB*TEMP*TT)
    DO iq=1,nq
      do imode=1,nmodes
        dEa2_dQ2 = 0.0
        !do ifre=1,nfre
        !  if(iefre /= iesurface) then
        !    dEa2_dQ2 = dEa2_dQ2 + (EE(iefre)-EE(iesurface))*DD(iefre,iesurface,imode,iq)*DD(iesurface,iefre,imode,iq)
        !  endif
        !  if(iefre /= ihsurface) then
        !    dEa2_dQ2 = dEa2_dQ2 + (ee(iefre)-ee(ihsurface))*dd(iefre,ihsurface,imode,iq)*dd(ihsurface,iefre,imode,iq)
        !  endif
        !enddo
        wwf2 = wf(imode,iq)**2+dEa2_dQ2
      
   
        !KK=K
        !DO JSITE=1,NSITE
        !  IF(JSITE.NE.ISURFACE) KK=KK+2.0D0*ALPHA*DD(JSITE,ISURFACE,ISITE)*PP(ISITE,ISURFACE)*PP(ISITE,JSITE)
        !ENDDO
    
        R1=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        R2=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        R3=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        R4=GAUSSIAN_RANDOM_NUMBER_FAST(0.0D0,SIGMAR)
        Z1=R1    ! V
        Z2=TT*(R1/2.0D0+R2/SQRT3/2.0D0)  !V*T
        Z3=TT**2*(R1/6.0D0+R2*SQRT3/12.0D0+R3/SQRT5/12.0D0) ! V*T**2
        Z4=TT**3*(R1/24.0D0+R2*SQRT3/40.0D0+R3/SQRT5/24.0D0+R4/SQRT7/120.0D0) ! V*T**3
        XX(imode,iq)=XX(imode,iq)+(Z2-GAMMA*Z3+(-wwf2+GAMMA**2)*Z4)
        VV(imode,iq)=VV(imode,iq)+(Z1-GAMMA*Z2+(-wwf2+GAMMA**2)*Z3+(2.0*gamma*wwf2-GAMMA**3)*Z4)
        
      enddo
    enddo
    ph_T = 0.5*VV**2
    ph_U = 0.5*(wf**2)*(XX**2)
    SUM_ph_U=SUM(ph_U)
    SUM_ph_T=SUM(ph_T)
    SUM_ph_E=SUM_ph_T+SUM_ph_U
    
  ENDSUBROUTINE  
  
  
end module dynamics