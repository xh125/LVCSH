module fssh
  use kinds,only : dp,dpc
  !use constants,only : cone,czero

  implicit none

  contains
  subroutine nonadiabatic_transition_fssh(lfeedback,nfre,nfre_sh,nq,nmodes,isurface,EE,P,DD,GG,VV)
    use elph2,only : iminusq
    use randoms,only : more_random
    use dynamics,only : test_v
    implicit none
    logical , intent(in)     :: lfeedback
    integer , intent(in)     :: nfre,nfre_sh,nq,nmodes
    integer , intent(inout)  :: isurface
    real(kind=dp),intent(in) :: EE(nfre)
    complex(kind=dpc),intent(in) :: P(nfre,nfre)
    complex(kind=dpc),intent(in) :: DD(2,nfre_sh,nmodes,nq)
    real(kind=dp),intent(in) :: GG(nfre_sh)
    complex(kind=dp),intent(inout) :: VV(nmodes,nq)
    
    real(kind=dp) :: sumvd
    real(kind=dp) :: sumdd,sumgg,flagr,flagd,detaE
    integer :: ifre,jfre,imode,iq,iq_
    real(kind=dp) :: SUM_E0,SUM_E1,dSUM_E
    complex(kind=dpc) :: DDq_q
    
    SUM_E0 = 0.5*SUM(VV*CONJG(VV))+EE(isurface)
    
    call more_random()
    call random_number(flagr)
    sumgg = 0.0d0
    do ifre=1,nfre_sh
      if(ifre /= isurface) then
        sumgg = sumgg + GG(ifre)
        if(flagr < sumgg) then
          ! ref:1 L. Wang, and D. Beljonne, Journal of Physical Chemistry Letters 4 (2013) 1888.
          ! eq(16)
          sumvd = 0.0
          sumdd = 0.0
          do iq=1,nq
            iq_ = iminusq(iq)
            do imode=1,nmodes
              !DDq_q = DD(1,ifre,imode,iq)
              DDq_q = (DD(1,ifre,imode,iq)-DD(2,ifre,imode,iq))/2.0
              sumvd = sumvd + REAL(CONJG(VV(imode,iq))*DDq_q) ! B
              sumdd = sumdd + ABS(DDq_q)**2  ! 2A                   
            enddo
          enddo
          detaE = EE(isurface)-EE(ifre) !C
          ! 2.0*detaE = 2.0*sumvd*dt+sumdd*(dt2)
          ! ref: 1 H. Xie et al., J. Chem. Phys. 156 (2022) 154116.
          ! Eq: (29)
          ! A*dt2 + B*dt - C = 0          
            
          ! b2+4ac
          flagd = sumvd**2 + 2.0*sumdd*detaE
          !flagd = 1.0+2.0*detaE*sumdd/(sumvd**2)  
          
          
          if(flagd > 0.0) then
            if(sumvd > 0.0) then
              flagd = (-1.0*sumvd+dsqrt(flagd))/sumdd  ! effective scratting time
            else
              flagd = (-1.0*sumvd-dsqrt(flagd))/sumdd
            endif
            if(lfeedback) then
              do iq=1,nq
                do imode=1,nmodes
                  !DDq_q = DD(1,ifre,imode,iq)
                  DDq_q = (DD(1,ifre,imode,iq)-DD(2,ifre,imode,iq))/2.0                
                  VV(imode,iq) = VV(imode,iq) + flagd*DDq_q
                enddo
              enddo
              
              ! Let VV CONJG again
              do iq=1,nq
                iq_=iminusq(iq)
                if(iq <= iq_) then
                do imode = 1, nmodes
                  if (iq == iq_) then
                    VV(imode,iq) = SIGN(ABS(VV(imode,iq)),REAL(VV(imode,iq)))
                  else
                    VV(imode,iq) = (VV(imode,iq)+CONJG(VV(imode,iq_)))/2.0
                    VV(imode,iq_)= CONJG(VV(imode,iq))
                  endif
                enddo
                endif
              enddo
              
              
            endif
            isurface = ifre   
            call test_v(nmodes,nq,VV)    
          endif
          
          exit
        
        endif
      endif
    enddo
    
    SUM_E1 = 0.5*SUM(VV**2)+EE(isurface)
    dSUM_E = SUM_E1-SUM_E0
    
  end subroutine nonadiabatic_transition_fssh
  
end module fssh