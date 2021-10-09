module fssh
  use kinds,only : dp,dpc
  !use constants,only : cone,czero

  implicit none

  contains
  subroutine nonadiabatic_transition_fssh(lfeedback,nfre,nfre_sh,nq,nmodes,isurface,EE,P,DD,GG,VV)
    use randoms,only :more_random
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
    integer :: ifre,jfre,imode,iq
    real(kind=dp) :: SUM_E0,SUM_E1,dSUM_E
    
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
            do imode=1,nmodes
              sumvd = sumvd + REAL(CONJG(VV(imode,iq))*DD(1,ifre,imode,iq)) ! A
              sumdd = sumdd + ABS(DD(1,ifre,imode,iq))**2  ! B                   
            enddo
          enddo
          detaE = EE(isurface)-EE(ifre) !C
          ! 2.0*detaE = 2.0*sumvd*dt+sumdd*(dt2)
          
          flagd = 1.0+2.0*detaE*sumdd/(sumvd**2)  
          
          
          if(flagd > 0.0) then
            flagd = (sumvd/sumdd)*(-1.0+dsqrt(flagd))  ! effective scratting time
            do iq=1,nq
              do imode=1,nmodes
                if(lfeedback) VV(imode,iq) = VV(imode,iq) + flagd*dd(1,ifre,imode,iq)
              enddo
            enddo
            isurface = ifre          
          endif
          
          exit
        
        endif
      endif
    enddo
    
    SUM_E1 = 0.5*SUM(VV**2)+EE(isurface)
    dSUM_E = SUM_E1-SUM_E0
    
  end subroutine nonadiabatic_transition_fssh
  
end module fssh