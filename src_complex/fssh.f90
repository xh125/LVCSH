module fssh
  use kinds,only : dp,dpc
  implicit none

  contains
  subroutine nonadiabatic_transition_fssh(lfeedback,nfre,nfre_sh,nq,nmodes,isurface,EE,P,DD,GG,VV)
    use randoms,only :more_random
    implicit none
    logical , intent(in)     :: lfeedback
    integer , intent(in)     :: nfre,nfre_sh,nq,nmodes
    integer , intent(inout)  :: isurface
    real(kind=dp),intent(in) :: EE(nfre)
    complex(kind=dp),intent(in) :: P(nfre,nfre)
    complex(kind=dp),intent(in) :: DD(nfre_sh,nfre_sh,nmodes,nq)
    real(kind=dp),intent(in) :: GG(nfre_sh)
    real(kind=dp),intent(inout) :: VV(nmodes,nq)
    
		complex(kind=dp) :: sumvd,sumdd
    real(kind=dp) :: sumgg,flagr,flagd,detaE
    integer :: ifre,jfre,imode,iq
    real(kind=dp) :: SUM_E0,SUM_E1,dSUM_E
    
    SUM_E0 = 0.5*SUM(VV**2)+EE(isurface)
    
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
              sumvd = sumvd + VV(imode,iq)*DD(isurface,ifre,imode,iq) ! A
              sumdd = sumdd + DD(isurface,ifre,imode,iq)**2           ! B
            enddo
          enddo
          detaE = EE(isurface)-EE(ifre)
          flagd = 1.0+2.0*detaE*sumdd/(sumvd**2)  
          
          if(flagd > 0.0) then
            flagd = (sumvd/sumdd)*(-1.0+dsqrt(flagd))
            do iq=1,nq
              do imode=1,nmodes
                if(lfeedback) VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,ifre,imode,iq)
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