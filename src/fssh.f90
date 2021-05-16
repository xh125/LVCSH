module fssh
  use kinds,only : dp,dpc
  implicit none


  contains
  subroutine nonadiabatic_transition_fssh(nfre,nq,nmodes,isurface,EE,P,DD,GG,VV)
    use randoms,only :more_random
    implicit none
    integer , intent(in)     :: nfre,nq,nmodes
    integer , intent(inout)  :: isurface
    real(kind=dp),intent(in) :: EE(nfre)
    real(kind=dp),intent(in) :: P(nfre,nfre)
    real(kind=dp),intent(in) :: DD(nfre,nfre,nmodes,nq)
    real(kind=dp),intent(in) :: GG(nfre)
    real(kind=dp),intent(inout) :: VV(nmodes,nq)
    
    real(kind=dp) :: sumvd,sumdd,sumgg,flagr,flagd,detaE
    integer :: ifre,jfre,imode,iq
    !real(kind=dp) :: SUM_S,max_Sbi(1),SUM_G
    
    call more_random()
    call random_number(flagr)
    sumgg = 0.0d0
    do ifre=1,nfre
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
          !if(.not. lelec) detaE = -1.0*detaE  ! 针对空穴修改能量
          flagd = 1.0+2.0*detaE*sumdd/sumvd**2  
          
          if(flagd > 0.0) then
            flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
            do iq=1,nq
              do imode=1,nmodes
                VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,ifre,imode,iq)
              enddo
            enddo
            isurface = ifre          
          endif
          exit
        endif
      endif
    enddo
  
  end subroutine nonadiabatic_transition_fssh
  
end module fssh