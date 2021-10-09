module sc_fssh
  use kinds,only : dp,dpc
  use constants,only : cone,czero
  implicit none

  real(kind=dp) :: minde_e,minde_h
  
  contains
  subroutine get_G_SC_FSSH(isurface,nfre,nfre_sh,e0,w0,w,g1,g)
    use constants ,only   : ryd2eV
    implicit none
    integer,intent(in) :: isurface
    integer,intent(in) :: nfre,nfre_sh
    real(kind=dp),intent(in)     :: e0(nfre)
    complex(kind=dpc),intent(in) :: w0(nfre),w(nfre)
    real(kind=dp),intent(in)     :: g1(nfre_sh)
    real(kind=dp),intent(inout)  :: g(nfre_sh)
    
    integer :: ifre
    real(kind=dp) :: sumg0,sumg1
    real(kind=dp) :: minde
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% calculate sumg0,sumg1,minde %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    ! sumg0 总的跃迁几率(透热表象下计算得出)
    ! ref: 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
    ! eq(4) ,eq(5)
    ! ref: 1 L. Wang, and O. V. Prezhdo, Journal of Physical Chemistry Letters 5 (2014) 713.
    ! eq(4)-eq(12)
    ! Assumption at most one trivial crossing is encountered during a time step.    
    
    sumg0=(abs(w0(isurface))**2-abs(w(isurface))**2)/abs(w0(isurface))**2
    sumg1=sum(g1)
    if(isurface == 1) then
      minde=(e0(isurface+1)-e0(isurface))*ryd2eV
    elseif(isurface == nfre) then
      minde=(e0(isurface)-e0(isurface-1))*ryd2eV
    elseif((e0(isurface+1)-e0(isurface)) < (e0(isurface)-e0(isurface-1))) then
      minde=(e0(isurface+1)-e0(isurface))*ryd2eV
    else
      minde=(e0(isurface)-e0(isurface-1))*ryd2eV
    endif    

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% fixed g                         %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    if(isurface == 1) then
      ifre=isurface+1
    elseif(isurface == nfre_sh) then
      ifre=isurface-1
    elseif((e0(isurface+1)-e0(isurface)) < (e0(isurface)-e0(isurface-1))) then
      ifre=isurface+1
    else
      ifre=isurface-1
    endif
    g(ifre)=sumg0-(sum(g1)-g1(ifre))
    if(g(ifre).lt.0.0d0) g(ifre)=0.0d0
    if(sum(g).ge.1.0d0) g=g/sum(g)    
    
  end subroutine get_G_SC_FSSH  
  
  subroutine nonadiabatic_transition_scfssh(lfeedback,nfre,nfre_sh,nq,nmodes,isurface,EE,P,DD,GG,VV)
    use randoms,only :more_random
    implicit none
    logical , intent(in)     :: lfeedback
    integer , intent(in)     :: nfre,nfre_sh,nq,nmodes
    integer , intent(inout)  :: isurface
    real(kind=dp),intent(in) :: EE(nfre)
    complex(kind=dpc),intent(in) :: P(nfre,nfre)
    complex(kind=dpc),intent(in) :: DD(2,nfre_sh,nmodes,nq)
    real(kind=dp),intent(in) :: GG(nfre_sh)
    complex(kind=dpc),intent(inout) :: VV(nmodes,nq)
    
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
          detaE = EE(isurface)-EE(ifre) ! C
          ! 2.0*detaE = 2.0*sumvd*dt+sumdd*(dt2)

          flagd = 1.0+2.0*detaE*sumdd/(sumvd**2)  
          
          if(flagd > 0.0) then
            flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
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
    dSUM_E = SUM_E1 - SUM_E0
    
  end subroutine nonadiabatic_transition_scfssh
  
end module sc_fssh