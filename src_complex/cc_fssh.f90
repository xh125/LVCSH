module cc_fssh
  use kinds,only : dp,dpc
  use constants,only : cone,czero
  implicit none
  
  integer :: iesurface_j,ihsurface_j
  real(kind=dp),allocatable :: S_ai_e(:),S_ai_h(:)
  real(kind=dp),allocatable :: S_bi_e(:),S_bi_h(:)
 
  contains
  
  subroutine allocate_ccfssh(lelecsh,lholesh,nefre,nefre_sh,nhfre,nhfre_sh)
    implicit none
    logical,intent(in) :: lelecsh,lholesh
    integer,intent(in) :: nefre,nhfre,nefre_sh,nhfre_sh
    
    allocate(S_ai_e(nefre_sh))
    allocate(S_ai_h(nhfre_sh))
    allocate(S_bi_e(nefre_sh))
    allocate(S_bi_h(nhfre_sh))
  
  end subroutine allocate_ccfssh
    
  subroutine get_G_CC_FSSH(nfre,nfre_sh,isurface,isurface_j,p0,p,w0,w,S_ai,g1,g)
    implicit none
    integer,intent(in) :: nfre,nfre_sh,isurface
    integer,intent(out):: isurface_j
    complex(kind=dpc),intent(in) :: p0(nfre,nfre)
    complex(kind=dpc),intent(in) :: p(nfre,nfre)
    complex(kind=dpc),intent(in) :: w0(nfre),w(nfre)
    real(kind=dp),intent(inout)  :: S_ai(nfre_sh)
    real(kind=dp),intent(in)  	 :: g1(nfre_sh)
    real(kind=dp),intent(inout)  :: g(nfre_sh)
    
    integer :: ifre
    integer :: isurface_a,isurface_b,isurface_k
    real(kind=dp) :: S_aa,SUM_S
    integer :: max_Sai(1)
    real(kind=dp) :: sumg0,sumg1
    
    ! ref : 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
    !lallocate = allocated(S_ai)
    !if(.not. lallocate) allocate(S_ai(nfre))
    S_ai = 0.0
    isurface_a = isurface
    S_aa = SUM(p0(:,isurface_a)*p(:,isurface_a))
    
    !S_aa = SUM(CONJG(p0(:,isurface_a))*p(:,isurface_a))
    if(S_aa**2 > 0.5) then
      !type(1) or type(3)
      isurface_j = isurface_a
    else
      !type(2) or type(4)
      do ifre = 1,nfre_sh
        S_ai(ifre) = ABS(SUM(p0(:,isurface_a)*p(:,ifre)))**2
        !S_ai(ibasis) = SUM(CONJG(pp0(:,isurface_a))*pp(:,ibasis))
      enddo
      !S_ai = S_ai**2
      SUM_S = SUM(S_ai)   ! = 1.00
      max_Sai =  MAXLOC(S_ai)
      isurface_j = max_Sai(1)
    endif
    
    sumg0 = (ABS(W0(ISURFACE))**2-ABS(W(ISURFACE))**2)/ABS(W0(ISURFACE))**2
    sumg1 = SUM(g1)      
    if(isurface_j == isurface) then
      ! type(1) or type(3)
      g = g
    else
      g(isurface_j) = sumg0 - (SUM(G1)-G1(isurface_j))
      if(G(isurface_j) < 0.0d0) G(isurface_j) = 0.0d0
      if(SUM(G) > 1.0d0) G=G/SUM(G)
    endif
    !    
  end subroutine get_G_CC_FSSH  
  
  subroutine nonadiabatic_transition_ccfssh(lfeedback,nfre,nfre_sh,nmodes,nq,&
              &isurface,isurface_j,surface_type,EE,P,P0,DD,S_bi,GG,VV)
    use randoms,only :more_random
    use constants,only : ryd2ev,eps10
    use io,only : stdout 
    implicit none
    logical , intent(in)     :: lfeedback
    integer , intent(in)     :: nfre,nfre_sh,nmodes,nq
    integer , intent(inout)  :: isurface
    integer , intent(in)     :: isurface_j
    integer , intent(out)    :: surface_type
    real(kind=dp),intent(in) :: EE(nfre)
    complex(kind=dpc),intent(in) :: P(nfre,nfre),P0(nfre,nfre)
    complex(kind=dpc),intent(in) :: DD(2,nfre_sh,nmodes,nq)
    real(kind=dp),intent(out):: S_bi(nfre_sh)
    real(kind=dp),intent(in) :: GG(nfre_sh)
    complex(kind=dpc),intent(inout) :: VV(nmodes,nq)
    
    real(kind=dp) :: sumvd
    real(kind=dp) :: sumdd,sumgg,flagr,flagd,detaE
    integer :: ifre,jfre,imode,iq,isurface_a,isurface_b,isurface_k
    real(kind=dp) :: max_Sbi(1)
    real(kind=dp) :: SUM_E0,SUM_E1,dSUM_E
    
    SUM_E0 = 0.5*SUM(VV*CONJG(VV))+EE(isurface)
    
    isurface_a = isurface
    call more_random()
    call random_number(flagr)
    sumgg = 0.0d0
    do ifre=1,nfre_sh
      if(ifre /= isurface_a) then
        sumgg = sumgg + GG(ifre)
        if(flagr < sumgg) then
          isurface_b = ifre
          
          S_bi = 0.0
          do jfre = 1,nfre_sh
            S_bi(jfre) = ABS(SUM(p0(:,isurface_b)*CONJG(p(:,jfre))))**2
          enddo
          max_Sbi =  MAXLOC(S_bi)
          isurface_k = max_Sbi(1)
          
          if(isurface_j == isurface_a) then   ! type (1) and (3)
            if(isurface_k == isurface_b) then
              surface_type = 1
            else
              surface_type = 3
            endif
          else !(j/=a) type (2) and (4)
            if(isurface_k == isurface_b .or. isurface_b == isurface_j ) then
              surface_type = 3
            else!(isurface_k /= isurface_b .and. isurface_b /= isurface_j )
              surface_type = 4
            endif
          endif

          
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
          detaE = EE(isurface_a)-EE(ifre) ! C
          ! 2.0*detaE = 2.0*sumvd*dt+sumdd*(dt2)

          flagd = 1.0+2.0*detaE*sumdd/sumvd**2  
          
          if(isurface_j == isurface_a) then ! type (1) (3)
            if(flagd > 0.0) then
              flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
              do iq=1,nq
                do imode=1,nmodes
                  if(lfeedback) VV(imode,iq) = VV(imode,iq) + flagd*dd(1,ifre,imode,iq)
                enddo
              enddo
              isurface = isurface_k
              
              !
              SUM_E1 = 0.5*SUM(VV*CONJG(VV))+EE(isurface_b)
              dSUM_E = SUM_E1 - SUM_E0
              if(lfeedback) then
                if(ABS(dSUM_E)>eps10) then
                  write(stdout,"(5X,A49,F12.4,A4)") "In feedback, the energy is not conservied,dSUM_E=",dSUM_E*ryd2ev," eV."
                endif
              endif              
              !
              
            else
              isurface = isurface_a
            endif
          elseif(isurface_j /= isurface_a) then !type (2)(4)
            if(isurface_b /= isurface_j) then!b/=j
              if(flagd > 0.0) then
                flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
                do iq=1,nq
                  do imode=1,nmodes
                    if(lfeedback) VV(imode,iq) = VV(imode,iq) + flagd*dd(1,ifre,imode,iq)
                  enddo
                enddo
                isurface = isurface_k
                
                !
                SUM_E1 = 0.5*SUM(VV*CONJG(VV))+EE(isurface_b)
                dSUM_E = SUM_E1 - SUM_E0
                if(lfeedback) then
                  if(ABS(dSUM_E)>eps10) then
                    write(stdout,"(5X,A49,F12.4,A4)") "In feedback, the energy is not conservied,dSUM_E=",dSUM_E*ryd2ev," eV."
                  endif
                endif
                !
                
              else
                isurface = isurface_j
              endif
            else  !b=j
              isurface = isurface_j
            endif
          
          endif
          
          exit
        endif
      endif
    enddo
  
    
  end subroutine nonadiabatic_transition_ccfssh
  
end module cc_fssh