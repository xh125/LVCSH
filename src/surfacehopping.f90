module surfacehopping
  use kinds, only : dp,dpc
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,epcq,nktotf,nbndfst,nqtotf,ibndmin,ibndmax
  use hamiltonian,only : nphfre,nefre
  use parameters, only : nsnap,naver
  use surfacecom
  implicit none
  !integer :: iaver
  !integer :: isnap,istep
  !integer :: iesurface,ihsurface
  !integer :: ierr
  !
  !! phonons normal mode coordinate,and phonons P
  !real(kind=dp),allocatable :: phQ(:,:),phP(:,:),phQ0(:,:),phP0(:,:)
  !real(kind=dp),allocatable :: e(:),p(:,:),p_nk(:,:,:),d(:,:,:,:),ge(:),gh(:)
  !real(kind=dp),allocatable :: e0(:),p0(:,:),d0(:,:,:,:),ge1(:),gh1(:)  
  !real(kind=dp),allocatable :: pes(:,:,:),inf(:,:,:),csit(:,:),wsit(:,:),&
  !                             psit(:,:),xsit(:,:),ksit(:,:) 
  !real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  !complex(kind=dpc),allocatable :: celec_nk(:,:),w_e(:),w0_e(:)
  !complex(kind=dpc),allocatable :: chole_nk(:,:),w_h(:),w0_h(:)
  !real(kind=dp) :: minde_e,minde_h,sumg0_e,sumg0_h,sumg1_e,sumg1_h
  contains
  
  subroutine allocatesh(nmodes)
    implicit none
    integer,intent(in):: nmodes
    allocate(phQ(nmodes,nqtotf),stat=ierr)  ! x(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ',1)
    phQ = 0.0
    ! phonons normal mode coordinates
    allocate(phP(nmodes,nqtotf),stat=ierr)  ! v(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP',1)    
    phP = 0.0
    !phonons normal mode verlosity
    allocate(e(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating e',1)
    allocate(p(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating p',1)
    allocate(p_nk(nbndfst,nktotf,nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating p_nk',1)
    allocate(d(nefre,nefre,nmodes,nqtotf),stat=ierr) !d_ijk
    if(ierr /=0) call errore('surfacehopping','Error allocating d',1)
    allocate(ge(1:nefre),stat=ierr)  !g_ij
    if(ierr /=0) call errore('surfacehopping','Error allocating ge',1)
    allocate(gh(1:nefre),stat=ierr)  !g_ij
    if(ierr /=0) call errore('surfacehopping','Error allocating gh',1)    
    allocate(phQ0(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ0',1)
    allocate(phP0(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP0',1)    
    allocate(ph_U(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ph_U',1)
    allocate(ph_T(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ph_T',1)            
    allocate(e0(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating e0',1)
    allocate(p0(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating p0',1)
    allocate(d0(nefre,nefre,nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating d0',1)
    allocate(ge1(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ge1',1) 
    allocate(gh1(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating gh1',1) 
    allocate(pes(0:nefre,1:nsnap,1:naver),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating pes',1)
    pes = 0.0
    allocate(inf(1:3,1:nsnap,1:naver),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating inf',1)
    inf = 0.0
    allocate(csit(nefre,nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating csit',1)
    csit = 0.0
    allocate(wsit(nefre,nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating wsit',1)
    wsit = 0.0 
    allocate(psit(nefre,nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating psit',1)   
    psit = 0.0
    allocate(xsit(nefre,nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating xsit',1)  
    xsit = 0.0
    allocate(ksit(nefre,nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ksit',1)    
    ksit = 0.0
    allocate(msd(nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating msd',1)
    msd = 0.0
    allocate(ipr(nsnap),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ipr',1)
    ipr = 0.0
    allocate(msds(nsnap,naver),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating csit',1)
    msds = 0.0
    !allocate(c(nefre),stat=ierr)
    !if(ierr /=0) call errore('surfacehopping','Error allocating c',1)
    allocate(celec_nk(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating celec_nk',1)  
    allocate(chole_nk(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating chole_nk',1)      
    !allocate(w(nefre),stat=ierr)
    !if(ierr /=0) call errore('surfacehopping','Error allocating w',1)
    allocate(w_e(nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w_e',1)
    allocate(w_h(nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w_h',1)      
    allocate(w0_e(nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w0_e',1) 
    allocate(w0_h(nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w0_h',1)     
  end subroutine allocatesh
  
  

  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% convert wavefunction from diabatix to adiabatic basis %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  subroutine convert_diabatic_adiabatic(pp_nk,c_nk,ww)
    use f95_precision
    use blas95
    implicit none
    real(kind=dp),intent(in) :: pp_nk(nbndfst,nktotf,nefre)
    complex(kind=dpc),intent(in) :: c_nk(nbndfst,nktotf)
    complex(kind=dpc),intent(out):: ww(nefre)
    
    real(kind=dp) :: pp(nefre,nefre)
    complex(kind=dpc) :: cc(nefre) 
    
    integer :: iefre,jefre,ik,jk,iband,jband 
    
    pp = reshape(pp_nk,(/nefre,nefre/))
    cc = reshape(c_nk,(/nefre/))
    ww= 0.0d0
    call gemv(pp,cc,ww,trans='T')
    
    !ww=0.0d0
    !do ik=1,nktotf
    !  do iband=1,nbndfst
    !    iefre = (ik-1)*nbndfst + iband
    !    do jk=1,nktotf
    !      do jband=1,nbndfst
    !        ww(iefre) = ww(iefre)+pp_nk(jband,jk,iefre)*c_nk(jband,jk)
    !      enddo
    !    enddo
    !  enddo
    !enddo
    !
    !do iefre=1,nefre
    !  if(ww(iefre)/=ww_(iefre)) then
    !    lgemv=.false.
    !    exit
    !  endif
    !enddo
      
    
  
  end subroutine convert_diabatic_adiabatic
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% calculate nonadiabatic coupling %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  ! ref : PPT-91
  subroutine calculate_nonadiabatic_coupling(nmodes,ee,p_nk,dd)
    use kinds,only :  dp
    implicit none
    integer, intent(in) :: nmodes
    real(kind=dp),intent(in) :: ee(nefre)
    !real(kind=dp),intent(in) :: epcq(nbndfst,nbndfst,nktotf,nmodes,nqtotf)
    real(kind=dp),intent(in) :: p_nk(nbndfst,nktotf,nefre)
    real(kind=dp),intent(out):: dd(nefre,nefre,nmodes,nqtotf)
    integer :: iefre,jefre,iq,imode 
    integer :: ik,ikq,iband1,iband2
    
    dd=0.0d0
    do iefre=1,nefre-1
      do jefre=iefre+1,nefre
        do iq =1, nqtotf
          do imode=1,nmodes
            do ik=1,nktotf
              ikq = kqmap(ik,iq)
              do iband1=1,nbndfst
                do iband2=1,nbndfst
                  dd(iefre,jefre,imode,iq) = dd(iefre,jefre,imode,iq) + &
                  p_nk(iband1,ik,iefre)*p_nk(iband2,ikq,jefre)*epcq(iband1,iband2,ik,imode,iq)
                enddo
              enddo
            enddo
            dd(iefre,jefre,imode,iq) = sqrt(2.0*wf(imode,iq)/nqtotf)*dd(iefre,jefre,imode,iq)/(ee(jefre)-ee(iefre))
          enddo
        enddo
        !dd(iefre,jefre,:,:) = dd(iefre,jefre,:,:)/(ee(jefre)-ee(iefre))
        dd(jefre,iefre,:,:) = - dd(iefre,jefre,:,:)
      enddo
    enddo    
    
  end subroutine calculate_nonadiabatic_coupling
  
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% calculate eigenenergy and eigenstate %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  
  function bolziman(womiga,temp)
    use kinds ,only : dp
    implicit none
    real(kind=dp)::womiga,temp,bolziman

    bolziman=1.0/(exp(womiga/(temp))-1.0)
    !<nb>=1/(exp{hw/kbT}-1)
  end function bolziman
 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% CALCULATE HOPPING PROBABILITY %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% REF: NOTEBOOK PAGE 631        %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  !call calculate_hopping_probability(w0_e,phP0,d0,dt,g,g1) using FSSH in adiabatic representation
  !ref : 1 J. C. Tully, J. Chem. Phys. 93 (1990) 1061.
  !ref : 2 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
  !ref : eq(2)
  ! 其中gg1为原始跃迁几率，gg为采用FSSH方法，令跃迁几率小于0的部分等于0
  ! if(gg(iefre) < 0.0) gg(iefre) = 0.0
  subroutine calculate_hopping_probability(isurface,WW,VV,dd,tt,gg,gg1)
    use modes,only : nmodes
    implicit none
    integer,intent(in)           :: isurface
    complex(kind=dpc),intent(in) :: WW(nefre)
    real(kind=dp),intent(in)     :: VV(nmodes,nqtotf)
    real(kind=dp),intent(in)     :: dd(nefre,nefre,nmodes,nqtotf)
    real(kind=dp),intent(in)     :: tt
    real(kind=dp),intent(out)     :: gg(nefre)
    real(kind=dp),intent(out)     :: gg1(nefre)
    
    real(kind=dp) :: sumvd
    integer :: iefre,iq,imode
    
    gg = 0.0
    gg1= 0.0
    ! FSSH
    ! ref: 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
    do iefre=1,nefre
      if(iefre /= isurface) then
        sumvd = 0.0
        do iq=1,nqtotf
          do imode=1,nmodes
            sumvd = sumvd+VV(imode,iq)*dd(isurface,iefre,imode,iq)
          enddo
        enddo
        ! in adiabatic representation：the switching probabilities from the active surface isurface to another surface iefre 
        gg(iefre)=2.0*tt*Real(CONJG(WW(isurface))*WW(iefre))*sumvd/REAL(CONJG(WW(isurface))*WW(isurface))
        gg1(iefre) = gg(iefre)  ! 绝热表象原始的跃迁几率
        !FSSH if g_ij<0,reset to g_ij=0
        if(gg(iefre) < 0.0) gg(iefre) = 0.0
      endif
    enddo
      
  end subroutine calculate_hopping_probability

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% CALCULATE SUMG0,SUMG1,MINDE %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  subroutine calculate_sumg_pes(sumg0,sumg1,w0,w,gg1,gg,isurface,isurface_j,minde)
    use constants,only : cmplx_0
    implicit none
    real(kind=dp),intent(out):: sumg0,sumg1
    complex(kind=dpc),intent(in) :: w0(nefre),w(nefre)
    real(kind=dp),intent(in) :: gg1(nefre)
    real(kind=dp),intent(out):: gg(nefre)
    integer,intent(in)       :: isurface 
    integer,intent(out)      :: isurface_j
    real(kind=dp),intent(out):: minde 
    real(kind=dp),allocatable :: S_ai(:)
    real(kind=dp) :: S_aa,SUM_S
    integer :: max_Sai(1)
    integer :: iefre,isurface_a
    logical :: lallocate
    
    ! sumg0 总的跃迁几率(透热表象下计算得出)
    ! SC_FSSH
    if(methodsh == "SC-FSSH") then
      ! ref: 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
      ! eq(4) ,eq(5)
      ! ref: 1 L. Wang, and O. V. Prezhdo, Journal of Physical Chemistry Letters 5 (2014) 713.
      ! eq(4)-eq(12)
      ! Assumption at most one trivial crossing is encountered during a time step.
      sumg0 = (ABS(W0(ISURFACE))**2-ABS(W(ISURFACE))**2)/ABS(W0(ISURFACE))**2
      sumg1 = SUM(gg1)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !% CHANGE POTENTIAL ENERGY SURFACE %!
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!      
      IF(ISURFACE == 1) THEN
        MINDE=(E0(ISURFACE+1)-E0(ISURFACE))
        iefre=isurface+1
      ELSEIF(ISURFACE == nefre) THEN
        MINDE=(E0(ISURFACE)-E0(ISURFACE-1))
        iefre = isurface-1
      ELSEIF((E0(ISURFACE+1)-E0(ISURFACE)) < (E0(ISURFACE)-E0(ISURFACE-1))) THEN
        MINDE=(E0(ISURFACE+1)-E0(ISURFACE))
        iefre= isurface+1
      ELSE
        MINDE=(E0(ISURFACE)-E0(ISURFACE-1))
        iefre= isurface-1
      ENDIF 
      
      ! eq(13)
      gg(iefre) = sumg0 - (sumg1-gg1(iefre))
      if(gg(iefre) < 0.0) gg(iefre)=0.0
      if(SUM(GG) > 1.0 ) GG=GG/SUM(GG)
    elseif(methodsh == "CC-FSSH") then
      ! ref : 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
      lallocate = allocated(S_ai)
      if(.not. lallocate) allocate(S_ai(nefre))
      S_ai = 0.0
      isurface_a = isurface
      S_aa = SUM(p0(:,isurface_a)*p(:,isurface_a))
      
      !S_aa = SUM(CONJG(p0(:,isurface_a))*p(:,isurface_a))
      if(S_aa**2 > 0.5) then
        !type(1) or type(3)
        isurface_j = isurface_a
      else
        !type(2) or type(4)
        do iefre = 1,nefre
          S_ai(iefre) = SUM(p0(:,isurface_a)*p(:,iefre))
          !S_ai(ibasis) = SUM(CONJG(pp0(:,isurface_a))*pp(:,ibasis))
        enddo
        S_ai = S_ai**2
        SUM_S = SUM(S_ai)   ! = 1.00
        !S_ai = CONJG(S_ai)*S_ai
        max_Sai =  MAXLOC(S_ai)
        isurface_j = max_Sai(1)
      endif
      
      sumg0 = (ABS(W0(ISURFACE))**2-ABS(W(ISURFACE))**2)/ABS(W0(ISURFACE))**2
      sumg1 = SUM(gg1)      
      if(isurface_j == isurface) then
        ! type(1) or type(3)
        gg = gg
      else
        gg(isurface_j) = sumg0 - (SUM(GG1)-GG1(isurface_j))
        if(GG(isurface_j) < 0.0d0) GG(isurface_j) = 0.0d0
        if(SUM(GG) > 1.0d0) GG=GG/SUM(GG)
      endif
    !
    endif    

    
  end subroutine calculate_sumg_pes
  
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% NONADIABATIC TRANSITION %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% REF: NOTEBOOK PAGE 635  %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  subroutine nonadiabatic_transition(lelec,isurface,isurface_j,surface_type,EE,PP,DD,GG,WW,VV)
    use randoms,only :more_random
    use modes,only : nmodes
    implicit none
    logical , intent(in) :: lelec
    integer , intent(inout)  :: isurface
    integer , intent(in)     :: isurface_j
    integer , intent(out)    :: surface_type
    real(kind=dp),intent(in) :: EE(nefre)
    real(kind=dp),intent(in) :: PP(nefre,nefre)
    real(kind=dp),intent(in) :: DD(nefre,nefre,nmodes,nqtotf)
    real(kind=dp),intent(in) :: GG(nefre)
    real(kind=dp),intent(inout) :: VV(nmodes,nqtotf)
    complex(kind=dpc),intent(in) :: WW(nefre)
    
    real(kind=dp),allocatable :: S_bi(:)
    real(kind=dp) :: sumvd,sumdd,sumgg,flagr,flagd,detaE
    integer :: iefre,jefre,imode,iq,isurface_a,isurface_b,isurface_k
    logical :: lallocate
    real(kind=dp) :: SUM_S,max_Sbi(1),SUM_G
    
    isurface_a = isurface
    call more_random()
    call random_number(flagr)
    sumgg = 0.0d0
    SUM_G = SUM(GG)
    do iefre=1,nefre
      if(iefre /= isurface) then
        sumgg = sumgg + GG(iefre)
        if(flagr < sumgg) then
          isurface_b = iefre
          if(methodsh == "CC-FSSH") then
            !if(isurface_b /= isurface_j) then
            lallocate = allocated(S_bi)
            if(.not. lallocate) allocate(S_bi(nefre))
            S_bi = 0.0

            do jefre = 1,nefre
              S_bi(jefre) = SUM(p0(:,isurface_b)*p(:,jefre))
              !S_bi(ibasis) = SUM(CONJG(pp0(:,isurface_b))*pp(:,ibasis))
            enddo
            SUM_S = SUM(S_bi) 
            S_bi = S_bi**2
            SUM_S = SUM(S_bi)
            !S_ai = CONJG(S_ai)*S_ai
            max_Sbi =  MAXLOC(S_bi)
            isurface_k = max_Sbi(1)
            !endif
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
          endif
          
          ! ref:1 L. Wang, and D. Beljonne, Journal of Physical Chemistry Letters 4 (2013) 1888.
          ! eq(16)
          sumvd = 0.0
          sumdd = 0.0
          do iq=1,nqtotf
            do imode=1,nmodes
              sumvd = sumvd + VV(imode,iq)*DD(isurface_a,iefre,imode,iq) ! A
              sumdd = sumdd + DD(isurface_a,iefre,imode,iq)**2           ! B
            enddo
          enddo
          detaE = EE(isurface_a)-EE(iefre)
          if(.not. lelec) detaE = -1.0*detaE  ! 针对空穴修改能量
          flagd = 1.0+2.0*detaE*sumdd/sumvd**2  
          
          if(methodsh == "CC-FSSH") then
            if(isurface_j == isurface_a) then ! type (1) (3)
              if(flagd > 0.0) then
                flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
                do iq=1,nqtotf
                  do imode=1,nmodes
                    VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,iefre,imode,iq)
                  enddo
                enddo
                isurface = isurface_k
              else
                isurface = isurface_a
              endif
            elseif(isurface_j /= isurface_a) then !type (2)(4)
              if(isurface_b /= isurface_j) then!b/=j
                if(flagd > 0.0) then
                  flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
                  do iq=1,nqtotf
                    do imode=1,nmodes
                      VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,iefre,imode,iq)
                    enddo
                  enddo
                  isurface = isurface_k
                else
                  isurface = isurface_j
                endif
              else  !b=j
                isurface = isurface_j
              endif
            
            endif
            
            
          else
            if(flagd > 0.0) then
              flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
              do iq=1,nqtotf
                do imode=1,nmodes
                  VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,iefre,imode,iq)
                enddo
              enddo
              isurface = iefre          
            endif
          endif
          
          exit
        endif
      endif
    enddo
  
  end subroutine nonadiabatic_transition
  
end module surfacehopping