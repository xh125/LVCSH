module surfacehopping
  use kinds, only : dp,dpc
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,nktotf,nbndfst,nqtotf,ibndmin,ibndmax
  use hamiltonian,only : nphfre,neband,nhband,nefre,nhfre,&
                         E_e,P_e,P_e_nk,E0_e,P0_e,P0_e_nk,&
                         E_h,P_h,P_h_nk,E0_h,P0_h,P0_h_nk
  use parameters, only : nsnap,naver
  use surfacecom, only : MethodSH,iesurface,ihsurface,esurface_type,hsurface_type,&
                         phQ,phP,phQ0,phP0,ph_T,ph_U,SUM_ph_U,SUM_ph_T,SUM_ph_E,&
                         d_e,ge,ge1,c_e_nk,w_e,w0_e,&
                         d0_e,&
                         d_h,gh,gh1,c_h_nk,w_h,w0_h,&
                         d0_h,&
                         minde_e,minde_h,sumg0_e,sumg0_h,sumg1_e,sumg1_h
  implicit none
  
  complex(kind=dpc),allocatable :: cc0_e(:),dc1_e(:),dc2_e(:),dc3_e(:),dc4_e(:)
  complex(kind=dpc),allocatable :: cc0_h(:),dc1_h(:),dc2_h(:),dc3_h(:),dc4_h(:)
  real(kind=dp),allocatable :: n_e(:),n_h(:)
  
  contains
  
  subroutine allocatesh(lelecsh,lholesh,nmodes)
    implicit none
    logical,intent(in):: lelecsh,lholesh
    integer,intent(in):: nmodes
    integer :: ierr
    
    allocate(phQ(nmodes,nqtotf),stat=ierr)  ! x(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ',1)
    phQ = 0.0
    ! phonons normal mode coordinates
    allocate(phP(nmodes,nqtotf),stat=ierr)  ! v(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP',1)    
    phP = 0.0
    !phonons normal mode verlosity
    allocate(phQ0(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ0',1)
    allocate(phP0(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP0',1)    
    allocate(ph_U(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ph_U',1)
    allocate(ph_T(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating ph_T',1)            

    
    if(lelecsh) then
      allocate(E_e(1:nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_e',1)
      allocate(P_e(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_e',1)
      allocate(P_e_nk(neband,nktotf,nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_e_nk',1)
      allocate(d_e(nefre,nefre,nmodes,nqtotf),stat=ierr) !d_ijk
      if(ierr /=0) call errore('surfacehopping','Error allocating d_e',1)

      allocate(c_e_nk(neband,nktotf),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_e_nk',1)
      allocate(cc0_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating cc0_e',1)
      allocate(dc1_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc1_e',1)
      allocate(dc2_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc2_e',1)
      allocate(dc3_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc3_e',1)
      allocate(dc4_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc4_e',1)
      allocate(n_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating n_e',1)            
      
      allocate(w_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w_e',1)  
      allocate(w0_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w0_e',1)
      
      allocate(ge(1:nefre),stat=ierr)  !g_ij
      if(ierr /=0) call errore('surfacehopping','Error allocating ge',1)
      allocate(ge1(1:nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating ge1',1) 
      
      allocate(E0_e(1:nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E0_e',1)
      allocate(P0_e(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_e',1)
      allocate(d0_e(nefre,nefre,nmodes,nqtotf),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating d0_e',1)

    endif
    
    if(lholesh) then
      allocate(E_h(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_h',1)
      allocate(P_h(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_h',1)
      allocate(P_h_nk(nhband,nktotf,nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_h_nk',1)
      allocate(d_h(nhfre,nhfre,nmodes,nqtotf),stat=ierr) !d_ijk
      if(ierr /=0) call errore('surfacehopping','Error allocating d_h',1)

      allocate(c_h_nk(nhband,nktotf),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_h_nk',1)  
      allocate(cc0_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating cc0_h',1)
      allocate(dc1_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc1_h',1)
      allocate(dc2_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc2_h',1)
      allocate(dc3_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc3_h',1)
      allocate(dc4_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating dc4_h',1)
      allocate(n_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating n_h',1)                  
      allocate(w_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w_h',1)  
      allocate(w0_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w0_h',1)
      
      allocate(gh(1:nhfre),stat=ierr)  !g_ij
      if(ierr /=0) call errore('surfacehopping','Error allocating gh',1)
      allocate(gh1(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating gh1',1) 
      
      allocate(E0_h(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E0_h',1)
      allocate(P0_h(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_h',1)
      allocate(d0_h(nhfre,nhfre,nmodes,nqtotf),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating d0_h',1)
    endif
    
      !allocate(pes(0:nefre,1:nsnap,1:naver),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating pes',1)
      !pes = 0.0
      !allocate(inf(1:3,1:nsnap,1:naver),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating inf',1)
      !inf = 0.0
      !allocate(csit(nefre,nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating csit',1)
      !csit = 0.0
      !allocate(wsit(nefre,nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating wsit',1)
      !wsit = 0.0 
      !allocate(psit(nefre,nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating psit',1)   
      !psit = 0.0
      !allocate(xsit(nefre,nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating xsit',1)  
      !xsit = 0.0
      !allocate(ksit(nefre,nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating ksit',1)    
      !ksit = 0.0
      !allocate(msd(nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating msd',1)
      !msd = 0.0
      !allocate(ipr(nsnap),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating ipr',1)
      !ipr = 0.0
      !allocate(msds(nsnap,naver),stat=ierr)
      !if(ierr /=0) call errore('surfacehopping','Error allocating csit',1)
      !msds = 0.0
    
    
  end subroutine allocatesh
  
  

  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% convert wavefunction from diabatix to adiabatic basis %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  subroutine convert_diabatic_adiabatic(nfre,pp,cc,ww)                                        
    use f95_precision
    use blas95
    implicit none
    integer,intent(in) :: nfre 
    real(kind=dp),intent(in) :: pp(nfre,nfre)
    complex(kind=dpc),intent(in) :: cc(nfre)
    complex(kind=dpc),intent(out):: ww(nfre)
    
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
  subroutine calculate_nonadiabatic_coupling(nmodes,nq,nband,nk,wf,ee,p_nk,epcq,dd)
    use kinds,only :  dp
    implicit none
    integer, intent(in) :: nmodes,nq,nband,nk
    real(kind=dp),intent(in) :: wf(nmodes,nq),ee(nband*nk)
    real(kind=dp),intent(in) :: epcq(nband,nband,nk,nmodes,nq)
    real(kind=dp),intent(in) :: p_nk(nband,nk,nband*nk)
    real(kind=dp),intent(out):: dd(nband*nk,nband*nk,nmodes,nq)
    integer :: nfre,ifre,jfre,iq,imode 
    integer :: ik,ikq,iband1,iband2
    
    nfre = nband*nk
    
    dd=0.0d0
    do ifre=1,nfre-1
      do jfre=ifre+1,nfre
        do iq =1, nq
          do imode=1,nmodes
            do ik=1,nk
              ikq = kqmap(ik,iq)
              do iband1=1,nband
                do iband2=1,nband
                  dd(ifre,jfre,imode,iq) = dd(ifre,jfre,imode,iq) + &
                  p_nk(iband1,ik,ifre)*p_nk(iband2,ikq,jfre)*epcq(iband1,iband2,ik,imode,iq)
                enddo
              enddo
            enddo
            dd(ifre,jfre,imode,iq) = sqrt(2.0*wf(imode,iq)/nq)*dd(ifre,jfre,imode,iq)/(ee(jfre)-ee(ifre))
          enddo
        enddo
        !dd(iefre,jefre,:,:) = dd(iefre,jefre,:,:)/(ee(jefre)-ee(iefre))
        dd(jfre,ifre,:,:) = - dd(ifre,jfre,:,:)
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
  subroutine calculate_hopping_probability(isurface,nfre,nmodes,nq,WW,VV,dd,tt,gg,gg1)
    implicit none
    integer,intent(in)           :: isurface,nfre,nmodes,nq
    complex(kind=dpc),intent(in) :: WW(nfre)
    real(kind=dp),intent(in)     :: VV(nmodes,nq)
    real(kind=dp),intent(in)     :: dd(nfre,nfre,nmodes,nq)
    real(kind=dp),intent(in)     :: tt
    real(kind=dp),intent(out)    :: gg(nfre)
    real(kind=dp),intent(out)    :: gg1(nfre)
    
    real(kind=dp) :: sumvd
    integer :: ifre,iq,imode
    
    gg = 0.0d0
    gg1= 0.0d0
    ! FSSH
    ! ref: 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
    do ifre=1,nfre
      if(ifre /= isurface) then
        sumvd = 0.0
        do iq=1,nq
          do imode=1,nmodes
            sumvd = sumvd+VV(imode,iq)*dd(isurface,ifre,imode,iq)
          enddo
        enddo
        ! in adiabatic representation：the switching probabilities from the active surface isurface to another surface iefre 
        gg(ifre)=2.0*tt*Real(CONJG(WW(isurface))*WW(ifre))*sumvd/REAL(CONJG(WW(isurface))*WW(isurface))
        gg1(ifre) = gg(ifre)  ! 绝热表象原始的跃迁几率
        !FSSH if g_ij<0,reset to g_ij=0
        if(gg(ifre) < 0.0d0) gg(ifre) = 0.0d0
      endif
    enddo
      
  end subroutine calculate_hopping_probability

  subroutine get_G_SC_FSSH(isurface,nfre,e0,w0,w,g1,g)
    use constants ,only   : ryd2eV
    implicit none
    integer,intent(in) :: isurface
    integer,intent(in) :: nfre
    real(kind=dp),intent(in)     :: e0(nfre)
    complex(kind=dpc),intent(in) :: w0(nfre),w(nfre)
    real(kind=dp),intent(in)     :: g1(nfre)
    real(kind=dp),intent(inout)  :: g(nfre)
    
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
    elseif(isurface == nfre) then
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

  subroutine get_G_CC_FSSH(nfre,isurface,p0,p,w0,w,S_ai,g1,g)
    implicit none
    integer,intent(in) :: nfre,isurface
    real(kind=dp),intent(in) :: p0(nfre,nfre)
    real(kind=dp),intent(in) :: p(nfre,nfre)
    complex(kind=dpc),intent(in):: w0(nfre),w(nfre)
    real(kind=dp),intent(inout) :: S_ai(nfre)
    real(kind=dp),intent(in)  :: g1(nfre)
    real(kind=dp),intent(inout) :: g(nfre)
    
    integer :: ifre
    integer :: isurface_a,isurface_b,isurface_j,isurface_k
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
      do ifre = 1,nfre
        S_ai(ifre) = SUM(p0(:,isurface_a)*p(:,ifre))
        !S_ai(ibasis) = SUM(CONJG(pp0(:,isurface_a))*pp(:,ibasis))
      enddo
      S_ai = S_ai**2
      SUM_S = SUM(S_ai)   ! = 1.00
      !S_ai = CONJG(S_ai)*S_ai
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

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% CALCULATE SUMG0,SUMG1,MINDE %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  subroutine calculate_sumg_pes(sumg0,sumg1,nfre,E0,P,P0,w0,w,gg1,gg,isurface,isurface_j,minde)
    use constants,only : cmplx_0
    implicit none
    real(kind=dp),intent(out):: sumg0,sumg1
    integer , intent(in)     :: nfre
    complex(kind=dpc),intent(in) :: w0(nfre),w(nfre)
    real(kind=dp),intent(in) :: gg1(nfre)
    real(kind=dp),intent(out):: gg(nfre)
    integer,intent(in)       :: isurface 
    integer,intent(out)      :: isurface_j
    real(kind=dp),intent(out):: minde 
    real(kind=dp),intent(in) :: E0(nfre)
    real(kind=dp),intent(in) :: P(nfre,nfre),P0(nfre,nfre)
    real(kind=dp),allocatable :: S_ai(:)
    real(kind=dp) :: S_aa,SUM_S
    integer :: max_Sai(1)
    integer :: ifre,isurface_a
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
        ifre=isurface+1
      ELSEIF(ISURFACE == nfre) THEN
        MINDE=(E0(ISURFACE)-E0(ISURFACE-1))
        ifre = isurface-1
      ELSEIF((E0(ISURFACE+1)-E0(ISURFACE)) < (E0(ISURFACE)-E0(ISURFACE-1))) THEN
        MINDE=(E0(ISURFACE+1)-E0(ISURFACE))
        ifre= isurface+1
      ELSE
        MINDE=(E0(ISURFACE)-E0(ISURFACE-1))
        ifre= isurface-1
      ENDIF 
      
      ! eq(13)
      gg(ifre) = sumg0 - (sumg1-gg1(ifre))
      if(gg(ifre) < 0.0) gg(ifre)=0.0
      if(SUM(GG) > 1.0 ) GG=GG/SUM(GG)
    elseif(methodsh == "CC-FSSH") then
      ! ref : 1 J. Qiu, X. Bai, and L. Wang, The Journal of Physical Chemistry Letters 9 (2018) 4319.
      lallocate = allocated(S_ai)
      if(.not. lallocate) allocate(S_ai(nfre))
      S_ai = 0.0
      isurface_a = isurface
      S_aa = SUM(p0(:,isurface_a)*p(:,isurface_a))
      
      !S_aa = SUM(CONJG(p0(:,isurface_a))*p(:,isurface_a))
      if(S_aa**2 > 0.5) then
        !type(1) or type(3)
        isurface_j = isurface_a
      else
        !type(2) or type(4)
        do ifre = 1,nfre
          S_ai(ifre) = SUM(p0(:,isurface_a)*p(:,ifre))
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
  subroutine nonadiabatic_transition(lelec,nfre,isurface,isurface_j,surface_type,EE,P,P0,DD,GG,WW,VV)
    use randoms,only :more_random
    use modes,only : nmodes
    implicit none
    logical , intent(in) :: lelec
    integer , intent(in)     :: nfre
    integer , intent(inout)  :: isurface
    integer , intent(in)     :: isurface_j
    integer , intent(out)    :: surface_type
    real(kind=dp),intent(in) :: EE(nfre)
    real(kind=dp),intent(in) :: P(nfre,nfre),P0(nfre,nfre)
    real(kind=dp),intent(in) :: DD(nfre,nfre,nmodes,nqtotf)
    real(kind=dp),intent(in) :: GG(nfre)
    real(kind=dp),intent(inout) :: VV(nmodes,nqtotf)
    complex(kind=dpc),intent(in) :: WW(nfre)
    
    real(kind=dp),allocatable :: S_bi(:)
    real(kind=dp) :: sumvd,sumdd,sumgg,flagr,flagd,detaE
    integer :: ifre,jfre,imode,iq,isurface_a,isurface_b,isurface_k
    logical :: lallocate
    real(kind=dp) :: SUM_S,max_Sbi(1),SUM_G
    
    isurface_a = isurface
    call more_random()
    call random_number(flagr)
    sumgg = 0.0d0
    SUM_G = SUM(GG)
    do ifre=1,nfre
      if(ifre /= isurface) then
        sumgg = sumgg + GG(ifre)
        if(flagr < sumgg) then
          isurface_b = ifre
          if(methodsh == "CC-FSSH") then
            !if(isurface_b /= isurface_j) then
            lallocate = allocated(S_bi)
            if(.not. lallocate) allocate(S_bi(nfre))
            S_bi = 0.0

            do jfre = 1,nfre
              S_bi(jfre) = SUM(p0(:,isurface_b)*p(:,jfre))
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
              sumvd = sumvd + VV(imode,iq)*DD(isurface_a,ifre,imode,iq) ! A
              sumdd = sumdd + DD(isurface_a,ifre,imode,iq)**2           ! B
            enddo
          enddo
          detaE = EE(isurface_a)-EE(ifre)
          if(.not. lelec) detaE = -1.0*detaE  ! 针对空穴修改能量
          flagd = 1.0+2.0*detaE*sumdd/sumvd**2  
          
          if(methodsh == "CC-FSSH") then
            if(isurface_j == isurface_a) then ! type (1) (3)
              if(flagd > 0.0) then
                flagd = sumvd/sumdd*(-1.0+dsqrt(flagd))
                do iq=1,nqtotf
                  do imode=1,nmodes
                    VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,ifre,imode,iq)
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
                      VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,ifre,imode,iq)
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
                  VV(imode,iq) = VV(imode,iq) + flagd*dd(isurface,ifre,imode,iq)
                enddo
              enddo
              isurface = ifre          
            endif
          endif
          
          exit
        endif
      endif
    enddo
  
  end subroutine nonadiabatic_transition
  
end module surfacehopping