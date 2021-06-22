module surfacehopping
  use kinds, only : dp,dpc
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,nktotf,nbndfst,ibndmin,ibndmax
  use hamiltonian,only : nphfre,neband,nhband,nefre,nhfre,&
                         E_e,P_e,P_e_nk,E0_e,P0_e,P0_e_nk,&
                         E_h,E_h_,P_h,P_h_nk,E0_h,P0_h,P0_h_nk
  use parameters, only : nsnap,naver,ncore,nnode
  use surfacecom, only : iesurface,ihsurface,esurface_type,hsurface_type,&
                         phQ,phP,phQ0,phP0,phK,phU,SUM_phU,SUM_phK,SUM_phE,&
                         phQsit,phPsit,phKsit,phUsit,ld_gamma,&
                         dEa_dQ,dEa_dQ_e,dEa_dQ_h,dEa2_dQ2,dEa2_dQ2_e,dEa2_dQ2_h,&
                         d_e,g_e,g1_e,c_e_nk,w_e,w0_e,&
                         d0_e,&
                         d_h,g_h,g1_h,c_h_nk,w_h,w0_h,&
                         d0_h,&
                         pes_e,csit_e,wsit_e,psit_e,mskds_e,mskd_e,ipr_e,&
                         pes_h,csit_h,wsit_h,psit_h,mskds_h,mskd_h,ipr_h,&
												 apes_e,apes_sum_e,apes_h,apes_sum_h
                         
  use cc_fssh,only : S_ai_e,S_ai_h,S_bi_e,S_bi_h
                  
  implicit none
  
  complex(kind=dpc),allocatable :: cc0_e(:),dc1_e(:),dc2_e(:),dc3_e(:),dc4_e(:)
  complex(kind=dpc),allocatable :: cc0_h(:),dc1_h(:),dc2_h(:),dc3_h(:),dc4_h(:)
  real(kind=dp)    ,allocatable :: n_e(:),n_h(:)
  
  contains
  
  subroutine allocatesh(methodsh,lelecsh,lholesh,nmodes,nq)
    implicit none
    character(len=*),intent(in) :: methodsh
    logical,intent(in):: lelecsh,lholesh
    integer,intent(in):: nmodes,nq
    integer :: ierr
    
    
    allocate(ld_gamma(nmodes,nq))
    ld_gamma = 0.0
    allocate(phQ(nmodes,nq),stat=ierr)  ! x(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ',1)
    phQ = 0.0
    ! phonons normal mode coordinates
    allocate(phQsit(nmodes,nq,0:nsnap))
    allocate(phPsit(nmodes,nq,0:nsnap))
    allocate(phKsit(nmodes,nq,0:nsnap))
    allocate(phUsit(nmodes,nq,0:nsnap))
    phQsit = 0.0d0
    phPsit = 0.0d0
    phKsit = 0.0d0
    phUsit = 0.0d0
    
    allocate(phP(nmodes,nq),stat=ierr)  ! v(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP',1)    
    phP = 0.0
    !phonons normal mode verlosity
    allocate(phQ0(nmodes,nq),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ0',1)
    allocate(phP0(nmodes,nq),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP0',1)    
    allocate(phU(nmodes,nq),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phU',1)
    allocate(phK(nmodes,nq),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phK',1)            
    allocate(dEa_dQ(nmodes,nq))
    allocate(dEa2_dQ2(nmodes,nq))
    
    if(lelecsh) then
      allocate(E_e(1:nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_e',1)
      allocate(P_e(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_e',1)
      allocate(P_e_nk(neband,nktotf,nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_e_nk',1)
      allocate(d_e(nefre,nefre,nmodes,nq),stat=ierr) !d_ijk
      if(ierr /=0) call errore('surfacehopping','Error allocating d_e',1)
      allocate(dEa_dQ_e(nmodes,nq))
      allocate(dEa2_dQ2_e(nmodes,nq))
      
      allocate(c_e_nk(neband,nktotf),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_e_nk',1)
      
      allocate(pes_e(0:nefre,0:nsnap,1:naver))
			allocate(apes_e(0:nsnap,1:naver),apes_sum_e(0:nsnap,1:naver*ncore*nnode))
      allocate(csit_e(nefre,0:nsnap))
      allocate(wsit_e(nefre,0:nsnap))
      allocate(psit_e(nefre,0:nsnap))
			allocate(mskds_e(0:nsnap,naver),mskd_e(0:nsnap),ipr_e(0:nsnap))
      pes_e  = 0.0d0
      csit_e = 0.0d0
      wsit_e = 0.0d0
      psit_e = 0.0d0
      mskds_e= 0.0
			mskd_e = 0.0
			ipr_e  = 0.0
			
			
      if(methodsh == "SC-FSSH" .OR. methodsh == "CC-FSSH") then
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
      endif
      
      allocate(w_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w_e',1)  
      allocate(w0_e(nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w0_e',1)
      
      allocate(g_e(1:nefre),stat=ierr)  !g_ij
      if(ierr /=0) call errore('surfacehopping','Error allocating g_e',1)
      allocate(g1_e(1:nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating g1_e',1) 
      
      allocate(E0_e(1:nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E0_e',1)
      allocate(P0_e(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_e',1)
      allocate(d0_e(nefre,nefre,nmodes,nq),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating d0_e',1)
      
      if(methodsh == "CC-FSSH") then
        allocate(S_ai_e(nefre))
        allocate(S_bi_e(nefre))
      endif
      
    endif
    
    if(lholesh) then
      allocate(E_h(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_h',1)
      allocate(P_h(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_h',1)
      allocate(P_h_nk(nhband,nktotf,nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P_h_nk',1)
      allocate(d_h(nhfre,nhfre,nmodes,nq),stat=ierr) !d_ijk
      if(ierr /=0) call errore('surfacehopping','Error allocating d_h',1)
      allocate(dEa_dQ_h(nmodes,nq))
      allocate(dEa2_dQ2_h(nmodes,nq))

      allocate(c_h_nk(nhband,nktotf),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating c_h_nk',1)  
      
      allocate(pes_h(0:nefre,0:nsnap,1:naver))
			allocate(apes_h(0:nsnap,1:naver),apes_sum_h(0:nsnap,1:naver*ncore*nnode))
      allocate(csit_h(nhfre,0:nsnap))
      allocate(wsit_h(nhfre,0:nsnap))
      allocate(psit_h(nhfre,0:nsnap))
		  allocate(mskds_h(0:nsnap,naver),mskd_h(0:nsnap),ipr_h(0:nsnap))
      pes_h  = 0.0
      csit_h = 0.0
      wsit_h = 0.0
      psit_h = 0.0
      mskds_h= 0.0
			mskd_h = 0.0
			ipr_h  = 0.0
      if(methodsh == "SC-FSSH" .OR. methodsh == "CC-FSSH") then
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
      endif
      
      allocate(w_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w_h',1)  
      allocate(w0_h(nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating w0_h',1)
      
      allocate(g_h(1:nhfre),stat=ierr)  !g_ij
      if(ierr /=0) call errore('surfacehopping','Error allocating g_h',1)
      allocate(g1_h(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating g1_h',1) 
      
      allocate(E0_h(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E0_h',1)
      allocate(E_h_(1:nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating E_h_',1)      
      
      allocate(P0_h(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating P0_h',1)
      allocate(d0_h(nhfre,nhfre,nmodes,nq),stat=ierr)
      if(ierr /=0) call errore('surfacehopping','Error allocating d0_h',1)
      
      if(methodsh == "CC-FSSH") then
        allocate(S_ai_h(nhfre))
        allocate(S_bi_h(nhfre))
      endif
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
    real(kind=dp) :: sum_cc2,sum_ww2
    
    sum_cc2 = REAL(SUM(cc*CONJG(cc)))
    
    ww= 0.0d0
    !call gemv(pp,cc,ww)
    call gemv(pp,cc,ww,trans='T')
    
    sum_ww2 = REAL(SUM(ww*CONJG(ww)))
    
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
                  dd(ifre,jfre,imode,iq) = dd(ifre,jfre,imode,iq)+&
                  &p_nk(iband1,ik,ifre)*p_nk(iband2,ikq,jfre)*epcq(iband1,iband2,ik,imode,iq)
                enddo
              enddo
            enddo
            dd(ifre,jfre,imode,iq) = sqrt(2.0*wf(imode,iq)/nq)*dd(ifre,jfre,imode,iq)/(ee(jfre)-ee(ifre))
          enddo
        enddo
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



end module surfacehopping