module surfacehopping
  use kinds, only : dp,dpc
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,epcq,nktotf,nbndfst,nqtotf,ibndmin,ibndmax
  use hamiltonian,only : nphfre,nefre
  use parameters, only : nsnap,naver
  implicit none
  integer :: iaver
  integer :: isnap,istep
  integer :: iesurface,ihsurface
  integer :: ierr

  ! phonons normal mode coordinate,and phonons P
  real(kind=dp),allocatable :: phQ(:,:),phP(:,:),phQ0(:,:),phP0(:,:)
  real(kind=dp),allocatable :: e(:),p(:,:),p_nk(:,:,:),d(:,:,:,:),g(:)
  real(kind=dp),allocatable :: e0(:),p0(:,:),d0(:,:,:,:),g1(:)  
  real(kind=dp),allocatable :: pes(:,:,:),inf(:,:,:),csit(:,:),wsit(:,:),&
                               psit(:,:),xsit(:,:),ksit(:,:) 
  real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  complex(kind=dpc),allocatable :: celec_nk(:,:),w_e(:),w0_e(:)
  complex(kind=dpc),allocatable :: chole_nk(:,:),w_h(:),w0_h(:)
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
    allocate(g(1:nefre),stat=ierr)  !g_ij
    if(ierr /=0) call errore('surfacehopping','Error allocating g',1)
    allocate(phQ0(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ0',1)
    allocate(phP0(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP0',1)    
    allocate(e0(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating e0',1)
    allocate(p0(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating p0',1)
    allocate(d0(nefre,nefre,nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating d0',1)
    allocate(g1(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating g1',1) 
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
      !dE_dqv(iefre,:,:) = dd(iefre,iefre,:,:)
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
  !call calculate_hopping_probability(w0_e,phP0,d0,dt,g,g1)
  subroutine calculate_hopping_probability(WW,VV,dd,tt,gg,gg1)
    implicit none
    complex(kind=dpc),intent(in) :: WW(nefre)
    real(kind=dp),intent(in)     :: VV(nmodes,nqtotf)
    real(kind=dp),intent(in)     :: dd(nefre,nefre,nmodes,nqtotf)
    real(kind=dp),intent(in)     :: tt
    real(kind=dp),intent(out)     :: gg(nefre)
    real(kind=dp),intent(out)     :: gg1(nefre)
    
    real(kind=dp) :: sumvd
    
    gg = 0.0
    gg1= 0.0
    do iefre=1,nefre
      if(iefre /= iesurface) then
        sumvd = 0.0
        do iq=1,nqtotf
          do imode=1,nmodes
            sumvd = sumvd+VV(imode,iq)*dd(iesurface,iefre,imode,iq)
          enddo
        enddo
        gg(iefre)=2.0*tt*Real(CONJG(WW(iesurface))*WW(iefre))*sumvd/REAL(CONJG(WW(iesurface))*WW(iesurface))
        gg1(iefre) = gg(iefre)
        if(gg(iefre) < 0.0) gg(iefre) = 0.0
      endif
    enddo
      
  end subroutine calculate_hopping_probability
  
  
end module surfacehopping