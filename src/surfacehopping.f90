module surfacehopping
  use kinds, only : dp,dpc
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,epcq,nktotf,nbndfst,nqtotf,ibndmin,ibndmax
  use hamiltonian,only : nphfre,nefre
  use parameters, only : nsnap,naver
  implicit none
  integer :: iaver
  integer :: ierr
  real(kind=dp),allocatable :: phQ(:),phP(:)
  ! phonons normal mode coordinate,and phonons P
  real(kind=dp),allocatable :: e(:),p(:,:),d(:,:,:),g(:)
  real(kind=dp),allocatable :: phQ0(:),phP0(:)
  real(kind=dp),allocatable :: e0(:),p0(:,:),d0(:,:,:),g1(:)  
  real(kind=dp),allocatable :: pes(:,:,:),inf(:,:,:),csit(:,:),wsit(:,:),&
                               psit(:,:),xsit(:,:),ksit(:,:) 
  real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  real(kind=dp),allocatable :: c(:),w(:),w0(:)
  real(kind=dp),allocatable :: c_nk(:,:),w_nk(:,:),w0_nk(:,:)
  contains
  
  subroutine allocatesh()
    implicit none
    allocate(phQ(1:nphfre),stat=ierr)  ! x(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ',1)
    phQ = 0.0
    allocate(phP(1:nphfre),stat=ierr)  ! v(1:nphfre)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP',1)    
    phP=0.0
    allocate(e(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating e',1)
    allocate(p(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating p',1)
    !allocate(d(nefre,nefre,nphfre),stat=ierr) !d_ijk
    !if(ierr /=0) call errore('surfacehopping','Error allocating d',1)
    allocate(g(1:nefre),stat=ierr)  !g_ij
    if(ierr /=0) call errore('surfacehopping','Error allocating g',1)
    allocate(phQ0(1:nphfre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phQ0',1)
    allocate(phP0(1:nphfre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating phP0',1)    
    allocate(e0(1:nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating e0',1)
    allocate(p0(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating p0',1)
    !allocate(d0(nefre,nefre,nphfre),stat=ierr)
    !if(ierr /=0) call errore('surfacehopping','Error allocating d0',1)
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
    allocate(c(nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating c',1)
    allocate(c_nk(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating c_nk',1)    
    allocate(w(nefre),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w',1)
    allocate(w_nk(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w_nk',1)
    allocate(w0_nk(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore('surfacehopping','Error allocating w0_nk',1) 

    
  end subroutine allocatesh
  
  !=============================================!
  != init coordinate and Normal mode velocitie =!
  !=============================================!
  subroutine init_normalmode_coordinate_velocity(nq,nmodes,ph_l,ph_p,wf,T)
    use kinds,only   : dp
    use randoms,only : gaussian_random_number
    use phdisp, only : ph_lqv,ph_pqv
    implicit none
    integer      ,intent(in) :: nq,nmodes
    real(kind=dp),intent(in) :: T
    real(kind=dp),intent(in) :: wf(nmodes,nq)
    real(kind=dp),intent(out):: ph_l(nmodes,nq),ph_P(nmodes,nq)
    
    integer :: iq,imode
    
    ! ph_l= Q_qv/ sqrt(hbar/(2*w_qv))
    ! ph_P= P_qv/ sqrt(hbar/(2*w_qv))
    
    !ph_l = gaussian_random_number(0.0d0,ph_lqv)
    !ph_P = gaussian_random_number(0.0d0,ph_pqv)
    
    do iq=1,nq
      do imode=1,nmodes
        ph_l(imode,iq) = gaussian_random_number(0.0d0,ph_lqv(imode,iq))
        ph_P(imode,iq) = gaussian_random_number(0.0d0,ph_pqv(imode,iq))
      enddo
    enddo
    
  end subroutine init_normalmode_coordinate_velocity
  

  function bolziman(womiga,temp)
    use kinds ,only : dp
    implicit none
    real(kind=dp)::womiga,temp,bolziman

    bolziman=1.0/(exp(womiga/(temp))-1.0)
    !<nb>=1/(exp{hw/kbT}-1)
  end function bolziman
 
 

end module surfacehopping