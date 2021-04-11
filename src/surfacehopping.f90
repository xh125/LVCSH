module surfacehopping
  use kinds, only : dp,dpc
  use epwcom,only : nkf1,nkf2,nkf3,nqf1,nqf2,nqf3,kqmap
  use elph2,only  : wf,epcq,nktotf,nbndfst,nqtotf,ibndmin,ibndmax
  use hamiltonian,only : nphfre,nefre
  use parameters, only : nsnap,naver
  implicit none
  integer :: iaver
  integer :: isnap,istep
  integer :: isurface
  integer :: ierr

  ! phonons normal mode coordinate,and phonons P
  real(kind=dp),allocatable :: phQ(:,:),phP(:,:),phQ0(:,:),phP0(:,:)
  real(kind=dp),allocatable :: e(:),p(:,:),p_nk(:,:,:),d(:,:,:,:),g(:)
  real(kind=dp),allocatable :: e0(:),p0(:,:),d0(:,:,:,:),g1(:)  
  real(kind=dp),allocatable :: pes(:,:,:),inf(:,:,:),csit(:,:),wsit(:,:),&
                               psit(:,:),xsit(:,:),ksit(:,:) 
  real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  complex(kind=dpc),allocatable :: c(:),w(:),w0(:)
  complex(kind=dpc),allocatable :: c_nk(:,:),w_nk(:,:),w0_nk(:,:)
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
  subroutine init_normalmode_coordinate_velocity(nq,nmodes,ph_l,ph_p,ph_w,T)
    use kinds,only   : dp
    use randoms,only : gaussian_random_number
    use disp, only : ph_lqv,ph_pqv
    use constants,only : hbar_SI
    implicit none
    integer      ,intent(in) :: nq,nmodes
    real(kind=dp),intent(in) :: T
    real(kind=dp),intent(in) :: ph_w(nmodes,nq)    !rad/s
    real(kind=dp),intent(out):: ph_l(nmodes,nq),ph_p(nmodes,nq)
    
    integer :: iq,imode
    
    ! ph_l= Q_qv/ sqrt(hbar/(2*w_qv))
    ! ph_P= P_qv/ sqrt(hbar/(2*w_qv))
    
    !ph_l = gaussian_random_number(0.0d0,ph_lqv)
    !ph_P = gaussian_random_number(0.0d0,ph_pqv)
    
    do iq=1,nq
      do imode=1,nmodes
        ph_l(imode,iq) = gaussian_random_number(0.0d0,ph_lqv(imode,iq))
        ph_p(imode,iq) = gaussian_random_number(0.0d0,ph_pqv(imode,iq))
        !phQ(imode,iq) = ph_l(imode,iq) * sqrt(hbar_SI/2.0*ph_w(imode,iq))
      enddo
    enddo
    
    phQ = ph_l * sqrt(hbar_SI/2.0*ph_w)
    phP = ph_p * sqrt(hbar_SI/2.0*ph_w)
    
  end subroutine init_normalmode_coordinate_velocity
  
  !=============================================!
  != init dynamical varibale                   =!
  !=============================================!  
  subroutine init_dynamical_variable(nq,nmodes,ph_l,c_nk,ee,pp,ww)
    use parameters, only : init_cband,init_vband,init_ik
    use hamiltonian,only : H0_nk,H_nk,set_H_nk,calculate_eigen_energy_state
    implicit none
    integer,intent(in) :: nq,nmodes
    real(kind=dp),intent(in) :: ph_l(nmodes,nq)
    real(kind=dp),intent(out) :: ee(nbndfst*nktotf),pp(nbndfst*nktotf,nbndfst*nktotf)
    complex(kind=dpc),intent(out) :: c_nk(nbndfst,nktotf),ww(nbndfst*nktotf)
    integer :: iefre
    real(kind=dp) :: flagr,flagd
    
    c_nk = 0.0d0
    c_nk(init_cband,init_ik) = 1.0d0
    call set_H_nk(nq,nmodes,ph_l,nbndfst,nktotf,epcq,kqmap,H0_nk,H_nk)
    call calculate_eigen_energy_state(nktotf,nbndfst,H_nk,ee,pp)
    call convert_diabatic_adiabatic(pp,c_nk,ww)
    
    call random_number(flagr)
    flagd = 0.0d0
    do iefre = 1,nefre
      flagd = flagd + pp((init_ik-1)*nbndfst+init_ik,iefre)**2
      if(flagr <= flagd) then
        isurface = iefre
        exit
      endif
    enddo
    
  end subroutine init_dynamical_variable
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% convert wavefunction from diabatix to adiabatic basis %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  subroutine convert_diabatic_adiabatic(pp,c_nk,ww)
    implicit none
    real(kind=dp),intent(in) :: pp(nefre,nefre)
    complex(kind=dpc),intent(in) :: c_nk(nbndfst,nktotf)
    complex(kind=dpc),intent(out):: ww(nefre)
    
    integer :: iefre,jefre,ik,jk,iband,jband 
    
    ww=0.0d0
    do ik=1,nktotf
      do iband=1,nbndfst
        iefre = (ik-1)*nbndfst + iband
        do jk=1,nktotf
          do jband=1,nbndfst
            jefre = (jk-1)*nbndfst + jband 
            ww(iefre) = ww(iefre)+pp(jefre,iefre)*c_nk(jband,jk)
          enddo
        enddo
      enddo
    enddo
  
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
    do iefre=1,nefre
      do jefre=1,iefre
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
          enddo
        enddo
        dd(iefre,jefre,:,:) = dd(iefre,jefre,:,:)/(ee(jefre)-ee(iefre))
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
 
 

end module surfacehopping