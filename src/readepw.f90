module readepw
  use kinds ,only :dp
  use io, only : io_file_unit,open_file,close_file,findkword,findkline,stdout
  use klist, only : nkstot,xk,wk 
  use klist_epw, only : xk_all
  use epwcom, only : nkc1,nkc2,nkc3,nqf1,nqf2,nqf3,nkf1,nkf2,nkf3,nbndsub,kqmap
  use elph2, only : nkqtotf,wf,epcq,nktotf,nqtotf,ibndmin,ibndmax
  implicit none
  
  integer ::  lower_bnd,           &!
              upper_bnd,           &!
              ik       ,           &! K-point index
              ikk      ,           &! K-point index
              ikq      ,           &! K+q-point index
              ibnd     ,           &! Band index
              jbnd     ,           &! Band index
              pbnd     ,           &! Band index
              mu       ,           &! Mode index
              nu       ,           &! Mode index
              totq     ,           &! Total number of q-points within the fsthick window.
              iq,iq_                ! Current q-point index
              
  real(kind=dp), allocatable :: xkf_all(:,:)     
  !  Collect k-point coordinate from all pools in parallel case
  real(kind=dp), allocatable :: etf_all(:,:)
  !  Collect eigenenergies from all pools in parallel case
  real(kind=dp) ::  wq          ,&! Phonon frequency
                    ekk         ,&! Eigenenergies at k
                    ekq           ! Eigenenergies at k+q
  
  real(kind=dp),allocatable :: epc(:,:,:,:)
  !! g vectex accross all pools
  real(kind=dp),allocatable :: epc_sym(:,:,:)
  !! Temporary g-vertex for each pool
  !real(kind=dp),allocatable :: epmatq(:,:,:,:,:)
  real(kind=dp) :: ebndmax,ebndmin
  real(kind=dp),allocatable :: Enk(:,:),Emkq(:,:)
  
  contains
  
  subroutine readepwout(fepwout)
    use elph2,    only : nqtotf,xqf,nbndfst
    use grid,     only : loadqmesh,kq2k_map,loadkmesh_fullBZ
    use modes,    only : nmodes
    implicit none
    character(len=*) ,intent(in) :: fepwout
    integer :: unitepwout
    integer :: ipol,ik_,ibnd_,jbnd_,nu_
    real(kind=dp) :: xiq(3)
    !real(kind=dp) :: enk,emkq
    integer :: ierr
    !! error status
    
    unitepwout = io_file_unit()
    call open_file(fepwout,unitepwout)
    
    call findkline(unitepwout,"number of k points=",6,24)
    read(unitepwout,"(24X,I5)") nkstot
    !if(.not. allocated(xk_all)) allocate(xk_all(3,nkstot))
    !if(.not. allocated(wk)) allocate(wk(nkstot))
    allocate(xk_all(3,nkstot),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating xk_all',1)
    allocate(wk(nkstot),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating wk(nkstot)',1)    
    
    read(unitepwout,*)
    do ik =1,nkstot
      read(unitepwout,"(20X,3f12.7,7X,f12.7)") (xk_all(ipol,ik),ipol=1,3),wk(ik)
    enddo
    
    call findkline(unitepwout,"     Wannierization on ",1,23) 
    read(unitepwout,"(23X,3(i2,3X))") nkc1,nkc2,nkc3
    
    !WRITE(stdout, '("      - Number of total bands is (", i3, ")")') nbnd
    call findkline(unitepwout,"      - Number of total bands is (",1,34)
    read(unitepwout,"(34X,i3)") nbndsub
    
    call findkline(unitepwout,"     Using uniform q-mesh: ",1,27)
    read(unitepwout,"(27X,3i4)") nqf1,nqf2,nqf3
    nqtotf = nqf1 * nqf2 * nqf3
    
    !! set xqf(3,nqtotf) wqf(nqtotf)
    call loadqmesh()      
    
    !WRITE(stdout, '(5x,"Size of q point mesh for interpolation: ",i10)') nqtotf
    read(unitepwout,"(45X,i10)") nqtotf
    
    
    read(unitepwout,"(27X,3i4)") nkf1,nkf2,nkf3
    nktotf = nkf1 * nkf2 * nkf3
    nkqtotf = 2 * nktotf
    allocate(xkf_all(3,nkqtotf),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating xkf_all',1)

    allocate(wf(nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating wf',1)
    !allocate(epmatq(nbndfst,nbndfst,nktotf,nmodes,nqtotf))
    
    read(unitepwout,"(45X,i10)") nkqtotf
    
    ! set xkf_bz(3,nktotf)
    call loadkmesh_fullBZ()
    call kq2k_map()
    
    
    
    call findkword(unitepwout,"ibndmin")
    read(unitepwout,"(14X,10X,i5,2x,10X,f9.3)") ibndmin, ebndmin
    read(unitepwout,"(14X,10X,i5,2x,10X,f9.3)") ibndmax, ebndmax
    nbndfst = ibndmax-ibndmin + 1
    allocate(etf_all(nbndfst,nkqtotf),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating etf_all',1)
    allocate(Enk(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore("readepw",'Error allocating Enk',1)
    allocate(Emkq(nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore("readepw",'Error allocating Emkq',1)
    allocate(epcq(nbndfst,nbndfst,nktotf,nmodes,nqtotf),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating epmatq',1)
  
    call findkline(unitepwout,"We only need to compute",6,28)
    read(unitepwout,"(29X,i8)") totq  ! totq = nqf
    
    

    do iq=1,nqtotf
      !WRITE(stdout, '(5x, a)') ' Electron-phonon vertex |g| (meV)'   printing.f90
      call findkline(unitepwout," Electron-phonon vertex |g| (meV)",6,38)  
      read(unitepwout,"(//,10x,i7,9x, 3f12.7)") iq_,(xiq(ipol),ipol=1,3) !(xqf(ipol,iq),ipol=1,3)
      if((xiq(1) /= xqf(1,iq)) .or. (xiq(2) /= xqf(2,iq)) .or. (xiq(3) /= xqf(3,iq)) ) then
        write(stdout,*) "Warning: xqf set is wrong"
      endif
      xqf(:,iq) = xiq
      do ik=1,nktotf
        ikk = 2*ik-1
        ikq = ikk + 1
        read(unitepwout,'(5x,5x,i7, 9x, 3f12.7)') ik_, xkf_all(:, ikk)
        read(unitepwout,*)
        read(unitepwout,*)
        do ibnd = 1,nbndfst
          do jbnd = 1, nbndfst
            do nu = 1, nmodes
              !ekq = etf_all(ibndmin - 1 + jbnd, ikq)
              !ekk = etf_all(ibndmin - 1 + ibnd, ikk)
              !WRITE(stdout, '(3i9, 2f12.4, 1f20.10, 1e20.10)') ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, &
              !nu, ryd2ev * ekk, ryd2ev * ekq, ryd2mev * wf(nu, iq), ryd2mev * epc(ibnd, jbnd, nu, ik)
              read(unitepwout,'(3i9, 2f12.4, 1f20.10, 1e20.10)') ibnd_,jbnd_,nu_,&
                   etf_all(ibnd,ikk),etf_all(jbnd,ikq),wf(nu,iq),epcq(ibnd,jbnd,ik,nu,iq)
              Enk(ibnd,ik) = etf_all(ibnd,ikk)
              Emkq(jbnd,kqmap(ik,iq)) = etf_all(jbnd,ikq)
              !read(unitepwout,'(3i9, 2f12.4, 1f20.10, 1e20.10)') ibnd_,jbnd_,nu_,&
                   !enk,emkq,wf(nu,iq),epcq(ibnd,jbnd,ik,nu,iq)
            enddo
          enddo
        enddo
        read(unitepwout,"(/)")
      enddo

    enddo
    
    do ik=1,nktotf
      do ibnd=1,nbndfst
        if (Enk(ibnd,ik) /= Emkq(ibnd,ik)) then
          write(stdout,*) "Enk /= Emkq","ibnd=",ibnd,"ik=",ik
        endif
      enddo
    enddo
    
  end subroutine readepwout
  

end module readepw