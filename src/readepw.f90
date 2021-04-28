module readepw
  use kinds ,only :dp
  use constants,only : maxlen,amu_ry,rytoev,ryd2mev
  use io, only : io_file_unit,open_file,close_file,findkword,findkline,stdout,io_error
  use klist, only : lgauss, degauss, ngauss, nkstot, wk
  use klist_epw, only : xk_all
  use epwcom, only : nkc1,nkc2,nkc3,nqf1,nqf2,nqf3,nkf1,nkf2,nkf3,nbndsub,kqmap,scdm_proj,vme
  use pwcom, only : ef
  use elph2, only : nkqf,nkqtotf,wf,wqf,xkf,wkf,etf,epcq,nkf,&
                    nktotf,nqf,nqtotf,ibndmin,ibndmax,efnew,vmef,dmef
  use cell_base,only : ibrav,alat,omega,at,bg,celldm
  use ions_base,only : nat,iat, ntyp,ityp,atm,zv,amass,iatm,tau,iamass 
  use symm_base,only : nrot,nsym
  use funct, only: iexch, icorr, igcx, igcc, inlc, imeta, imetac,exx_fraction,dft
  use spin_orb,only :lspinorb,domag
  use noncolin_module,only : noncolin
  use gvec, only : ecutrho,ecutwfc
  use wannierEPW,only : num_bands,n_wannier,center_w,l_w,mr_w,nexband
  use wvfct,only : nbnd
  implicit none
  
  character(len=maxlen) :: epw_info,ctmp
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
  INTEGER :: i,j,k,iw
  !! generic counter
  INTEGER :: ipol
  !! counter on polarizations
  INTEGER :: apol
  !! counter on polarizations
  INTEGER :: na
  !! counter on atoms
  
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
  real(kind=dp),allocatable :: E_nk(:,:),E_mkq(:,:)
  
  contains
  
  subroutine readepwout(fepwout)
    use elph2,    only : nqtotf,xqf,nbndfst
    use grid,     only : loadqmesh,kq2k_map,loadkmesh_fullBZ,get_ikq
    use modes,    only : nmodes
    use lr_symm_base,only : l_nsymq_le_1,minus_q,nsymq
    implicit none
    character(len=*) ,intent(in) :: fepwout
    integer :: unitepwout
    integer :: ipol,ik_,ibnd_,jbnd_,nu_
    real(kind=dp) :: epc_
    real(kind=dp) :: xiq(3),xik(3)
    !real(kind=dp) :: E_nk,E_mkq
    integer :: ierr
    !! error status
    logical :: alive
    
    inquire(file=trim(adjustl(fepwout)),exist=alive)
    if(.NOT. alive) then
      call io_error("Error:EPW output file "//trim(adjustl(fepwout))//" doesn't exist.")
    else   
      unitepwout = io_file_unit()
      call open_file(fepwout,unitepwout)
    endif
    
    call findkword(unitepwout,"Program")
    read(unitepwout,"(A)") epw_info
    write(stdout,*) epw_info
    
    call findkline(unitepwout,"bravais-lattice index     =",6,32)
    read(unitepwout,"(33X,i12)")   ibrav
    read(unitepwout,"(33X,f12.4)") alat
    read(unitepwout,"(33X,f12.4)") omega
    read(unitepwout,"(33X,i12)")   nat
    read(unitepwout,"(33X,i12)")   ntyp
    read(unitepwout,"(33X,f12.4)") ecutwfc   
    read(unitepwout,"(33X,f12.4)") ecutrho !ecutwfc * dual
    
    nmodes = 3*nat
    
    
    !call write_dft_name ( ) 
    read(unitepwout,"(27X,A)") dft
    read(unitepwout,"(28X,7I4)") iexch, icorr, igcx, igcc, inlc, imeta, imetac
    read(unitepwout,"(A)") ctmp
    backspace(unit=unitepwout)
    if(ctmp(6:32)=="EXX-fraction              =") then
      read(unitepwout,"(33X,1PE12.1)") exx_fraction
    endif               
    
    !  !  Here add a message if this is a noncollinear or a spin_orbit calculation
    !epw_summary.f90 line 79
    !IF (noncolin) THEN
    !  IF (lspinorb) THEN
    !    IF (domag) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    read(unitepwout,"(A)") ctmp    
    if(ctmp(6:45)=="Noncollinear calculation with spin-orbit") then
      domag = .true.
      lspinorb = .true.
      noncolin = .true.
    elseif(ctmp(6:45)=="Non magnetic calculation with spin-orbit") then
      domag = .false.
      lspinorb = .true.
      noncolin = .true.
    elseif(ctmp(6:48)=="Noncollinear calculation without spin-orbit") then
      lspinorb = .false.
      noncolin = .true.
    else
      noncolin = .false.
    endif    

    !
    ! Description of the unit cell
    !    
    read(unitepwout,*)
    read(unitepwout,"(2(3X,3(12X,F11.5),/))") (celldm(i),i=1,6)
    read(unitepwout,"(/,3(23X,3F8.4,/))") ((at(ipol,apol),ipol=1,3),apol=1,3)
    read(unitepwout,"(/,3(23X,3F8.4,/))") ((bg(ipol,apol),ipol=1,3),apol=1,3)
    read(unitepwout,"(/////)")

    !
    ! Description of the atoms inside the unit cell
    !    
    if(.not. allocated(iatm)) then 
      allocate(iatm(nat),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating iatm',1)    
    endif
    if(.not. allocated(iamass)) then 
      allocate(iamass(nat),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating iamass',1)    
    endif  
    if(.not. allocated(tau)) then 
      allocate(tau(3,nat),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating tau',1)    
    endif    
    do iat =1 ,nat
      read(unitepwout,"(17X,a3,F8.4,14X,3f11.5)") iatm(iat),iamass(iat),(tau(ipol,iat),ipol=1,3)
    enddo
    ! atoms mass in "Rydberg" atomic units
    iamass = iamass * amu_ry
    
    !
    ! Description of symmetries
    !
    read(unitepwout,*)
    read(unitepwout,"(A)") ctmp
    if(ctmp(6:17)=="No symmetry!") then
      l_nsymq_le_1 = .true.
      minus_q = .false.
    else
      if(ctmp(8:34)==" Sym.Ops. (with q -> -q+G )") then
        minus_q = .true.
        backspace(unitepwout)
        read(unitepwout,"(5X,i2)") nsymq
        nsymq = nsymq -1
      else
        minus_q = .false.
        backspace(unitepwout)
        read(unitepwout,"(5X,i2)") nsymq
      endif
    endif

    !IF (iverbosity == 1) THEN
    !WRITE(stdout, '(36x,"s",24x,"frac. trans.")')
    
    !
    !     Description of the reciprocal lattice vectors
    !    
    
    
    
    call findkline(unitepwout,"number of k points=",6,24)
    read(unitepwout,"(A)") ctmp
    backspace(unit=unitepwout)
    if(len(trim(adjustl(ctmp)))==24) then
      lgauss = .false.
      read(unitepwout,"(24X,I5)") nkstot
      !  !! total number of k points in pwscf calculation
    else
      lgauss = .true.
      read(unitepwout,"(24X,I5,23X,f8.4,14X,i3)") nkstot,degauss,ngauss
    endif
    if(.not. allocated(xk_all)) then 
      allocate(xk_all(3,nkstot),stat=ierr)
      !! List of all kpoints in cartesian coordinates
      !! cart. coord. in units 2pi/a_0
      if(ierr /=0) call errore('readepw','Error allocating xk_all',1)    
    endif
    
    if( allocated(wk)) deallocate(wk) 
    allocate(wk(nkstot),stat=ierr)
    !! weight of k points of nkstot
    if(ierr /=0) call errore('readepw','Error allocating wk',1)    
    
    read(unitepwout,"(A)") ctmp
    if(trim(adjustl(ctmp))=="cart. coord. in units 2pi/a_0") then
      do ik=1,nkstot
        read(unitepwout,"(20X,3f12.7,7X,f12.7)") &
             (xk_all(ipol, ik) , ipol = 1, 3), wk(ik)
      enddo
    endif
    
    
    call findkline(unitepwout,"     Wannierization on ",1,23) 
    read(unitepwout,"(23X,3(i2,3X))") nkc1,nkc2,nkc3
    !! kx,ky,kz sizes of the uniform electron coarse mesh to be used
    
    call findkline(unitepwout,"Initial Wannier",6,21)
    read(unitepwout,"(A)") ctmp
    
    if(len(trim(adjustl(ctmp)))==27) scdm_proj= .false.
    
    !WRITE(stdout, '(/, "      - Number of bands is (", i3, ")")') num_bands
    call findkline(unitepwout,"      - Number of bands is (",1,28)
    read(unitepwout,"(28X,i3)") num_bands
    read(unitepwout,"(34X,i3)") nbnd
    read(unitepwout,"(37X,i3)") nexband
    read(unitepwout,"(40X,i3)") n_wannier
    nbndsub = n_wannier
    
    if((nbnd-nexband)/=num_bands) &
    call errore('setup_nnkp', ' something wrong with num_bands', 1)
    
    
    allocate(center_w(3,n_wannier),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating center_w',1) 
    allocate(l_w(nbnd),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating l_w',1) 
    allocate(mr_w(nbnd),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating mr_w',1) 
    
    if (.not. scdm_proj) then
      do i=1,n_wannier+5
        backspace(unitepwout)
      enddo
      do iw=1,n_wannier
        read(unitepwout,"(6x,3f10.5,9x,i3,6x,i3)") &
                      (center_w(ipol,iw),ipol=1,3),l_w(iw),mr_w(iw)
      enddo
    endif
    
    call findkline(unitepwout,"Symmetries of Bravais lattice: ",6,36)
    read(unitepwout,"(36X,i3)") nrot
    call findkline(unitepwout,"Symmetries of crystal:         ",6,36)
    read(unitepwout,"(36X,i3)") nsym        
    
    call findkline(unitepwout,"     Using uniform q-mesh: ",1,27)
    read(unitepwout,"(27X,3i4)") nqf1,nqf2,nqf3
    !!! qx,qy,qz sizes of the uniform phonon fine mesh to be used
    nqtotf = nqf1 * nqf2 * nqf3
    !  total number of q points (fine grid)
    if(.not. allocated(xqf)) then 
      allocate(xqf(3,nqtotf),stat=ierr)
      !  fine q point grid
      if(ierr /=0) call errore('readepw','Error allocating xqf',1)    
    endif        
    if(.not. allocated(wqf)) then 
      allocate(wqf(nqtotf),stat=ierr)
      !  weights on the fine q grid
      if(ierr /=0) call errore('readepw','Error allocating wqf',1)    
    endif
    wqf = 1.0d0/(dble(nqtotf))
    
    do i = 1 ,nqf1
      do j=1,nqf2
        do k=1,nqf3
          iq = (i - 1) * nqf2 * nqf3 + (j - 1) * nqf3 + k
          xqf(1, iq) = DBLE(i - 1) / DBLE(nqf1)
          xqf(2, iq) = DBLE(j - 1) / DBLE(nqf2)
          xqf(3, iq) = DBLE(k - 1) / DBLE(nqf3)          
        enddo
      enddo
    enddo  
    
    nqf = nqtotf
    !! set xqf(3,nqtotf) wqf(nqtotf)
    !call loadqmesh()      
    
    !WRITE(stdout, '(5x,"Size of q point mesh for interpolation: ",i10)') nqtotf
    read(unitepwout,"(45X,i10)") nqtotf
    
    
    read(unitepwout,"(27X,3i4)") nkf1,nkf2,nkf3
    nktotf = nkf1 * nkf2 * nkf3
    nkqtotf = 2 * nktotf
    if(.not. allocated(xkf)) then 
      allocate(xkf(3,nkqtotf),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating xkf',1)    
    endif        
    if(.not. allocated(wkf)) then 
      allocate(wkf(nkqtotf),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating wkf',1)    
    endif
    wkf = 0.0d0
    do ik=1,nktotf
      wkf(2*ik-1) = 2.0d0/dble(nkqtotf/2) !
    enddo
    DO i = 1, nkf1
      DO j = 1, nkf2
        DO k = 1, nkf3
          ik = (i - 1) * nkf2 * nkf3 + (j - 1) * nkf3 + k
          ikk = 2 * ik - 1
          ikq = ikk + 1
          xkf(1, ikk) = DBLE(i - 1) / DBLE(nkf1)
          xkf(2, ikk) = DBLE(j - 1) / DBLE(nkf2)
          xkf(3, ikk) = DBLE(k - 1) / DBLE(nkf3)
          xkf(1, ikq) = xkf(1, ikk)
          xkf(2, ikq) = xkf(2, ikk)
          xkf(3, ikq) = xkf(3, ikk)
        ENDDO
      ENDDO
    ENDDO       
    
    read(unitepwout,"(45X,i10)") nkqtotf    
    nkf = nkqtotf/2
    nkqf= nkqtotf

    
    
    call findkline(unitepwout,"Fermi energy coarse grid = ",6,32)
    read(unitepwout,"(32X,f10.6)") ef
    ef = ef /rytoev  
    call findkline(unitepwout,"Fermi energy is calculated from the fine k-mesh: Ef =",6,58)
    read(unitepwout,"(58X,f10.6)") efnew   
    efnew = efnew /rytoev
    ef = efnew
    
    ! set xkf_bz(3,nktotf)
    !call loadkmesh_fullBZ()
    !call kq2k_map()
    
    
    !icbm = 1
    !IF (noncolin) THEN
    !  icbm = FLOOR(nelec / 1.0d0) + 1
    !ELSE
    !  icbm = FLOOR(nelec / 2.0d0) + 1
    !ENDIF       
    call findkword(unitepwout,"ibndmin")
    read(unitepwout,"(14X,10X,i5,2x,10X,f9.3)") ibndmin, ebndmin
    ebndmin = ebndmin/rytoev
    read(unitepwout,"(14X,10X,i5,2x,10X,f9.3)") ibndmax, ebndmax
    ebndmax = ebndmax/rytoev
    
    nbndfst = ibndmax-ibndmin + 1
    
    
    allocate(etf(nbndsub,nkqf),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating etf',1)
    etf = 0.0d0
    
    allocate(E_nk(nbndfst,nkf),stat=ierr)
    if(ierr /=0) call errore("readepw",'Error allocating E_nk',1)
    allocate(E_mkq(nbndfst,nkf),stat=ierr)
    if(ierr /=0) call errore("readepw",'Error allocating E_mkq',1)
    allocate(epcq(nbndfst,nbndfst,nkf,nmodes,nqf),stat=ierr)
    if(ierr /=0) call errore('readepw','Error allocating epcq',1)
    
    if(.not. allocated(wf)) then 
      allocate(wf(nmodes,nkqtotf),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating wf',1)    
    endif            

    if(.not. allocated(kqmap)) then 
      allocate(kqmap(nkf,nqf),stat=ierr)
      if(ierr /=0) call errore('readepw','Error allocating kqmap',1)    
    endif              
    
  
    call findkline(unitepwout,"We only need to compute",6,28)
    read(unitepwout,"(29X,i8)") totq  ! totq = nqf
    nqf=totq
    
    !call findkline(unitepwout," Electron-phonon vertex |g| (meV)",6,38)     

    ! 1160 DO iqq = iq_restart, totq
    do iq=1,nqf
      !WRITE(stdout, '(5x, a)') ' Electron-phonon vertex |g| (meV)'   printing.f90
      call findkline(unitepwout," Electron-phonon vertex |g| (meV)",6,38)  
      read(unitepwout,"(//,10x,i7,9x, 3f12.7)") iq_,(xiq(ipol),ipol=1,3) !(xqf(ipol,iq),ipol=1,3)
      if((xiq(1) /= xqf(1,iq)) .or. (xiq(2) /= xqf(2,iq)) .or. (xiq(3) /= xqf(3,iq)) ) then
        write(stdout,*) "Warning: xqf set is wrong"
      endif
      !xqf(:,iq) = xiq
      
      do ik=1,nkf
        ikk = 2*ik-1
        ikq = ikk + 1
        read(unitepwout,'(5x,5x,i7, 9x, 3f12.7)') ik_, (xik(ipol),ipol=1,3)
        if((xik(1) /= xkf(1,ikk)) .or. (xik(2) /= xkf(2,ikk)) .or. (xik(3) /= xkf(3,ikk)) ) then
          write(stdout,*) "Warning: xkf set is wrong"
        endif
        kqmap(ik,iq) = get_ikq(xik,xiq)
        
        read(unitepwout,*)
        read(unitepwout,*)
        do ibnd = 1,nbndfst
          do jbnd = 1, nbndfst
            do nu = 1, nmodes
              !WRITE(stdout, '(5x, a)') ' ibnd     jbnd     imode   enk[eV]    enk+q[eV]  omega(q)[meV]   |g|[meV]'
              !ekq = etf_all(ibndmin - 1 + jbnd, ikq)
              !ekk = etf_all(ibndmin - 1 + ibnd, ikk)
              !WRITE(stdout, '(3i9, 2f12.4, 1f20.10, 1e20.10)') ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, &
              !nu, ryd2ev * ekk, ryd2ev * ekq, ryd2mev * wf(nu, iq), ryd2mev * epc(ibnd, jbnd, nu, ik)
              read(unitepwout,'(3i9, 2f12.4, 1f20.10, 1e20.10)') ibnd_,jbnd_,nu_,&
                   ekk,ekq,wf(nu,iq),epcq(ibnd,jbnd,ik,nu,iq)
              ekk = ekk /rytoev
              ekq = ekq /rytoev
              
              E_nk(ibnd,ik) = ekk
              E_mkq(jbnd,kqmap(ik,iq)) = ekq
              etf(ibndmin-1+ibnd,ikk) = ekk
              !read(unitepwout,'(3i9, 2f12.4, 1f20.10, 1e20.10)') ibnd_,jbnd_,nu_,&
                   !E_nk,E_mkq,wf(nu,iq),epcq(ibnd,jbnd,ik,nu,iq)
            enddo
          enddo
        enddo
        read(unitepwout,"(/)")
      enddo
      
      do nu=1,nmodes
        if(wf(nu,iq)<0.0) then
          wf(nu,iq)=-1.0*wf(nu,iq) 
          write(stdout,*) "Carefully!!! the energy of phonon in iq=",iq,"modes=",nu,"=",-1.0*wf(nu,iq)
        endif
      enddo
    enddo
    
    wf = wf/ryd2mev
    do nu=1,3
      epcq(:,:,:,nu,1) = 0.0
    enddo
    epcq = epcq/ryd2mev
    !do iq=1,nqf
    !  do nu=1,nmodes
    !    epcq(:,:,:,nu,iq) = epcq(:,:,:,nu,iq) * sqrt(2.0 *wf(nu,iq))
    !  enddo
    !enddo
    
    ! Testing wrong in reading Energy
    !do ik=1,nktotf
    !  do ibnd=1,nbndfst
    !    if (E_nk(ibnd,ik) /= E_mkq(ibnd,ik)) then
    !      write(stdout,*) "E_nk /= E_mkq","ibnd=",ibnd,"ik=",ik
    !    endif
    !  enddo
    !enddo
    
    call findkline(unitepwout,"matrix elements",15,29)
    read(unitepwout,"(A)") ctmp
    if(ctmp(32:35)=="vmef") then
      vme = .true.
    else
      vme= .false.
      !allocate(dmef(3,nbndsub,nbndsub,nkf))
    endif
    allocate(vmef(3,nbndsub,nbndsub,nkf))
    do ik=1,nkf
      read(unitepwout,"(//)")
      do ibnd=1,nbndsub
        do jbnd=1,nbndsub
          read(unitepwout,"(31X,3(2E16.6))") (vmef(ipol,ibnd,jbnd,ik),ipol=1,3)
        enddo
      enddo
    enddo
    ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k) 
    ! vmef = 2 *dmef
    ! ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)   
    if (.not. vme) vmef = 2.0*vmef  !in unit of Ryd*bohr   
    
    write(stdout,"(5X,A)") "Reading epw.out successfull!"
    
  end subroutine readepwout
  

end module readepw