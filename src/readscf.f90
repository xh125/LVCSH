module readscf
  use kinds,only : dp
  implicit none
  logical :: print_ef_scf = .false.
  real(kind=dp) :: ef_scf_up,ef_scf_dw  
  contains
  
  subroutine readpwscf_out(pwscfout_name)
    use constants,only : maxlen,ryd2ev
    use io,only : stdout,io_file_unit,open_file,close_file,findkword,findkline,io_error
    use kinds,only : dp
    use cell_base,only : alat, ibrav, omega,celldm,at,bg,wmass
    use ions_base,only : nat,iat, ntyp,ityp,atm,zv,amass,iatm,tau,xau
    use symm_base,only : nsym,isym,sname,s, sr
    use klist,only : nelec,nelup,neldw,two_fermi_energies,nkstot,lgauss,ltetra,xk,xkg,&
                     wk,tot_charge
    use gvec, only : ecutrho,ecutwfc,dft_is_hybrid,ecutfock,ecfixed,qcutz,q2sigma,&
                    doublegrid,ngm_g,ngms_g
    use relax, only: epse,epsf,epsp
    use funct, only: iexch, icorr, igcx, igcc, inlc, imeta, imetac,exx_fraction,dft
    use pw_control_flags, only : nstep,nmix,imix,tr2,mixing_beta,mixing_style,lscf,conv_elec,&
                                  lbands,llondon,ldftd3,lxdm,ts_vdw,textfor,iverbosity,&
                                  conv_elec
    use ener,only : etot, hwf_energy, eband, deband, ehart, &
                    vtxc, etxc, etxcc, ewld, demet, epaw, &
                    elondon, edftd3, ef_up, ef_dw, exdm,ef
    USE tsvdw_module,         ONLY : EtsvdW 
    use spin_orb,only :lspinorb,domag,lforcet
    use noncolin_module,only : noncolin
    use vlocal,only : starting_charge
    use lsda_mod,only : lsda,starting_magnetization
    use wvfct,only : nbnd,et,wg
    use fixed_occ,only : one_atom_occupations
    USE paw_variables,ONLY : okpaw, ddd_paw, total_core_energy, only_paw 
    USE extfield,             ONLY : tefield, etotefield,gate, etotgatefield
    
    implicit none
    character(len=*),intent(in) :: pwscfout_name
    integer :: pwscfout_unit
    logical :: alive
    character(len=maxlen) :: pwscf_info
    integer :: nk_
    character(len=maxlen) :: ctmp
    integer :: itmp
    logical :: lfindkword 
    integer :: i,ipol,apol,nt,ik,ibnd
    REAL(DP) :: xkg_(3)
    
    real(kind=dp) :: ef_scf = 0.0,ehomo,elumo
    REAL(DP) :: eext=0.0_DP
    !! external forces contribution to the total energy
    REAL(DP) :: dr2
    !! the norm of the diffence between potential    
    !LOGICAL :: tefield
    !! if .TRUE. a finite electric field is added to the
    !! local potential
    !LOGICAL :: dipfield
    !! if .TRUE. the dipole field is subtracted
    ! TB
    
    ! coordinates of the k point in crystal axes
    ! ... a few local variables
    !  
    INTEGER, ALLOCATABLE :: &
        ngk_g(:)       ! number of plane waves summed on all nodes
    logical :: lda_plus_u
    real(kind=dp) :: eth
    REAL(DP) :: descf
    !! correction for variational energy    
    
    pwscfout_unit = io_file_unit()
    
    
    inquire(file=trim(adjustl(pwscfout_name)),exist=alive)
    if(.NOT. alive) then
      call io_error("Error:PWscf output file "//trim(adjustl(pwscfout_name))//" doesn't exist.")
    else
      call open_file(pwscfout_name,pwscfout_unit)
    endif    

    
    
    call findkword(pwscfout_unit,"Program")
    read(pwscfout_unit,"(A)") pwscf_info
    write(stdout,*) pwscf_info


    !100 FORMAT( /,/,5X, &
    ! &     'bravais-lattice index     = ',I12,/,5X, &
    ! &     'lattice parameter (alat)  = ',F12.4,'  a.u.',/,5X, &
    ! &     'unit-cell volume          = ',F12.4,' (a.u.)^3',/,5X, &
    ! &     'number of atoms/cell      = ',I12,/,5X, &
    ! &     'number of atomic types    = ',I12)    
    call findkword(pwscfout_unit,"bravais-lattice")
    read(pwscfout_unit,"(33X,I12)")   ibrav
    read(pwscfout_unit,"(33X,F12.4)") alat   ! in unit a.u.
    read(pwscfout_unit,"(33X,F12.4)") omega  ! in unit (a.u.)^3
    read(pwscfout_unit,"(33X,I12)")   nat
    read(pwscfout_unit,"(33X,I12)")   ntyp


    !101 FORMAT(5X, &
    !   &     'number of electrons       = ',F12.2,' (up:',f7.2,', down:',f7.2,')')    
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(len(trim(adjustl(ctmp)))==66) then
      !101
      two_fermi_energies = .true.
      read(pwscfout_unit,"(33X,F12.2,5X,f7.2,7X,f7.2)") nelec,nelup,neldw
    else
      !102
      two_fermi_energies = .false.
      read(pwscfout_unit,"(33X,F12.2)") nelec
    endif
    !103
    read(pwscfout_unit,"(33X,I12)") nbnd
    read(pwscfout_unit,"(33X,F12.4)") ecutwfc
    read(pwscfout_unit,"(33X,F12.4)") ecutrho
    
    !104
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:32)=="cutoff for Fock operator  =") then
      dft_is_hybrid = .true.
      read(pwscfout_unit,"(33X,F12.4)") ecutfock
    endif
    
    !105
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:32)=="scf convergence threshold =") then
      lscf = .true.
      read(pwscfout_unit,"(33X,1PE12.1)") tr2
      read(pwscfout_unit,"(33X,0PE12.4)") mixing_beta
      read(pwscfout_unit,"(33X,I12,2X,A)") nmix,mixing_style
      select case(trim(adjustl(mixing_style)))
        case('plain')
          imix = 0
        case('TF')
          imix = 1
        case('local-TF')
          imix = 2
      end select    
      
    endif
    
    !106
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:32)=="energy convergence thresh.=") then
      read(pwscfout_unit,"(33X,1PE12.1)") epse
      read(pwscfout_unit,"(33X,1PE12.1)") epsf
    endif    

    !107
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:32)=="press convergence thresh. =") then
      read(pwscfout_unit,"(33X,1PE12.1)") epsp
    endif        
    
    !call write_dft_name ( ) 
    read(pwscfout_unit,"(27X,A)") dft
    read(pwscfout_unit,"(28X,7I4)") iexch, icorr, igcx, igcc, inlc, imeta, imetac
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:32)=="EXX-fraction              =") then
      read(pwscfout_unit,"(33X,1PE12.1)") exx_fraction
    endif           
    
    !IF ( lmd .OR. lbfgs )
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:32)=="nstep                     =") then
      read(pwscfout_unit,"(33X,I12)") nstep
    endif
    
    read(pwscfout_unit,"(A)") ctmp    
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
      backspace(unit=pwscfout_unit)
    endif
    
    ! IF ( qcutz > 0.D0 ) THEN    
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(ctmp(6:49)=="A smooth kinetic-energy cutoff is imposed at") then
      read(pwscfout_unit,"(50X,f12.4)") ecfixed
      read(pwscfout_unit,"(41X,f21.4)") qcutz
      read(pwscfout_unit,"(41X,f21.4)") q2sigma
    endif
    

    !
    ! ... and here more detailed information. Description of the unit cell
    !
    call findkword(pwscfout_unit,"celldm(1)=")   
    read(pwscfout_unit,"(2(3X,3(12X,F11.6),/))") (celldm(i),i=1,6)   
    read(pwscfout_unit,"(/,3(24X,3F11.6,/))") ((at(ipol,apol),ipol=1,3),apol=1,3)
    !cart. coord. in units of alat
    read(pwscfout_unit,"(/,3(24X,3F10.6,/))") ((bg(ipol,apol),ipol=1,3),apol=1,3)    
    !cart. coord. in units 2 pi/alat
    
    !  CALL print_ps_info ( )
    
    !WRITE( stdout, '(/5x, "atomic species   valence    mass     pseudopotential")')
    call findkword(pwscfout_unit,"atomic")
    read(pwscfout_unit,*)
    do nt =1 ,ntyp
      read(pwscfout_unit,"(5x,a6,6x,f10.2,2x,f10.5)") atm(nt),zv(nt),amass(nt)
    enddo
    
    !  IF (calc.EQ.'cd' .OR. calc.EQ.'cm' )
    read(pwscfout_unit,"(/5X,A)") ctmp
    backspace(unit=pwscfout_unit)
    backspace(unit=pwscfout_unit)
    if(ctmp(1:12)==" cell mass =") then
      read(pwscfout_unit,"(/5x,12x,f10.5)") wmass
      if(ctmp(23:36)==" AMU/(a.u.)^2 ") wmass = wmass
    endif

    read(pwscfout_unit,"(/5X,A)") ctmp
    backspace(unit=pwscfout_unit)
    backspace(unit=pwscfout_unit)
    if(ctmp(1:26)=="Starting charge structure ") then
      read(pwscfout_unit,"(///,A)") ctmp
      DO nt = 1, ntyp
        read( pwscfout_unit, '(5x,a6,9x,f6.3)') atm(nt), starting_charge(nt)
      ENDDO      
    endif    

    !IF (lsda) THEN
    read(pwscfout_unit,"(/5X,A)") ctmp
    backspace(unit=pwscfout_unit)
    backspace(unit=pwscfout_unit)
    if(ctmp(1:28)=="Starting magnetic structure ") then
      lsda = .true.
      read(pwscfout_unit,"(///,A)") ctmp
      DO nt = 1, ntyp
        WRITE( pwscfout_unit, '(5x,a6,9x,f6.3)') atm(nt), starting_magnetization(nt)
      ENDDO      
    endif        
    
    !!
    !!   description of symmetries
    !!
    !!CALL  print_symmetries ( iverbosity, noncolin, domag )        
    !call findkline(pwscfout_unit,"Sym. Ops.",9,17)
    !read(pwscfout_unit,"(5X,i2)") nsym
    !call findkword(pwscfout_unit,"isym")
    !backspace(pwscfout_unit)
    !do isym=1 , nsym
    !  read(pwscfout_unit,"(/20X,a45/)") sname(isym)
    !  read(pwscfout_unit,"(19X,3(i6,5X))") (s(1,ipol,isym),ipol=1,3)
    !  read(pwscfout_unit,"(19X,3(i6,5X))") (s(2,ipol,isym),ipol=1,3)
    !  read(pwscfout_unit,"(19X,3(i6,5X)/)") (s(3,ipol,isym),ipol=1,3)
    !  read(pwscfout_unit,"(19X,3F11.7)") (sr(1,ipol,isym),ipol=1,3)
    !  read(pwscfout_unit,"(19X,3F11.7)") (sr(2,ipol,isym),ipol=1,3)
    !  read(pwscfout_unit,"(19X,3F11.7/)") (sr(3,ipol,isym),ipol=1,3)
    !  
    !enddo 
    
    
    call findkword(pwscfout_unit,"Cartesian")
    allocate(ityp(nat))
    allocate(iatm(nat))
    allocate(tau(3,nat))
    ityp = 0
    tau  = 0.0
    
    read(pwscfout_unit,"(//)")
    !positions (alat units)
    do iat=1,nat
      read(pwscfout_unit,"(21X,A3,15X,3F12.7)") iatm(iat),(tau(ipol,iat),ipol=1,3)
    enddo 
    
    !positions (cryst. coord.)
    allocate(xau(3,nat))
    xau = 0.0    
    read(pwscfout_unit,"(/,A)") ctmp
    if(trim(adjustl(ctmp))=="Crystallographic axes") then
    ! iverbosity > 0
      read(pwscfout_unit,"(//)")
      do iat=1,nat
        read(pwscfout_unit,"(21X,A3,15X,3F11.7)") iatm(iat),(xau(ipol,iat),ipol=1,3)
      enddo 
    else
      backspace(unit=pwscfout_unit)
      backspace(unit=pwscfout_unit)
      do iat=1,nat
        do ipol=1,3
          xau(ipol,iat)=bg(1,ipol)*tau(1,iat)+bg(2,ipol)*tau(2,iat)+bg(3,ipol)*tau(3,iat)
        enddo 
      enddo
    endif
       
    !call findkline(pwscfout_unit,"number of k points=",6,24)
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit=pwscfout_unit)
    if(len(trim(adjustl(ctmp)))==25) then
      lgauss = .false.
      ltetra = .false.
      read(pwscfout_unit,"(24X,I6)") nk_
    elseif(ctmp(31:51)==" (tetrahedron method)") then
      ltetra = .true.
      read(pwscfout_unit,"(24X,I6)") nk_
    else
      lgauss = .true.
      read(pwscfout_unit,"(24X,I6)") nk_
    endif
    
    IF ( lsda ) THEN
      !
      ! ... LSDA case: do not print replicated k-points
      !
      nkstot = 2*nk_
      !nk_ = nkstot/2
    ELSE
      nkstot = nk_
      !nk_ = nkstot
    END IF      
    allocate(xk(3,nk_),wk(nk_))
    allocate(xkg(3,nk_))
    xk = 0.0
    wk = 0.0
    
    read(pwscfout_unit,"(A)") ctmp
    if(trim(adjustl(ctmp))=="cart. coord. in units 2pi/alat") then
      ! iverbosity > 0 .OR. nk_ < 100
      do ik=1,nk_
        read(pwscfout_unit,"(20X,3F12.7,7X,F12.7)") (xk(ipol,ik),ipol=1,3),wk(ik)
      enddo
    else  
      read(pwscfout_unit,"(A)") ctmp
      if(trim(adjustl(ctmp))=="Number of k-points >= 100: set verbosity='high' to print them.") then
        write(stdout,"(A)") ctmp
      endif
    endif
    
    read(pwscfout_unit,"(/,A)") ctmp
    if(trim(adjustl(ctmp))=="cryst. coord.") then
      !iverbosity>0
      do ik=1,nk_
        read(pwscfout_unit,"(20x,3f12.7,7x,f12.7)") (xkg(ipol,ik),ipol=1,3),wk(ik)
      enddo
    else
      backspace(unit=pwscfout_unit)
      backspace(unit=pwscfout_unit)
    endif
    
    read(pwscfout_unit,"(/5x,13x,i8)") ngm_g
    
    read(pwscfout_unit,"(/,A)") ctmp
    backspace(unit=pwscfout_unit)
    backspace(unit=pwscfout_unit)
    if(ctmp(6:17)=="Smooth grid:") then
      doublegrid = .true.
      read(pwscfout_unit,"(/5x,13x,i8)") ngms_g
    endif
    
    
    !CALL create_scf_type( rhoin )
    call findkline(pwscfout_unit,"Self-consistent Calculation",6,32)
    
    call findkline(pwscfout_unit,"End of self-consistent calculation",6,39) ! line 847 of electrons.f90
    conv_elec = .true.
    
    ! call print_ks_engergies()
    read(pwscfout_unit,"(/,5x,a)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp=="Number of k-points >= 100: set verbosity='high' to print the bands.") then
      write(stdout,*) "Don't print_ks_engergies for the reason of nkstot >=100 and iverbosity <=0"
      write(stdout,*) "need to set iverbosity = 1 in scf.in"
    else
      allocate ( ngk_g(nkstot))
      read(pwscfout_unit,"(/,A)") ctmp
      if(ctmp == "------") then
        lforcet = .true.
        read(pwscfout_unit,"(///)")
      else
        lforcet = .false.
        backspace(unit=pwscfout_unit)
        backspace(unit=pwscfout_unit)
      endif                  
    
      allocate(et(nbnd,nk_),wg(nbnd,nk_))
    
      do ik=1,nkstot
        if (lsda) read(pwscfout_unit,*) 
        read(pwscfout_unit,"(A)") ctmp
        backspace(unit=pwscfout_unit)
        if(ctmp(36:37)=="(") then
          conv_elec = .true.
          read(pwscfout_unit,"(/,13X,3F7.4,2X,I6)") (xk(i,ik),i=1,3),ngk_g(ik)
          ! FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
          ! FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   bands (ev):'/ )
        else
          conv_elec = .false.
          read(pwscfout_unit,"(/,13X,3F7.4)") (xk(i,ik),i=1,3)
        endif
        read(pwscfout_unit,"(2X,8F9.4)") (et(ibnd,ik),ibnd=1,nbnd)
        read(pwscfout_unit,"(/,A)") ctmp
        if(ctmp == "     occupation numbers ") then
          lbands = .false.
          read(pwscfout_unit,"(2X,8F9.4)") (wg(ibnd,ik),ibnd=1,nbnd)
        else
          backspace(unit=pwscfout_unit)
          backspace(unit=pwscfout_unit)
        endif
        
      enddo
    endif
    
    !CALL print_ks_ef_homolumo ( .false., 0.0_dp, 0.0_dp, 0.0_dp )
    if(.not. lbands) then !IF ( .NOT. lbands) 
      if(lgauss .or. ltetra) then
        if(two_fermi_energies) then
          read(pwscfout_unit,"(/,39X,2F10.4)") ef_up,ef_dw
          ef_up = ef_up/ryd2ev
          ef_dw = ef_dw/ryd2ev
          !FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev' )
          if(print_ef_scf) then
            read(pwscfout_unit,"(20X,2F10.4)") ef_scf_up,ef_scf_dw
            ef_scf_dw= ef_scf_dw/ryd2ev
            ef_scf_up= ef_scf_up/ryd2ev
            !FORMAT( '     (compare with: ',2F10.4,' eV, computed in scf)' )
          endif
        else
          read(pwscfout_unit,"(/,25X,f10.4)") ef
          ef = ef/ryd2ev
          !FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
          if(print_ef_scf) then
            read(pwscfout_unit,"(20X,F10.4)") ef_scf
            ef_scf = ef_scf/ryd2ev
            !FORMAT( '     (compare with: ', F10.4,' eV, computed in scf)' )
          endif
        endif
      else if(.not. one_atom_occupations) then
        ! ... presumably not a metal: print HOMO (and LUMO if available)
        read(pwscfout_unit,"(/,A)") ctmp
        backspace(unit=pwscfout_unit)
        if(ctmp(6:33)=="highest occupied level (ev):") then
          !FORMAT(/'     highest occupied level (ev): ',F10.4 )
          read(pwscfout_unit,"(34X,F10.4)") ehomo
          ehomo = ehomo/ryd2ev
        else
          read(pwscfout_unit,"(53X,2F10.4)")ehomo,elumo
          !FORMAT(/'     highest occupied, lowest unoccupied level (ev): ',2F10.4 ) 
          ehomo =ehomo/ryd2ev
          elumo = elumo/ryd2ev
        endif
      endif  
    endif
    
    !CALL print_energies ( printout )
    read(pwscfout_unit,"(/,32X,0PF17.8)") etot
    !9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' )
    read(pwscfout_unit,"(/,A)") ctmp
    !9085 FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )
    if(ctmp(6:32)=="total all-electron energy =") then
      only_paw = .true.
      backspace(unit=pwscfout_unit)
      read(pwscfout_unit,"(32X,0PF17.6)") total_core_energy
      total_core_energy = total_core_energy -etot
    else
      backspace(unit=pwscfout_unit)
      backspace(unit=pwscfout_unit)
    endif
    
    !9082 FORMAT( '     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
    read(pwscfout_unit,"(A)") ctmp
    if(ctmp(6:32)=="Harris-Foulkes estimate   =") then
      backspace(unit=pwscfout_unit)
      read(pwscfout_unit,"(32X,0PF17.8)") hwf_energy
    else
      backspace(unit=pwscfout_unit)
    endif
    
    !9083 FORMAT( '     estimated scf accuracy    <',0PF17.8,' Ry' )
    !9084 FORMAT( '     estimated scf accuracy    <',1PE17.1,' Ry' )
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(45:46)=="E") then
      read(pwscfout_unit,"(32X,1PE17.1)") dr2
    else
      read(pwscfout_unit,"(32X,0PF17.8)") dr2
    endif
    
    if(lgauss) then
      read(pwscfout_unit,"(32X,F17.8)") demet
      read(pwscfout_unit,"(32X,0PF17.8)") etot
      etot = etot+demet
      read(pwscfout_unit,*)
    else
      read(pwscfout_unit,*)
    endif
    
    !9062 FORMAT( '     one-electron contribution =',F17.8,' Ry' &
    !        /'     hartree contribution      =',F17.8,' Ry' &
    !        /'     xc contribution           =',F17.8,' Ry' &
    !        /'     ewald contribution        =',F17.8,' Ry' )
    read(pwscfout_unit,"(32X,/,32X,F17.8,32X,/,32X,F17.8)") ehart,ewld
    
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="Dispersion Correction     =") then
      llondon = .true.
      read(pwscfout_unit,"(32X,F17.8)") elondon
    endif
    
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="DFT-D3 Dispersion         =") then
      ldftd3 = .true.
      read(pwscfout_unit,"(32X,F17.8)") edftd3
    endif    
    
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="Dispersion XDM Correction =") then
      lxdm = .true.
      read(pwscfout_unit,"(32X,F17.8)") exdm
    endif        

    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="Dispersion T-S Correction =") then
      ts_vdw = .true.
      read(pwscfout_unit,"(32X,F17.8)") EtsvdW
      EtsvdW = EtsvdW/ 2.0
    endif   
    
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="External forces energy    =") then
      textfor = .true.
      read(pwscfout_unit,"(32X,F17.8)") eext
    endif         
    
    
    
    !IF ( llondon ) WRITE ( stdout , 9074 ) elondon
    !IF ( ldftd3 )  WRITE ( stdout , 9078 ) edftd3
    !IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
    !IF ( ts_vdw )  WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
    !IF ( textfor)  WRITE ( stdout , 9077 ) eext
    !IF ( tefield )            WRITE( stdout, 9064 ) etotefield
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="electric field correction =") then
      tefield = .true.
      read(pwscfout_unit,"(32X,F17.8)") etotefield
    endif    
    !IF ( gate )               WRITE( stdout, 9065 ) etotgatefield
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="gate field correction     =") then
      gate = .true.
      read(pwscfout_unit,"(32X,F17.8)") etotgatefield
    endif             
    !IF ( lda_plus_u )         WRITE( stdout, 9066 ) eth
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="Hubbard energy            =") then
       lda_plus_u= .true.
      read(pwscfout_unit,"(32X,F17.8)") eth
    endif        
    !IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf  
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="scf correction            =") then
   
      read(pwscfout_unit,"(32X,F17.8)") descf 
    endif        
    !9064 FORMAT( '     electric field correction =',F17.8,' Ry' )
    !9065 FORMAT( '     gate field correction     =',F17.8,' Ry' ) ! TB
    !9066 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
    !9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
    !9068 FORMAT( '      -> PAW hartree energy AE =',F17.8,' Ry' &
    !        /'      -> PAW hartree energy PS =',F17.8,' Ry' &
    !        /'      -> PAW xc energy AE      =',F17.8,' Ry' &
    !        /'      -> PAW xc energy PS      =',F17.8,' Ry' &
    !        /'      -> total E_H with PAW    =',F17.8,' Ry' &
    !        /'      -> total E_XC with PAW   =',F17.8,' Ry' )
    !9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
    !9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
    !9075 FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
    !9076 FORMAT( '     Dispersion T-S Correction =',F17.8,' Ry' )
    !9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
    !9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="one-center paw contrib.   =") then
      okpaw = .true.
      read(pwscfout_unit,"(32X,F17.8)") epaw 
      if(iverbosity > 0 ) then
        read(pwscfout_unit,"(/////)")
      endif
    endif            
    
    !9072 FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
    read(pwscfout_unit,"(A)") ctmp
    backspace(unit = pwscfout_unit)
    if(ctmp(6:32)=="pot.stat. contrib. (-muN) =") then
      read(pwscfout_unit,"(32X,F17.8)") tot_charge 
      tot_charge = tot_charge/ef
    endif            
    
    call close_file(pwscfout_name,pwscfout_unit)
    write(stdout,*) "Read PWscf outfile Succesful."
    
  end subroutine readpwscf_out
  
end module readscf