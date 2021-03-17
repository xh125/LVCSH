module readscf
implicit none
  
  contains
  
  subroutine readpwscf_out(pwscfout_name)
    use constants,only : maxlen
    use io,only : stdout,io_file_unit,open_file,close_file,findkword,findkline,io_error
    use kinds,only : dp
    use cell_base,only : alat, ibrav, omega,celldm,at,bg,wmass
    use ions_base,only : nat,iat, ntyp,ityp,atm,zv,amass,iatm,tau,xau
    use symm_base,only : nsym,isym,sname,s, sr
    use klist,only : nelec,nelup,neldw,two_fermi_energies,nkstot,lgauss,ltetra,xk,wk
    use gvec, only : ecutrho,ecutwfc,dft_is_hybrid,ecutfock,ecfixed,qcutz,q2sigma,&
                    doublegrid,ngm_g,ngms_g
    use relax, only: epse,epsf,epsp
    use funct, only: iexch, icorr, igcx, igcc, inlc, imeta, imetac,exx_fraction,dft
    use pw_control_flags, only : nstep,nmix,imix,tr2,mixing_beta,mixing_style,lscf
    use spin_orb,only :lspinorb,domag
    use noncolin_module,only : noncolin
    use vlocal,only : starting_charge
    use lsda_mod,only : lsda,starting_magnetization
    use wvfct,only : nbnd,et,wg
    
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
    ! coordinates of the k point in crystal axes
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
    
    call findkword(pwscfout_unit,"bravais-lattice")
    read(pwscfout_unit,"(33X,I12)")   ibrav
    read(pwscfout_unit,"(33X,F12.4)") alat   ! in unit a.u.
    read(pwscfout_unit,"(33X,F12.4)") omega  ! in unit (a.u.)^3
    read(pwscfout_unit,"(33X,I12)")   nat
    read(pwscfout_unit,"(33X,I12)")   ntyp

    call findkline(pwscfout_unit,"number of electrons       =",6,32)
    read(pwscfout_unit,"(A)") ctmp
    itmp =len(trim(adjustl(ctmp)))
    backspace(unit=pwscfout_unit)
    if(itmp> 40) then
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
        WRITE( pwscfout_unit, '(5x,a6,9x,f6.3)') atm(nt), starting_charge(nt)
      ENDDO      
    endif    

    !IF (lsda) THEN
    read(pwscfout_unit,"(/5X,A)") ctmp
    backspace(unit=pwscfout_unit)
    backspace(unit=pwscfout_unit)
    if(ctmp(1:28)=="Starting magnetic structure ") then
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
    do iat=1,nat
      do ipol=1,3
        xau(ipol,iat)=bg(1,ipol)*tau(1,iat)+bg(2,ipol)*tau(2,iat)+bg(3,ipol)*tau(3,iat)
      enddo 
    enddo
    
    call findkline(pwscfout_unit,"number of k points=",6,24)
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
    xk = 0.0
    wk = 0.0
    
    read(pwscfout_unit,"(A)") ctmp
    if(trim(adjustl(ctmp))=="cart. coord. in units 2pi/alat") then
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
        read(pwscfout_unit,"(20x,3f12.7,7x,f12.7)") (xkg_(ipol),ipol=1,3),wk(ik)
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
      read(pwscfout_unit,"(/5x,12x,i8)") ngms_g
    endif
    
    
    
    allocate(et(nbnd,nk_),wg(nbnd,nk_))
    call findkline(pwscfout_unit,"          k =",1,13)
    backspace(pwscfout_unit)
    do ik=1,nk_
      read(pwscfout_unit,"(/,A,/)") ctmp
      read(pwscfout_unit,"(2X,8F9.4)") (et(ibnd,ik),ibnd=1,nbnd)
      read(pwscfout_unit,"(A,/)") ctmp
      read(pwscfout_unit,"(2X,8F9.4)") (wg(ibnd,ik),ibnd=1,nbnd)
    enddo
    
    
    
    call close_file(pwscfout_name,pwscfout_unit)
    
  end subroutine readpwscf_out
  
end module readscf