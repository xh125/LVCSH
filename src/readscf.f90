module readscf
implicit none
  
  contains
  
  subroutine readpwscf_out(pwscfout_name)
    use constants,only : maxlen
    use io,only : stdout,io_file_unit,open_file,close_file,findkword,findkline
    use kinds,only : dp
    use cell_base
    use ions_base
    use symm_base,only : nsym,isym,sname,s, sr
    use kpoints
    use wavefct
    implicit none
    character(len=*),intent(in) :: pwscfout_name
    integer :: pwscfout_unit
    character(len=maxlen) :: pwscf_info
    character(len=maxlen) :: ctmp
    logical :: lfindkword 
    !integer :: ibrav , nat, ntyp
    !real(kind=dp) :: omega, alat
    integer :: i,ipol,apol,nt,ik,iband
    
    pwscfout_unit = io_file_unit()
    
    call open_file(pwscfout_name,pwscfout_unit)
    
    
    call findkword(pwscfout_unit,"Program")
    read(pwscfout_unit,"(A)") pwscf_info
    write(stdout,*) pwscf_info
    
    call findkword(pwscfout_unit,"bravais-lattice")
    read(pwscfout_unit,"(33X,I12)")   ibrav
    read(pwscfout_unit,"(33X,F12.4)") alat   ! in unit a.u.
    read(pwscfout_unit,"(33X,F12.4)") omega  ! in unit (a.u.)^3
    read(pwscfout_unit,"(33X,I12)")   nat
    read(pwscfout_unit,"(33X,I12)")   ntyp

    call findkline(pwscfout_unit,"     number of Kohn-Sham states=",1,32)
    read(pwscfout_unit,"(33X,I12)") nband
    
    call findkword(pwscfout_unit,"celldm(1)=")   
    read(pwscfout_unit,"(2(3X,3(12X,F11.6),/))") (celldm(i),i=1,6)   
    read(pwscfout_unit,"(/,3(24X,3F11.6,/))") ((at(ipol,apol),ipol=1,3),apol=1,3)
    !cart. coord. in units of alat
    read(pwscfout_unit,"(/,3(24X,3F10.6,/))") ((bg(ipol,apol),ipol=1,3),apol=1,3)    
    !cart. coord. in units 2 pi/alat
    
    call findkword(pwscfout_unit,"atomic")
    read(pwscfout_unit,"(/,3(5x,a6,6x,f10.2,2x,f10.5,/))") ((atm(nt),zv(nt),amass(nt)),nt=1,ntyp)
    
    call findkline(pwscfout_unit,"Sym. Ops.",9,17)
    read(pwscfout_unit,"(5X,i2)") nsym
    call findkword(pwscfout_unit,"isym")
    backspace(pwscfout_unit)
    do isym=1 , nsym
      read(pwscfout_unit,"(/20X,a45/)") sname(isym)
      read(pwscfout_unit,"(19X,3(i6,5X))") (s(1,ipol,isym),ipol=1,3)
      read(pwscfout_unit,"(19X,3(i6,5X))") (s(2,ipol,isym),ipol=1,3)
      read(pwscfout_unit,"(19X,3(i6,5X)/)") (s(3,ipol,isym),ipol=1,3)
      read(pwscfout_unit,"(19X,3F11.7)") (sr(1,ipol,isym),ipol=1,3)
      read(pwscfout_unit,"(19X,3F11.7)") (sr(2,ipol,isym),ipol=1,3)
      read(pwscfout_unit,"(19X,3F11.7/)") (sr(3,ipol,isym),ipol=1,3)
      
    enddo 
    
    
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
    
    call findkword(pwscfout_unit,"number")
    read(pwscfout_unit,"(24X,I6,/)") nirkps
    allocate(xk(3,nirkps),wk(nirkps))
    xk = 0.0
    wk = 0.0
    do ik=1,nirkps
      read(pwscfout_unit,"(20X,3F12.7,7X,F12.7)") (xk(ipol,ik),ipol=1,3),wk(ik)
    enddo
    
    allocate(et(nband,nirkps),wg(nband,nirkps))
    call findkline(pwscfout_unit,"          k =",1,13)
    backspace(pwscfout_unit)
    do ik=1,nirkps
      read(pwscfout_unit,"(/,A,/)") ctmp
      read(pwscfout_unit,"(2X,8F9.4)") (et(iband,ik),iband=1,nband)
      read(pwscfout_unit,"(A,/)") ctmp
      read(pwscfout_unit,"(2X,8F9.4)") (wg(iband,ik),iband=1,nband)
    enddo
    
    
    
    call close_file(pwscfout_name,pwscfout_unit)
    
  end subroutine readpwscf_out
  
end module readscf