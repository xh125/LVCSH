module readphout
  use io,only : io_file_unit,open_file,close_file,findkword,findkline,stdout
  use kinds,only :dp
  use constants,only : maxlen
  use phdisp,only : nq1,nq2,nq3,nqs,iq,x_q,omega_cm,omega_thz,omega_disp,phid, &
                    disdyn
  use modes, only : nirr,irr,nmodes
  use efield_mod, only : epsil,zeu
  use ions_base, only : nat,ntyp,atm,amass,ityp,tau
  use cell_base, only : ibrav,celldm,at
  implicit none
  integer :: imode0
  contains
  
  
  subroutine readph_out(phoutname)
    use parameters,only : fildyn
    use utility, only : int_to_char
    implicit none
    character(len=*),intent(in) :: phoutname
    integer :: phoutunit,dynunit
    integer :: i,imode,j,nt,nu_i,na,nb,iat
    integer :: ipol
    !counter on polarization
    integer :: ios
    character(len=maxlen) :: phrun_info
    character(len=maxlen) :: dynname
    character(len=maxlen) :: atm1,ctmp,line,ctmp1,ctmp2
    real(kind=dp) :: ftmp1,ftmp2
    
    phoutunit = io_file_unit()
    
    call open_file(phoutname,phoutunit)
    
    call findkword(phoutunit,"Program") 
    read(phoutunit,"(A)") phrun_info
    write(stdout,*) phrun_info
    
    nq1=1
    nq2=1
    nq3=1
    call findkword(phoutunit,"Dynamical")
    read(phoutunit,"(29X,2(i2,X),i2)") nq1,nq2,nq3
    
    nqs=1
    call findkline(phoutunit,"q-points):",12,21)
    read(phoutunit,"(6X,i4)") nqs
    if(.not. allocated(x_q)) allocate(x_q(3,nqs))
    x_q = 0.0
    read(phoutunit,*)  
    do iq=1,nqs
      read(phoutunit,"(8X,3f14.9)") (x_q(ipol,iq),ipol=1,3)
    enddo
    
    nmodes = 3*nat
    if(.not. allocated(omega_disp)) allocate(omega_disp(nmodes,nqs))
    if(.not. allocated(omega_thz)) allocate(omega_thz(nmodes,nqs))
    if(.not. allocated(omega_cm)) allocate(omega_cm(nmodes,nqs))
    !allocate(omega_disp(nmodes,nqs),omega_thz(nmodes,nqs),omega_cm(nmodes,nqs))
    omega_disp = 0.0
    omega_cm   = 0.0
    omega_thz  = 0.0
    
    if(.not. allocated(phid)) allocate(phid(3,3,nat,nat,nqs))
    if(.not. allocated(disdyn)) allocate(disdyn(3,nat,nmodes,nqs))
    
    do iq=1,nqs
      call findkline(phoutunit,"     Calculation of q =",1,23)
      call findkline(phoutunit,"     Computing dynamical matrix for",1,35)
      call findkline(phoutunit,"     Atomic displacements:",1,26)
      call findkline(phoutunit," irreducible representations",20,47)
      !read(phoutunit,"(15X,i4)") nirr
      !allocate(npert(nirr))
      !npert = 1
      !allocate(name_rap_mode(3*nat))
      !name_rap_mode = ' '
      !imode0=0
      !do irr=1,nirr
      !  read(phoutunit,"(/,25X,i7,8X,A)") npert(irr), name_rap_mode(imode0+1)
      !  imode0 = imode0 + npert(irr)
      !enddo 
      call findkline(phoutunit,"Diagonalizing the dynamical matrix",6,39)
      read(phoutunit,"(//,11X,3F14.9,//)") (x_q(ipol,iq),ipol=1,3)
      
      do imode=1,nmodes
        read(phoutunit,"(19X,F15.6,8X,F15.6)") omega_thz(imode,iq),omega_cm(imode,iq)
      enddo
    enddo
    
    call close_file(phoutname,phoutunit)

    
    do iq=1,nqs
    
      dynname = trim(fildyn)//trim(int_to_char(iq))
      dynunit = io_file_unit()
      call open_file(dynname,dynunit)
      
      read(dynunit,*)
      read(dynunit,*)
      read(dynunit,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
      if (ibrav==0) then
        read(dynunit,'(a)') atm1
        read(dynunit,*) ((at(i,j),i=1,3),j=1,3)
      endif
      do nt=1,ntyp
        read(dynunit,*) i,atm(nt) !,amass(nt)
      enddo
      if(.not. allocated(ityp)) allocate(ityp(nat))
      if(.not. allocated(tau)) allocate(tau(3,nat))
      
      !allocate(ityp(nat),tau(3,nat))
      do na=1,nat
        read(dynunit,*) i,ityp(na),(tau(j,na),j=1,3)
      enddo
      
      read(dynunit,"(/,A)") line
      if(line(6:14) /= "Dynamical") then
        do while(line(6:15) /= "dielectric")
          read(dynunit,"(A)",iostat=ios) line
        enddo
        read(dynunit,*) ((epsil(i,j),j=1,3),i=1,3)
        read(dynunit,*)
        read(dynunit,*)
        read(dynunit,*)
        do na=1,nat
          read(dynunit,*)
          read(dynunit,*) ((zeu(i,j,na),j=1,3),i=1,3)
        enddo
      endif
      
      read(dynunit,"(/,11X,3F14.9,/)") (x_q(i,iq),i=1,3)
      
      do na=1,nat 
        do nb=1,nat
          read(dynunit,"(2I5)") i,j
          do i=1,3
            read(dynunit,"(3(2(f12.8,1X),2X))") (phid(i,j,na,nb,iq),j=1,3)
            !read(dynunit,*) (phir(j),phii(j),j=1,3)
            !do j=1,3
            !  phiq(i,j,na,nb,iq) = CMPLX(phir(j),phii(j),kind=dp)
            !enddo
          enddo
        enddo
      enddo
      
      call findkline(dynunit,"Diagonalizing the dynamical matrix",6,39)
      read(dynunit,"(////)")
      
      do imode=1,nmodes
        read(dynunit,"(19X,F15.6,8X,F15.6)") omega_thz(imode,iq),omega_cm(imode,iq)
        !read(dynunit,"(2X,3(2f10.6),3X)") ((disdyn(ipol,iat,imode,iq),ipol=1,3),iat=1,nat)
        do iat=1,nat
          read(dynunit,"(2X,3(2f10.6))") (disdyn(ipol,iat,imode,iq),ipol=1,3)
        enddo
      enddo
      
      call close_file(dynname,dynunit)
    enddo
    

    
  end subroutine readph_out
  
end module readphout