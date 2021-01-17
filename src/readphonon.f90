module readphonon
  use kinds,only:dp
  use parameters,only:nq1,nq2,nq3,ntotq,nmode,nfreem,ifreem
  use constants
  implicit none
    integer         :: natom,iatom
    real(kind=dp),allocatable ::  womiga(:)
    real(kind=dp),allocatable :: phQ(:),phP(:),phQ0(:),phP0(:),&
                                 phQ_max(:),phQ_sita(:)
    
  contains
  
  !!get the phonon energy and convert to Hatree unit
  subroutine get_phononen(dir)
    use readposcar
    implicit none
    character(len=maxlen),intent(in) :: dir
    ntotq = nq1*nq2*nq3
    call getPOSCAR(dir)
    natom = num_atoms
    nfreem = ntotq*natom*3-3
    allocate(womiga(nfreem))
    allocate(phQ(nfreem),phP(nfreem),phQ0(nfreem),phP0(nfreem),phQ_max(nfreem))
    
    if(ntotq == 1) then!phonon use gamma only
      call readwomiga(dir)
      womiga = womiga/Au2mev
    else
    
    endif

    
  end subroutine get_phononen
  

  
  !====================================================================!
  !  read phonon womiga(mev) !use the OUTCAR of phonon calculate       !
  !====================================================================!
  subroutine readwomiga(dir)
    use constants
    use io
    
    implicit none
    character(len=maxlen),intent(in) :: dir
    integer              :: wvecter_unit
    character(len=maxlen):: wvecter_name 
    character(len=maxlen):: ctmp
    logical              :: lexist
    
    wvecter_unit = io_file_unit()
    wvecter_name = trim(adjustl(dir))//"phonon-gamma/wvecter.txt"
    inquire(file=wvecter_name,exist=lexist)
    !inquire(file="./phonon-gamma/wvecter.txt",exist=lexist)
    if(.Not. lexist) then
      write(stdout,*) "Wrong: the file:",trim(adjustl(dir))//"phonon-gamma/wvecter.txt is not found!"
      stop
    endif
    
    call open_file(wvecter_name,wvecter_unit)
    read(wvecter_unit,'(//,A)') ctmp
    do ifreem=1,nfreem
      read(wvecter_unit,'(A)') ctmp
      read(wvecter_unit,'(T63,F12.7)') womiga(ifreem) !mev
      do iatom=1,natom+1
        read(wvecter_unit,*) ctmp
      enddo
    enddo
    
    call close_file(wvecter_name,wvecter_unit)
    write(stdout,"(/,1X,A67)") repeat("=",67)
    write(stdout,"(T5,A)")     "Womiga information:"
    write(stdout,"(1X,A67)")   repeat("=",67)    
    do ifreem=1,nfreem
      write(stdout,*) "womiga(",ifreem,")=",womiga(ifreem),"mev"
    enddo
    write(stdout,"(1X,A67)")   repeat("=",67) 
    
  end subroutine readwomiga
  

  
  
end module readphonon