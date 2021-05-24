module saveinf
  use kinds,only : dp
  use surfacecom,only : dt,nstep
  use io,only : io_file_unit,open_file,close_file
  use constants,only : ry_to_fs
  implicit none
  integer :: iaver,isnap,ifre,iq,imode
  
  character(len=9) :: pes_e_file ="pes_e.dat"
  character(len=9) :: pes_h_file ="pes_h.dat"
  character(len=10):: csit_e_file="csit_e.dat"
  character(len=10):: csit_h_file="csit_h.dat"
  character(len=10):: wsit_e_file="wsit_e.dat"
  character(len=10):: wsit_h_file="wsit_h.dat" 
  character(len=10):: psit_e_file="psit_e.dat"
  character(len=10):: psit_h_file="psit_h.dat"   
  
  contains 
  
  subroutine save_pes(nfre,nsnap,naver,pes,pes_filename)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: pes(0:nfre,1:nsnap,1:naver)
    character(len=*),intent(in) :: pes_filename
    
    integer :: pes_unit

    pes_unit = io_file_unit()
    call open_file(pes_filename,pes_unit)
    !write(pes_unit,"(A)") "dt(fs) active_energy pes(RYD)"
    do iaver =1 , 1
      write(pes_unit,"(A6,I8)") "iaver=",iaver
      do isnap=1,nsnap
        write(pes_unit,"(/5X,A5,F11.2,A5,A26,E12.5,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs).",&
        " Energy of active surface:",pes(0,isnap,iaver)," RYD"
        write(pes_unit,"(7(1X,E12.5))") (pes(ifre,isnap,iaver),ifre=1,nfre)
      enddo
    enddo
 
    call close_file(pes_filename,pes_unit)
    
  end subroutine save_pes
  
  subroutine save_csit(nfre,nsnap,naver,csit,csit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: csit(nfre,nsnap)
    character(len=*),intent(in) :: csit_file
    
    integer :: csit_unit
    
    csit_unit = io_file_unit()
    call open_file(csit_file,csit_unit)
    write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    do isnap=1,nsnap
      write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(csit_unit,"(7(1X,E12.5))") (csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file,csit_unit)
  
  end subroutine save_csit

  subroutine save_wsit(nfre,nsnap,naver,wsit,wsit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: wsit(nfre,nsnap)
    character(len=*),intent(in) :: wsit_file
    
    integer :: wsit_unit
    
    wsit_unit = io_file_unit()
    call open_file(wsit_file,wsit_unit)
    write(wsit_unit,"(A6,I8,A40)") "naver=",naver, " wsit(ifre)=REAL(w(ifre)*CONJG(w(ifre)))"
    do isnap=1,nsnap
      write(wsit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(wsit_unit,"(7(1X,E12.5))") (wsit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(wsit_file,wsit_unit)
  
  end subroutine save_wsit

  subroutine save_psit(nfre,nsnap,naver,psit,psit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: psit(nfre,nsnap)
    character(len=*),intent(in) :: psit_file
    
    integer :: psit_unit
    
    psit_unit = io_file_unit()
    call open_file(psit_file,psit_unit)
    write(psit_unit,"(A6,I8,A40)") "naver=",naver, " psit(ifre)=P(ifre,isurface)**2"
    do isnap=1,nsnap
      write(psit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(psit_unit,"(7(1X,E12.5))") (psit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(psit_file,psit_unit)
  
  end subroutine save_psit
  
  subroutine save_phQ(nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phQsit(nmodes,nq,nsnap)
    
    character(len=10) :: phQ_filename = "phQsit.dat"
    integer :: phq_unit
    
    phq_unit = io_file_unit()
    call open_file(phQ_filename,phq_unit)
    
    do isnap=1,nsnap
      write(phq_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(phq_unit,"(7(1X,E12.5))") ((phQsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phQ_filename,phq_unit)
   
  end subroutine save_phQ
  
  subroutine save_phP(nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phPsit(nmodes,nq,nsnap)
    
    character(len=10) :: phP_filename = "phPsit.dat"
    integer :: php_unit
    
    php_unit = io_file_unit()
    call open_file(phP_filename,php_unit)
    
    do isnap=1,nsnap
      do iq=1,nq
        write(php_unit,"(F11.2,7(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(phPsit(imode,iq,isnap),imode=1,nmodes)
      enddo
    enddo
    
    call close_file(phP_filename,php_unit)
   
  end subroutine save_phP  
  
  subroutine save_phK(nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phKsit(nmodes,nq,nsnap)
    
    character(len=10) :: phK_filename = "phKsit.dat"
    integer :: phK_unit
    
    phK_unit = io_file_unit()
    call open_file(phK_filename,phK_unit)
    
    do isnap=1,nsnap
      do iq=1,nq
        write(phK_unit,"(F11.2,7(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(phKsit(imode,iq,isnap),imode=1,nmodes)
      enddo
    enddo
    
    call close_file(phK_filename,phK_unit)
   
  end subroutine save_phK

  subroutine save_phU(nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phUsit(nmodes,nq,nsnap)
    
    character(len=10) :: phU_filename = "phUsit.dat"
    integer :: phU_unit
    
    phU_unit = io_file_unit()
    call open_file(phU_filename,phU_unit)
    
    do isnap=1,nsnap
      do iq=1,nq
        write(phU_unit,"(F11.2,7(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(phUsit(imode,iq,isnap),imode=1,nmodes)
      enddo
    enddo
    
    call close_file(phU_filename,phU_unit)
   
  end subroutine save_phU  
  
end module saveinf