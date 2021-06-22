module saveinf
  use kinds,only : dp
  use surfacecom,only : dt,nstep
  use io,only : io_file_unit,open_file,close_file
  use constants,only : ry_to_fs,maxlen
  implicit none
  integer :: iaver,isnap,ifre,iq,imode
  
  character(len=9) :: pes_e_file  ="pes_e.dat"
  character(len=9) :: pes_h_file  ="pes_h.dat"
  character(len=10):: apes_e_file ="apes_e.dat"
  character(len=10):: apes_h_file ="apes_h.dat"	
  character(len=10):: csit_e_file ="csit_e.dat"
  character(len=10):: csit_h_file ="csit_h.dat"
  character(len=10):: wsit_e_file ="wsit_e.dat"
  character(len=10):: wsit_h_file ="wsit_h.dat" 
  character(len=10):: psit_e_file ="psit_e.dat"
  character(len=10):: psit_h_file ="psit_h.dat"
  character(len=10):: mskd_e_file ="mskd_e.dat"
  character(len=10):: mskd_h_file ="mskd_h.dat"
  character(len=11):: mskds_e_file="mskds_e.dat"
  character(len=11):: mskds_h_file="mskds_h.dat"	
  character(len=9):: ipr_e_file   ="ipr_e.dat"
  character(len=9):: ipr_h_file   ="ipr_h.dat"  
  contains 

  subroutine save_apes(nsnap,naver,apes,apes_filename)
    implicit none
    integer , intent(in) :: nsnap,naver
    real(kind=dp),intent(in) :: apes(0:nsnap,1:naver)
    character(len=*),intent(in) :: apes_filename
    
    integer :: apes_unit

    apes_unit = io_file_unit()
    call open_file(apes_filename,apes_unit)
    write(apes_unit,"(A6,I8,A)") "naver=",naver, " apes(isnap,iaver)"
    do isnap=0,nsnap
      write(apes_unit,"(5X,A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(apes_unit,"(7(1X,E12.5))") (apes(isnap,iaver),iaver=1,naver)
    enddo
 
    call close_file(apes_filename,apes_unit)
    
  end subroutine save_apes

  subroutine read_apes(inode,icore,nsnap,naver,apes,apes_filename)
    implicit none
    integer , intent(in) :: nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: apes(0:nsnap,1:naver)
    character(len=*),intent(in) :: apes_filename
    
    integer :: apes_unit
		character(len=maxlen) :: apes_filename_
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: apes_(:,:)
		
		if(.not. allocated(apes_)) allocate(apes_(0:nsnap,1:naver))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		apes_filename_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//apes_filename		

    apes_unit = io_file_unit()
    call open_file(apes_filename_,apes_unit)
    read(apes_unit,*)

    do isnap=0,nsnap
      read(apes_unit,*)
		  read(apes_unit,"(7(1X,E12.5))") (apes_(isnap,iaver),iaver=1,naver)
    enddo

		
		apes =  apes_
		
    call close_file(apes_filename_,apes_unit)
    
  end subroutine read_apes

 
  subroutine save_pes(nfre,nsnap,naver,pes,pes_filename)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: pes(0:nfre,0:nsnap,1:naver)
    character(len=*),intent(in) :: pes_filename
    
    integer :: pes_unit

    pes_unit = io_file_unit()
    call open_file(pes_filename,pes_unit)
    write(pes_unit,"(A6,I8,A)") "naver=",naver, " (pes(ifre,isnap,iaver=1),ifre=0,nfre)"
    do iaver =1 , 1
      do isnap=0,nsnap
        write(pes_unit,"(5X,A5,F11.2,A5,A26,E12.5,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs).",&
        " Energy of active surface:",pes(0,isnap,iaver)," RYD"
        write(pes_unit,"(7(1X,E12.5))") (pes(ifre,isnap,iaver),ifre=0,nfre)
      enddo
    enddo
 
    call close_file(pes_filename,pes_unit)
    
  end subroutine save_pes
  
  subroutine read_pes(inode,icore,nfre,nsnap,naver,pes,pes_filename)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: pes(0:nfre,0:nsnap,1:naver)
    character(len=*),intent(in) :: pes_filename
    
    integer :: pes_unit
		character(len=maxlen) :: pes_filename_
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: pes_(:,:,:)
		
		if(.not. allocated(pes_)) allocate(pes_(0:nfre,0:nsnap,1:naver))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		pes_filename_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//pes_filename		

    pes_unit = io_file_unit()
    call open_file(pes_filename_,pes_unit)
    read(pes_unit,*)
		do iaver =1 , 1
      do isnap=0,nsnap
        read(pes_unit,*)
				read(pes_unit,"(7(1X,E12.5))") (pes_(ifre,isnap,iaver),ifre=0,nfre)
      enddo
    enddo
		
		pes =  pes_
		
    call close_file(pes_filename_,pes_unit)
    
  end subroutine read_pes
	
  subroutine save_csit(nfre,nsnap,naver,csit,csit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: csit(nfre,0:nsnap)
    character(len=*),intent(in) :: csit_file
    
    integer :: csit_unit
    
    csit_unit = io_file_unit()
    call open_file(csit_file,csit_unit)
    write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    do isnap=0,nsnap
      write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(csit_unit,"(7(1X,E12.5))") (csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file,csit_unit)
  
  end subroutine save_csit
	
  subroutine read_csit(inode,icore,nfre,nsnap,naver,csit,csit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: csit(nfre,0:nsnap)
    character(len=*),intent(in) :: csit_file
    
    integer :: csit_unit
		character(len=maxlen) :: csit_file_
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: csit_(:,:)
		
		if(.not. allocated(csit_)) allocate(csit_(nfre,0:nsnap))

		write(ctmp1,*) inode
		write(ctmp2,*) icore
		csit_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//csit_file				
		
    csit_unit = io_file_unit()
    call open_file(csit_file_,csit_unit)
    !write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    read(csit_unit,*)
		do isnap=0,nsnap
      !write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(csit_unit,*)
			read(csit_unit,"(7(1X,E12.5))") (csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file,csit_unit)
  
  end subroutine read_csit	

  subroutine save_wsit(nfre,nsnap,naver,wsit,wsit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: wsit(nfre,0:nsnap)
    character(len=*),intent(in) :: wsit_file
    
    integer :: wsit_unit
    
    wsit_unit = io_file_unit()
    call open_file(wsit_file,wsit_unit)
    write(wsit_unit,"(A6,I8,A40)") "naver=",naver, " wsit(ifre)=REAL(w(ifre)*CONJG(w(ifre)))"
    do isnap=0,nsnap
      write(wsit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(wsit_unit,"(7(1X,E12.5))") (wsit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(wsit_file,wsit_unit)
  
  end subroutine save_wsit

  subroutine read_wsit(inode,icore,nfre,nsnap,naver,wsit,wsit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: wsit(nfre,0:nsnap)
    character(len=*),intent(in) :: wsit_file
    
    integer :: wsit_unit
		character(len=maxlen) :: wsit_file_
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: wsit_(:,:)
		
		if(.not. allocated(wsit_)) allocate(wsit_(nfre,0:nsnap))

		write(ctmp1,*) inode
		write(ctmp2,*) icore
		wsit_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//wsit_file				
		
    wsit_unit = io_file_unit()
    call open_file(wsit_file_,wsit_unit)
    !write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    read(wsit_unit,*)
		do isnap=0,nsnap
      !write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(wsit_unit,*)
			read(wsit_unit,"(7(1X,E12.5))") (wsit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(wsit_file,wsit_unit)
  
  end subroutine read_wsit	


  subroutine save_psit(nfre,nsnap,naver,psit,psit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: psit(nfre,0:nsnap)
    character(len=*),intent(in) :: psit_file
    
    integer :: psit_unit
    
    psit_unit = io_file_unit()
    call open_file(psit_file,psit_unit)
    write(psit_unit,"(A6,I8,A40)") "naver=",naver, " psit(ifre)=P(ifre,isurface)**2"
    do isnap=0,nsnap
      write(psit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(psit_unit,"(7(1X,E12.5))") (psit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(psit_file,psit_unit)
  
  end subroutine save_psit
	
  subroutine read_psit(inode,icore,nfre,nsnap,naver,psit,psit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: psit(nfre,0:nsnap)
    character(len=*),intent(in) :: psit_file
    
    integer :: psit_unit
		character(len=maxlen) :: psit_file_
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: psit_(:,:)
		
		if(.not. allocated(psit_)) allocate(psit_(nfre,0:nsnap))

		write(ctmp1,*) inode
		write(ctmp2,*) icore
		psit_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//psit_file				
		
    psit_unit = io_file_unit()
    call open_file(psit_file_,psit_unit)
    !write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    read(psit_unit,*)
		do isnap=0,nsnap
      !write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(psit_unit,*)
			read(psit_unit,"(7(1X,E12.5))") (psit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(psit_file,psit_unit)
  
  end subroutine read_psit		
	
  
  subroutine save_phQ(nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phQsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phQ_filename = "phQsit.dat"
    integer :: phq_unit
    
    phq_unit = io_file_unit()
    call open_file(phQ_filename,phq_unit)
    
		write(phq_unit,"(A40)")  "phQsit(imode,iq,isnap)"
    do isnap=0,nsnap
      write(phq_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(phq_unit,"(7(1X,E12.5))") ((phQsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phQ_filename,phq_unit)
   
  end subroutine save_phQ
	
  subroutine read_phQ(inode,icore,nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phQsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phQ_filename = "phQsit.dat"
    integer :: phq_unit
		character(len=maxlen) :: phQ_filename_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phQsit_(:,:,:)
		
		if(.not. allocated(phQsit_)) allocate(phQsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phQ_filename_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//phQ_filename
    
    phq_unit = io_file_unit()
    call open_file(phQ_filename_,phq_unit)
		
    do isnap=0,nsnap
      read(phq_unit,"(A)") ctmp1
      read(phq_unit,"(7(1X,E12.5))") ((phQsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phQsit =phQsit + phQsit_
		
    call close_file(phQ_filename_,phq_unit)
   
  end subroutine read_phQ	
		
  
  subroutine save_phP(nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phPsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phP_filename = "phPsit.dat"
    integer :: php_unit
    
    php_unit = io_file_unit()
    call open_file(phP_filename,php_unit)
    
		write(php_unit,"(A40)") "phPsit(imode,iq,isnap)"
    do isnap=0,nsnap
			write(php_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(php_unit,"(7(1X,E12.5))") ((phPsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phP_filename,php_unit)
   
  end subroutine save_phP  

  subroutine read_phP(inode,icore,nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phPsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phP_filename = "phPsit.dat"
    integer :: php_unit
		character(len=maxlen) :: phP_filename_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phPsit_(:,:,:)
		
		if(.not. allocated(phPsit_)) allocate(phPsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phP_filename_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//phP_filename
    
    php_unit = io_file_unit()
    call open_file(phP_filename_,php_unit)
		
		read(php_unit,*)
    do isnap=0,nsnap
      read(php_unit,"(A)") ctmp1
      read(php_unit,"(7(1X,E12.5))") ((phPsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phPsit =phPsit + phPsit_
		
    call close_file(phP_filename_,php_unit)
   
  end subroutine read_phP	
  
  subroutine save_phK(nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phKsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phK_filename = "phKsit.dat"
    integer :: phK_unit
    
    phK_unit = io_file_unit()
    call open_file(phK_filename,phK_unit)
    write(phK_unit,"(A40)")  "phKsit(imode,iq,isnap)"
    do isnap=0,nsnap
				write(phK_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
        write(phK_unit,"(7(1X,E12.5))") ((phKsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phK_filename,phK_unit)
   
  end subroutine save_phK

  subroutine read_phK(inode,icore,nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phKsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phK_filename = "phKsit.dat"
    integer :: phk_unit
		character(len=maxlen) :: phK_filename_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phKsit_(:,:,:)
		
		if(.not. allocated(phKsit_)) allocate(phKsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phK_filename_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//phK_filename
    
    phk_unit = io_file_unit()
    call open_file(phK_filename_,phk_unit)
		
    do isnap=0,nsnap
      read(phk_unit,"(A)") ctmp1
      read(phk_unit,"(7(1X,E12.5))") ((phKsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phKsit =phKsit + phKsit_
		
    call close_file(phK_filename_,phk_unit)
   
  end subroutine read_phK	

  subroutine save_phU(nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phUsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phU_filename = "phUsit.dat"
    integer :: phU_unit
    
    phU_unit = io_file_unit()
    call open_file(phU_filename,phU_unit)
    
		write(phU_unit,"(A40)")  "phUsit(imode,iq,isnap)"
    do isnap=0,nsnap
				write(phU_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
        write(phU_unit,"(7(1X,E12.5))") ((phUsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phU_filename,phU_unit)
   
  end subroutine save_phU  

  subroutine read_phU(inode,icore,nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phUsit(nmodes,nq,0:nsnap)
    
    character(len=10) :: phU_filename = "phUsit.dat"
    integer :: phu_unit
		character(len=maxlen) :: phU_filename_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phUsit_(:,:,:)
		
		if(.not. allocated(phUsit_)) allocate(phUsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phU_filename_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//phU_filename
    
    phu_unit = io_file_unit()
    call open_file(phU_filename_,phu_unit)
		
		read(phu_unit,*)
    do isnap=0,nsnap
      read(phu_unit,"(A)") ctmp1
      read(phu_unit,"(7(1X,E12.5))") ((phUsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phUsit =phUsit + phUsit_
		
    call close_file(phU_filename_,phu_unit)
   
  end subroutine read_phU	
  
	subroutine save_mskd(nsnap,mskd,mskd_file)
		implicit none
		integer,intent(in) :: nsnap
		real(kind=dp),intent(in) :: mskd(0:nsnap)
		character(len=10),intent(in) :: mskd_file
		
		integer :: mskd_unit
		
		mskd_unit = io_file_unit()
		call open_file(mskd_file,mskd_unit)
		
		write(mskd_unit,"(A40)")  "mskd(isnap)"
		do isnap=0,nsnap
			write(mskd_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(mskd_unit,"(7(1X,E12.5))") mskd(isnap)
    enddo  			
		
		call close_file(mskd_file,mskd_unit)
		
	end subroutine save_mskd

	subroutine read_mskd(inode,icore,nsnap,mskd,mskd_file)
		implicit none
		integer,intent(in) :: nsnap,inode,icore
		real(kind=dp),intent(inout) :: mskd(0:nsnap)
		character(len=10),intent(in) :: mskd_file
		
		integer :: mskd_unit
		character(len=maxlen) :: mskd_file_,ctmp1,ctmp2
		
		real(kind=dp),allocatable :: mskd_(:)
		if(.not. allocated(mskd_)) allocate(mskd_(0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		mskd_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//mskd_file
		
		mskd_unit = io_file_unit()
		call open_file(mskd_file_,mskd_unit)
		
		read(mskd_unit,*)
		!write(mskd_unit,"(A6,I8,A40)") "naver=",naver, "mskd(isnap)"
		do isnap=0,nsnap
			read(mskd_unit,*)
			!write(mskd_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(mskd_unit,"(7(1X,E12.5))") mskd(isnap)
    enddo  			
		
		mskd =mskd + mskd_
		
		call close_file(mskd_file_,mskd_unit)
		
	end subroutine read_mskd
	
	subroutine save_mskds(nsnap,naver,mskds,mskds_file)
		implicit none
		integer,intent(in) :: nsnap,naver
		real(kind=dp),intent(in) :: mskds(0:nsnap,naver)
		character(len=11),intent(in) :: mskds_file
		
		integer :: mskds_unit
		
		mskds_unit = io_file_unit()
		call open_file(mskds_file,mskds_unit)
		
		write(mskds_unit,"(A6,I8,A40)") "naver=",naver, "mskds(isnap,iaver)"
		do isnap=0,nsnap
      write(mskds_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
			write(mskds_unit,"(7(1X,E12.5))") (mskds(isnap,iaver),iaver=1,naver)
    enddo  			
		
		call close_file(mskds_file,mskds_unit)		
	
	end subroutine
	
	subroutine read_mskds(inode,icore,nsnap,naver,mskds,mskds_file)
		implicit none
		integer,intent(in) :: nsnap,naver,inode,icore
		real(kind=dp),intent(inout) :: mskds(0:nsnap,naver)
		character(len=11),intent(in) :: mskds_file
		
		integer :: mskds_unit
		character(len=maxlen) :: mskds_file_,ctmp1,ctmp2
		
		real(kind=dp),allocatable:: mskds_(:,:)
		if(.not. allocated(mskds_)) allocate(mskds_(0:nsnap,naver))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		mskds_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//mskds_file		
		
		
		mskds_unit = io_file_unit()
		call open_file(mskds_file_,mskds_unit)
		
		read(mskds_unit,*) 
		do isnap=0,nsnap
      read(mskds_unit,*)
			!write(mskds_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
			read(mskds_unit,"(7(1X,E12.5))") (mskds_(isnap,iaver),iaver=1,naver)
    enddo  			
		
		mskds = mskds + mskds_
		
		call close_file(mskds_file,mskds_unit)		
	
	end subroutine	
	
	subroutine save_ipr(nsnap,ipr,ipr_file)
		implicit none
		integer,intent(in) :: nsnap
		real(kind=dp),intent(in) :: ipr(0:nsnap)
		character(len=9),intent(in) :: ipr_file
		
		integer :: ipr_unit
		
		ipr_unit = io_file_unit()
		call open_file(ipr_file,ipr_unit)
		
		write(ipr_unit,"(A40)")  "ipr(isnap)"
		do isnap=0,nsnap
			write(ipr_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(ipr_unit,"(7(1X,E12.5))") ipr(isnap)
    enddo  			
		
		call close_file(ipr_file,ipr_unit)	
	end subroutine save_ipr
	
	subroutine read_ipr(inode,icore,nsnap,ipr,ipr_file)
		implicit none
		integer,intent(in) :: nsnap,inode,icore
		real(kind=dp),intent(inout) :: ipr(0:nsnap)
		character(len=9),intent(in) :: ipr_file
		
		integer :: ipr_unit
		character(len=maxlen) :: ipr_file_,ctmp1,ctmp2
		
		real(kind=dp),allocatable :: ipr_(:)
		if(.not. allocated(ipr_)) allocate(ipr_(0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		ipr_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//ipr_file
		
		ipr_unit = io_file_unit()
		call open_file(ipr_file_,ipr_unit)
		
		read(ipr_unit,*)
		!write(ipr_unit,"(A6,I8,A40)") "naver=",naver, "ipr(isnap)"
		do isnap=0,nsnap
			read(ipr_unit,*)
			!write(ipr_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(ipr_unit,"(7(1X,E12.5))") ipr(isnap)
    enddo  			
		
		ipr =ipr + ipr_
		
		call close_file(ipr_file_,ipr_unit)
		
	end subroutine read_ipr	
	
end module saveinf