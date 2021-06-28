module saveinf
  use kinds,only : dp
  use surfacecom,only : dt,nstep
  use io,only : io_file_unit,open_file,close_file,stdout
  use constants,only : ry_to_fs,maxlen,RYTOEV
	use parameters,only : outdir
  implicit none
  integer :: iaver,isnap,ifre,iq,imode,inode,icore
  
  character(len=maxlen) :: pes_e_file   ="pes_e.dat"
  character(len=maxlen) :: pes_h_file   ="pes_h.dat"
  character(len=maxlen) :: band_e_file  ="band_e.dat"
  character(len=maxlen) :: band_h_file  ="band_h.dat"	
  character(len=maxlen) :: eapes_e_file  ="eapes_e.dat"
  character(len=maxlen) :: eapes_h_file  ="eapes_h.dat"	
  character(len=maxlen) :: iapes_e_file  ="iapes_e.dat"
  character(len=maxlen) :: iapes_h_file  ="iapes_h.dat"		
  character(len=maxlen) :: csit_e_file  ="csit_e.dat"
  character(len=maxlen) :: csit_h_file  ="csit_h.dat"
  character(len=maxlen) :: wsit_e_file  ="wsit_e.dat"
  character(len=maxlen) :: wsit_h_file  ="wsit_h.dat" 
  character(len=maxlen) :: psit_e_file  ="psit_e.dat"
  character(len=maxlen) :: psit_h_file  ="psit_h.dat"
  character(len=maxlen) :: mskd_e_file  ="mskd_e.dat"
  character(len=maxlen) :: mskd_h_file  ="mskd_h.dat"
  character(len=maxlen) :: mskds_e_file ="mskds_e.dat"
  character(len=maxlen) :: mskds_h_file ="mskds_h.dat"	
  character(len=maxlen) :: ipr_e_file   ="ipr_e.dat"
  character(len=maxlen) :: ipr_h_file   ="ipr_h.dat"
  character(len=maxlen) :: phQ_file     = "phQsit.dat"  
  character(len=maxlen) :: phP_file     = "phPsit.dat"  
  character(len=maxlen) :: phU_file     = "phUsit.dat"  
  character(len=maxlen) :: phK_file     = "phKsit.dat"
	
  contains 

  
	!Write average phQ infortmation to the file: 	phQsit.dat.gnu	
  subroutine save_phQ(nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phQsit(nmodes,nq,0:nsnap)
		character(len=maxlen) :: phQ_file_    
    integer :: phq_unit

    phQ_file_ = trim(outdir)//trim(adjustl(phQ_file))    
    phq_unit = io_file_unit()
    call open_file(phQ_file_,phq_unit)

		
		write(phq_unit,"(A40)")  "phQsit(imode,iq,isnap)"
    do isnap=0,nsnap
      write(phq_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(phq_unit,"(7(1X,E12.5))") ((phQsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phQ_file_,phq_unit)
   
  end subroutine save_phQ
	
  subroutine read_phQ(inode,icore,nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phQsit(nmodes,nq,0:nsnap)
    
    integer :: phq_unit
		character(len=maxlen) :: phQ_file_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phQsit_(:,:,:)
		
		if(.not. allocated(phQsit_)) allocate(phQsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phQ_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phQ_file))
    
    phq_unit = io_file_unit()
    call open_file(phQ_file_,phq_unit)
		
		read(phq_unit,*)
    do isnap=0,nsnap
      read(phq_unit,"(A)") ctmp1
      read(phq_unit,"(7(1X,E12.5))") ((phQsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phQsit =phQsit + phQsit_
		
    call close_file(phQ_file_,phq_unit)
   
  end subroutine read_phQ	
	
  subroutine plot_phQ(nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phQsit(nmodes,nq,0:nsnap)
    
    integer :: phq_unit
		
		character(len=maxlen) :: phQ_file_
		phQ_file_ = trim(outdir)//trim(adjustl(phQ_file))//".gnu"
	
    phq_unit = io_file_unit()
    call open_file(phQ_file_,phq_unit)
    
		write(phq_unit,"(5X,A)") "Average of Normal mode coordinate for all trajecotry. ((phQ(imode,iq),imode=1,nmodes),iq=1,nq)"
		write(phq_unit,"(*(1X,A12))") "time ",(("phQ(ifre)",imode=1,nmodes),iq=1,nq)
		write(phq_unit,"(*(1X,A12))") " fs  ",((" a.u. ",imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
      write(phq_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,((phQsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phQ_file_,phq_unit)
		
		write(stdout,"(A,A)") "Write average phQ infortmation to the file: ",trim(phQ_file_)
   
  end subroutine plot_phQ  

	

	!Write average phP infortmation to the file: 	phPsit.dat.gnu	
  subroutine save_phP(nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phPsit(nmodes,nq,0:nsnap)
    
    integer :: php_unit
		character(len=maxlen) :: phP_file_    
    php_unit = io_file_unit()
		phP_file_ = trim(outdir)//trim(adjustl(phP_file))
    call open_file(phP_file_,php_unit)
    
		write(php_unit,"(A40)") "phPsit(imode,iq,isnap)"
    do isnap=0,nsnap
			write(php_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(php_unit,"(7(1X,E12.5))") ((phPsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phP_file_,php_unit)
   
  end subroutine save_phP  

  subroutine read_phP(inode,icore,nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phPsit(nmodes,nq,0:nsnap)
    
    integer :: php_unit
		character(len=maxlen) :: phP_file_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phPsit_(:,:,:)
		
		if(.not. allocated(phPsit_)) allocate(phPsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phP_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phP_file))
    
    php_unit = io_file_unit()
    call open_file(phP_file_,php_unit)
		
		read(php_unit,*)
    do isnap=0,nsnap
      read(php_unit,"(A)") ctmp1
      read(php_unit,"(7(1X,E12.5))") ((phPsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phPsit =phPsit + phPsit_
		
    call close_file(phP_file_,php_unit)
   
  end subroutine read_phP	
  
  subroutine plot_phP(nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phPsit(nmodes,nq,0:nsnap)
    
    integer :: php_unit

		character(len=maxlen) :: phP_file_  
		phP_file_ = trim(outdir)//trim(adjustl(phP_file))//".gnu"    
		
    php_unit = io_file_unit()
    call open_file(phP_file_,php_unit)
    
		write(php_unit,"(5X,A)") "Average of Normal mode verlocity for all trajecotry.((phP(imode,iq),imode=1,nmodes),iq=1,nq)"
		write(php_unit,"(*(1X,A12))") "time ",(("phP(ifre)",imode=1,nmodes),iq=1,nq)
		write(php_unit,"(*(1X,A12))") " fs  ",((" a.u. ",imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
      write(php_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,((phPsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phP_file_,php_unit)
  
		write(stdout,"(A,A)") "Write average phP infortmation to the file: ",trim(phP_file_)
	 
  end subroutine plot_phP  	

	
	!Write average phK infortmation to the file: phKsit.dat.gnu	
  subroutine save_phK(nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phKsit(nmodes,nq,0:nsnap)
    
    integer :: phK_unit

		character(len=maxlen) :: phK_file_    
    phK_unit = io_file_unit()
		phK_file_ = trim(outdir)//trim(adjustl(phK_file))
		
    call open_file(phK_file_,phK_unit)
    write(phK_unit,"(A40)")  "phKsit(imode,iq,isnap)"
    do isnap=0,nsnap
				write(phK_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
        write(phK_unit,"(7(1X,E12.5))") ((phKsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phK_file_,phK_unit)
   
  end subroutine save_phK

  subroutine read_phK(inode,icore,nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phKsit(nmodes,nq,0:nsnap)
    
    integer :: phk_unit
		character(len=maxlen) :: phK_file_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phKsit_(:,:,:)
		
		if(.not. allocated(phKsit_)) allocate(phKsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phK_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phK_file))
    
    phk_unit = io_file_unit()
    call open_file(phK_file_,phk_unit)
		
		read(phK_unit,*)
    do isnap=0,nsnap
      read(phk_unit,"(A)") ctmp1
      read(phk_unit,"(7(1X,E12.5))") ((phKsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phKsit =phKsit + phKsit_
		
    call close_file(phK_file_,phk_unit)
   
  end subroutine read_phK	
	
  subroutine plot_phK(nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phKsit(nmodes,nq,0:nsnap)
    
    integer :: phK_unit

		character(len=maxlen) :: phK_file_   
		phK_file_ = trim(outdir)//trim(adjustl(phK_file))//".gnu"    
    
    phK_unit = io_file_unit()
    call open_file(phK_file_,phK_unit)
		
		write(phK_unit,"(5X,A)") "Average of Normal mode kinetic energy for all trajecotry.SUM_phK,((phK(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(phK_unit,"(*(1X,A12))") "time ","SUM_phK",(("phK(mode,q)",imode=1,nmodes),iq=1,nq)
		write(phK_unit,"(*(1X,A12))") " fs  ","  eV   ",(("    eV     ",imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
        write(phK_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,SUM(phKsit(:,:,isnap))*RYTOEV,&
				((phKsit(imode,iq,isnap)*RYTOEV,imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phK_file_,phK_unit)
   
		write(stdout,"(A,A)") "Write average phK infortmation to the file: ",trim(phK_file_)
	 
  end subroutine plot_phK


	!Write average phU infortmation to the file: phUsit.dat.gnu
  subroutine save_phU(nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phUsit(nmodes,nq,0:nsnap)
    
    integer :: phU_unit
    
		character(len=maxlen) :: phU_file_
    phU_unit = io_file_unit()
		phU_file_ = trim(outdir)//trim(adjustl(phU_file))
    call open_file(phU_file_,phU_unit)

		
		write(phU_unit,"(A40)")  "phUsit(imode,iq,isnap)"
    do isnap=0,nsnap
				write(phU_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
        write(phU_unit,"(7(1X,E12.5))") ((phUsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phU_file_,phU_unit)
   
  end subroutine save_phU  

  subroutine read_phU(inode,icore,nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phUsit(nmodes,nq,0:nsnap)
    
    integer :: phu_unit
		character(len=maxlen) :: phU_file_
		character(len=maxlen) :: ctmp1,ctmp2
		real(kind=dp),allocatable :: phUsit_(:,:,:)
		
		if(.not. allocated(phUsit_)) allocate(phUsit_(nmodes,nq,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		phU_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phU_file))
    
    phu_unit = io_file_unit()
    call open_file(phU_file_,phu_unit)
		
		read(phu_unit,*)
    do isnap=0,nsnap
      read(phu_unit,"(A)") ctmp1
      read(phu_unit,"(7(1X,E12.5))") ((phUsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
		phUsit =phUsit + phUsit_
		
    call close_file(phU_file_,phu_unit)
   
  end subroutine read_phU	

  subroutine plot_phU(nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phUsit(nmodes,nq,0:nsnap)
    
    integer :: phU_unit

		character(len=maxlen) :: phU_file_  
		phU_file_ = trim(outdir)//trim(adjustl(phU_file))//".gnu"    
    
    phU_unit = io_file_unit()
    call open_file(phU_file_,phU_unit)
    
		write(phU_unit,"(5X,A)") "Average of Normal mode potential energy for all trajecotry.SUM_phU,((phU(imode,iq),imode=1,nmodes),iq=1,nq)"
		write(phU_unit,"(*(1X,A12))") "time ","SUM_phU",(("phU(mode,q)",imode=1,nmodes),iq=1,nq)
		write(phU_unit,"(*(1X,A12))") "  fs ","   eV  ",(("    eV     ",imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
        write(phU_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,SUM(phUsit(:,:,isnap))*RYTOEV,&
				((phUsit(imode,iq,isnap)*RYTOEV,imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phU_file_,phU_unit)
   
		write(stdout,"(A,A)") "Write average phU infortmation to the file: ",trim(phU_file_)
	 
  end subroutine plot_phU  



	!The electron(hole) population on different adiabatic wave function.
  subroutine save_wsit(nfre,nsnap,naver,wsit,wsit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: wsit(nfre,0:nsnap)
    character(len=*),intent(in) :: wsit_file
    character(len=maxlen) :: wsit_file_
    
		integer :: wsit_unit
    
    wsit_unit = io_file_unit()
		wsit_file_ = trim(outdir)//trim(adjustl(wsit_file))
		
    call open_file(wsit_file_,wsit_unit)
    write(wsit_unit,"(A6,I8,A40)") "naver=",naver, " wsit(ifre)=REAL(w(ifre)*CONJG(w(ifre)))"
    do isnap=0,nsnap
      write(wsit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(wsit_unit,"(7(1X,E12.5))") (wsit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(wsit_file_,wsit_unit)
  
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
		wsit_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(wsit_file))				
		
    wsit_unit = io_file_unit()
    call open_file(wsit_file_,wsit_unit)
    !write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    read(wsit_unit,*)
		do isnap=0,nsnap
      !write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(wsit_unit,*)
			read(wsit_unit,"(7(1X,E12.5))") (wsit_(ifre,isnap),ifre=1,nfre)   
    enddo
    
		wsit = wsit + wsit_
		
    call close_file(wsit_file,wsit_unit)
  
  end subroutine read_wsit	

  subroutine plot_wsit(nfre,nsnap,naver,wsit,wsit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: wsit(nfre,0:nsnap)
    character(len=*),intent(in) :: wsit_file
    
    integer :: wsit_unit
		character(len=maxlen) :: wsit_file_	
		wsit_file_ = trim(outdir)//trim(adjustl(wsit_file))//".gnu"    
    wsit_unit = io_file_unit()
    call open_file(wsit_file_,wsit_unit)

		write(wsit_unit,"(A)") "The electron(hole) population on different adiabatic wave function. "    
		write(wsit_unit,"(*(1X,A12))")          "time ",("wsit(ifre)",ifre=1,nfre)
    write(wsit_unit,"((1X,A12),*(1X,I12))") " fs  ",(ifre,ifre=1,nfre)
		do isnap=0,nsnap
      write(wsit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(wsit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(wsit_file_,wsit_unit)
		write(stdout,"(A,A)")"Write the average electron(hole) wavefction population on &
		          &different adiabatic wave function to file:",trim(wsit_file_)
  
  end subroutine plot_wsit



	! save the first trajecotry active PES and first trajecotry PES : pes_e.dat_f.gnu
	! save the average active PES and PES for all trajecotry. : pes_e.dat_average.gnu
  subroutine save_pes(nfre,nsnap,naver,pes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: pes(0:nfre,0:nsnap,1:naver)
    character(len=*),intent(in) :: pes_file
    character(len=maxlen) :: pes_file_
    
		integer :: pes_unit
		pes_file_= trim(outdir)//trim(adjustl(pes_file))

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    write(pes_unit,"(A6,I8,A)") "iaver=",1, " (pes(ifre,isnap,iaver=1),ifre=0,nfre)"
    do iaver =1 , 1
      do isnap=0,nsnap
        write(pes_unit,"(5X,A5,F11.2,A5,A26,E12.5,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs).",&
        " Energy of active surface:",pes(0,isnap,iaver)," RYD"
        write(pes_unit,"(7(1X,E12.5))") (pes(ifre,isnap,iaver),ifre=0,nfre)
      enddo
    enddo

    write(pes_unit,"(A,I8,A)") "Averager of naver=",naver, " SUM(pes(ifre,isnap,:))/naver,ifre=0,nfre)"
    do isnap=0,nsnap
      write(pes_unit,"(5X,A5,F11.2,A5,A26,E12.5,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs).",&
      " Averager Energy of PES."
      write(pes_unit,"(7(1X,E12.5))") (SUM(pes(ifre,isnap,:))/naver,ifre=0,nfre)
    enddo

 
    call close_file(pes_file_,pes_unit)
    
  end subroutine save_pes
  
  subroutine read_pes(inode,icore,nfre,nsnap,naver,pes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: pes(0:nfre,0:nsnap,1:naver)
    character(len=*),intent(in) :: pes_file
		character(len=maxlen) :: pes_file_
    
    integer :: pes_unit
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: pes_(:,:,:)
		
		if(.not. allocated(pes_)) allocate(pes_(0:nfre,0:nsnap,1:naver))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		pes_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(pes_file))		

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    read(pes_unit,*)
		do iaver =1 , 1
      do isnap=0,nsnap
        read(pes_unit,*)
				read(pes_unit,"(7(1X,E12.5))") (pes_(ifre,isnap,iaver),ifre=0,nfre)
      enddo
    enddo

		read(pes_unit,*)
    do isnap=0,nsnap
      read(pes_unit,*)
			read(pes_unit,"(7(1X,E12.5))") (pes_(ifre,isnap,2),ifre=0,nfre)
    enddo		
		
		
		pes(:,:,1) = pes_(:,:,1)
		pes(:,:,2) = pes(:,:,2) + pes_(:,:,2)
    call close_file(pes_file_,pes_unit)
    
  end subroutine read_pes
	
  subroutine plot_pes(nfre,nsnap,naver,pes,wsit,iapes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: pes(0:nfre,0:nsnap,1:naver),wsit(1:nfre,0:nsnap)
		integer,intent(in) :: iapes(0:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file
		character(len=maxlen) :: pes_file_
    
    integer :: pes_unit
		
		pes_file_ = trim(outdir)//trim(adjustl(pes_file))//"_f.gnu"

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
		write(pes_unit,"(A)") "USing the first trajecotry as an example to Plotting potential Energy Surface(PES),and the active PES."
    write(pes_unit,"(*(1X,A12))") "time ","active_fpes",("i_fpes",ifre=1,nfre)
		write(pes_unit,"(*(1X,A12))") " fs  ","   eV     ",("  eV ",ifre=1,nfre)
		do iaver =1 , 1
      do isnap=0,nsnap
				write(pes_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(pes(ifre,isnap,iaver)*RYTOEV,ifre=0,nfre)
      enddo
    enddo
		
    call close_file(pes_file_,pes_unit)
		
		write(stdout,"(A,A)") "Write the first trajecotry active PES and PES to the file:",trim(pes_file_)

		pes_file_ = trim(outdir)//trim(adjustl(pes_file))//"_average.gnu"

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
		write(pes_unit,"(A)") "Plotting averager of potential Energy Surface(PES),and the averager active PES and their weight."
    write(pes_unit,"(*(1X,A12))") "time ","average_PES","wsit      ","active"   ,"AVE_actPES"
		write(pes_unit,"(*(1X,A12))") " fs  ","   eV     " ,"w(ifre)**2" ," account "," eV "

		do ifre =1 , nfre
			write(pes_unit,"(A,I12,2(1X,A12))") "#isurface=",ifre,"wsit","weight"
			do isnap=0,nsnap
				if(ifre ==1 ) then
					write(pes_unit,"(3(1X,E12.5),1X,I12,1X,E12.5,1X,I12)") dt*nstep*isnap*ry_to_fs,pes(ifre,isnap,2)*RYTOEV,&
					wsit(ifre,isnap),iapes(ifre,isnap),pes(0,isnap,2)*RYTOEV,iapes(0,isnap)
				else
					write(pes_unit,"(3(1X,E12.5),1X,I12)") dt*nstep*isnap*ry_to_fs,pes(ifre,isnap,2)*RYTOEV,wsit(ifre,isnap),iapes(ifre,isnap)
				endif
				
			enddo
			write(pes_unit,*)
    enddo
		
    call close_file(pes_file_,pes_unit)
		
		write(stdout,"(A,A)") "Write average of active PES and PES for all trajecotry to the file:",trim(pes_file_)
    
  end subroutine plot_pes
	

	!save the number of active surface on evey PES for different time.
  subroutine save_iapes(nfre,nsnap,iapes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap
    integer,intent(in) :: iapes(0:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file
    character(len=maxlen) :: pes_file_
    
		integer :: pes_unit
		pes_file_= trim(outdir)//trim(adjustl(pes_file))

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    write(pes_unit,"(A)")  "(iapes(ifre,isnap),ifre=1,nfre)"
    do isnap=0,nsnap
      write(pes_unit,"(5X,A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(pes_unit,"(7(1X,I12))") (iapes(ifre,isnap),ifre=0,nfre)
    enddo
 
    call close_file(pes_file_,pes_unit)
    
  end subroutine save_iapes
  
  subroutine read_iapes(inode,icore,nfre,nsnap,iapes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap,inode,icore
    integer,intent(inout) :: iapes(0:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file
		character(len=maxlen) :: pes_file_
    
    integer :: pes_unit
		character(len=maxlen) :: ctmp1,ctmp2
		
		integer,allocatable :: iapes_(:,:)
		
		if(.not. allocated(iapes_)) allocate(iapes_(0:nfre,0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		pes_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(pes_file))		

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    read(pes_unit,*)

      do isnap=0,nsnap
        read(pes_unit,*)
				read(pes_unit,"(7(1X,I12))") (iapes_(ifre,isnap),ifre=0,nfre)
      enddo

		iapes = iapes +  iapes_
		
    call close_file(pes_file_,pes_unit)
    
  end subroutine read_iapes
	
  subroutine plot_iapes(nfre,nsnap,iapes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap
    integer , intent(in) :: iapes(0:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file
		character(len=maxlen) :: pes_file_
    
    integer :: pes_unit
		

		pes_file_ = trim(outdir)//trim(adjustl(pes_file))//"_.gnu"

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
		write(pes_unit,"(A)") "Plotting the number of potential Energy Surface(PES)as active PES."
    write(pes_unit,"(*(1X,A12))") "time ",("num_apes",ifre=1,nfre)
		write(pes_unit,"(*(1X,A12))") " fs  ",(" times ",ifre=1,nfre)
    do isnap=0,nsnap
				write(pes_unit,"(*(1X,I12))") dt*nstep*isnap*ry_to_fs,(iapes(ifre,isnap),ifre=0,nfre)
    enddo
		
    call close_file(pes_file_,pes_unit)
		
		write(stdout,"(A,A)") "Write the number of potential Energy Surface(PES)as active PES to the file:",trim(pes_file_)
    
  end subroutine plot_iapes
	

	
	! save the all active PES Energy and isurface for different trajecotry
	! USed to get the active PES DOS at different time
  subroutine save_apes(nsnap,naver,apes,apes_file)
    implicit none
    integer , intent(in) :: nsnap,naver
    real(kind=dp),intent(in) :: apes(0:nsnap,1:naver)
    character(len=*),intent(in) :: apes_file
    character(len=maxlen) :: apes_file_
		
    integer :: apes_unit
		apes_file_ = trim(outdir)//trim(adjustl(apes_file))
		
    apes_unit = io_file_unit()
    call open_file(apes_file_,apes_unit)
    write(apes_unit,"(A6,I8,A)") "naver=",naver, " apes(isnap,iaver)"
    do isnap=0,nsnap
      write(apes_unit,"(5X,A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(apes_unit,"(7(1X,E12.5))") (apes(isnap,iaver),iaver=1,naver)
    enddo
 
    call close_file(apes_file_,apes_unit)
    
  end subroutine save_apes

  subroutine read_apes(inode,icore,nsnap,naver,apes,apes_file)
    implicit none
    integer , intent(in) :: nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: apes(0:nsnap,1:naver)
    character(len=*),intent(in) :: apes_file
		character(len=maxlen) :: apes_file_
    
    integer :: apes_unit
		character(len=maxlen) :: ctmp1,ctmp2
		
		real(kind=dp),allocatable :: apes_(:,:)
		
		if(.not. allocated(apes_)) allocate(apes_(0:nsnap,1:naver))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		apes_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(apes_file))

    apes_unit = io_file_unit()
    call open_file(apes_file_,apes_unit)
    read(apes_unit,*)

    do isnap=0,nsnap
      read(apes_unit,*)
		  read(apes_unit,"(7(1X,E12.5))") (apes_(isnap,iaver),iaver=1,naver)
    enddo

		
		apes =  apes_
		
    call close_file(apes_file_,apes_unit)
    
  end subroutine read_apes

  subroutine plot_apes(nnode,ncore,nsnap,naver,apes,dsnap,apes_file)
    implicit none
    integer , intent(in) :: nsnap,naver,nnode,ncore,dsnap
    real(kind=dp),intent(in) :: apes(0:nsnap,1:naver*nnode*ncore)
    character(len=*),intent(in) :: apes_file
    character(len=maxlen) :: apes_file_
		
    integer :: apes_unit
		
		
		apes_file_ = trim(outdir)//trim(adjustl(apes_file))//".gnu"
    apes_unit = io_file_unit()
    call open_file(apes_file_,apes_unit)
		
		write(apes_unit,"(A)") "Plot the DOS of active potential Energy Surface."
    write(apes_unit,"(*(1X,A12))") "trajecotry ",("time(fs)",isnap=0,nsnap,dsnap)
		write(apes_unit,"((1X,A12),*(1X,F12.2))") " iaver  "   ,(dt*nstep*isnap*ry_to_fs,isnap=0,nsnap,dsnap)
    do ifre=1,naver*ncore*ncore
			write(apes_unit,"(I12,*(1X,E12.5))") ifre,(apes(isnap,ifre)*RYTOEV,isnap=0,nsnap,dsnap)
		enddo
		!do isnap=0,nsnap
    !  write(apes_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(apes(isnap,ifre)*RYTOEV,ifre=1,naver*ncore*ncore)
    !enddo
	
    call close_file(apes_file,apes_unit)
		
		write(stdout,"(A,A)") "Write the all avtive PES Energy for different trajecotry in all time to the file:",trim(apes_file_)
    
  end subroutine plot_apes



	!The electron(hole) population on different diabatic wave function.
  subroutine save_csit(nfre,nsnap,naver,csit,csit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: csit(nfre,0:nsnap)
    character(len=*),intent(in) :: csit_file
		character(len=maxlen) :: csit_file_
    
    integer :: csit_unit
    
    csit_unit = io_file_unit()
		csit_file_ = trim(outdir)//trim(adjustl(csit_file))
		
    call open_file(csit_file_,csit_unit)
    write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    do isnap=0,nsnap
      write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(csit_unit,"(7(1X,E12.5))") (csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file_,csit_unit)
  
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
		csit_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(csit_file))				
		
    csit_unit = io_file_unit()
    call open_file(csit_file_,csit_unit)
    !write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    read(csit_unit,*)
		do isnap=0,nsnap
      !write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(csit_unit,*)
			read(csit_unit,"(7(1X,E12.5))") (csit_(ifre,isnap),ifre=1,nfre)   
    enddo
    
		csit = csit + csit_
		
    call close_file(csit_file,csit_unit)
  
  end subroutine read_csit	

  subroutine plot_csit(nfre,nsnap,naver,csit,csit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: csit(nfre,0:nsnap)
    character(len=*),intent(in) :: csit_file
    
    integer :: csit_unit
		character(len=maxlen) :: csit_file_
		csit_file_ = trim(outdir)//trim(adjustl(csit_file))//".gnu"
    csit_unit = io_file_unit()
    call open_file(csit_file_,csit_unit)
    
		write(csit_unit,"(A)") "The electron(hole) population on different diabatic wave function. "
		write(csit_unit,"(*(1X,A12))") "time",("csit(ifre)",ifre=1,nfre)
    write(csit_unit,"((1X,A12),*(1X,I12))") " fs ",(ifre,ifre=1,nfre)
		do isnap=0,nsnap
      write(csit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file_,csit_unit)
		
		write(stdout,"(A,A)")"Write the average electron(hole) wavefction population on &
		          &different diabatic wave function to file:",trim(csit_file_)
		
  end subroutine plot_csit


	!The active adiabatic PES project to diabatic states. psit(ifre)=P(ifre,isurface)**2
  subroutine save_psit(nfre,nsnap,naver,psit,psit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: psit(nfre,0:nsnap)
    character(len=*),intent(in) :: psit_file
 		character(len=maxlen) :: psit_file_
		
    integer :: psit_unit
    
    psit_unit = io_file_unit()
		psit_file_ = trim(outdir)//trim(adjustl(psit_file))
		
    call open_file(psit_file_,psit_unit)
    write(psit_unit,"(A6,I8,A40)") "naver=",naver, " psit(ifre)=P(ifre,isurface)**2"
    do isnap=0,nsnap
      write(psit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(psit_unit,"(7(1X,E12.5))") (psit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(psit_file_,psit_unit)
  
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
		psit_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(psit_file))				
		
    psit_unit = io_file_unit()
    call open_file(psit_file_,psit_unit)
    !write(csit_unit,"(A6,I8,A40)") "naver=",naver, " csit(ifre)=REAL(c(ifre)*CONJG(c(ifre)))"
    read(psit_unit,*)
		do isnap=0,nsnap
      !write(csit_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      read(psit_unit,*)
			read(psit_unit,"(7(1X,E12.5))") (psit_(ifre,isnap),ifre=1,nfre)   
    enddo
    
		psit = psit + psit_
		
    call close_file(psit_file,psit_unit)
  
  end subroutine read_psit		

  subroutine plot_psit(nfre,nsnap,naver,psit,psit_file)
    implicit none
    integer,intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: psit(nfre,0:nsnap)
    character(len=*),intent(in) :: psit_file
    
    integer :: psit_unit
		character(len=maxlen) :: psit_file_		
		psit_file_ = trim(outdir)//trim(adjustl(psit_file))//".gnu"        
    psit_unit = io_file_unit()
    call open_file(psit_file_,psit_unit)
		
		write(psit_unit,"(2(1X,A12))") "time(fs) ","wsit(ifre)"
    do isnap=0,nsnap
      write(psit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(psit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(psit_file_,psit_unit)
  
	  write(stdout,"(A,A)")"Write the average electron(hole) active PES project to &
		          &different diabatic wave function to file:",trim(psit_file_)
	
  end subroutine plot_psit	


  subroutine plot_band_occupatin_withtime(nband,nk,Enk,xk,nsnap,psit,csit,dsnap,band_file)
		implicit none
		integer , intent(in) :: nband,nk,dsnap
		real(kind=dp),intent(in) :: Enk(nband,nk)
		real(kind=dp),intent(in) :: xk(3,nk)
		integer , intent(in) :: nsnap
		real(kind=dp),intent(in) :: psit(1:nband*nk,0:nsnap),csit(1:nband*nk,0:nsnap)
		character(len=*),intent(in) :: band_file
		
		character(len=maxlen) :: band_file_,ctmp1
		integer :: band_unit
		integer :: ipol,iband,ik,ifre
		
		

			!write(ctmp1,*)iband
			band_file_=trim(band_file)//".gnu"
			band_unit = io_file_unit()
			call open_file(band_file_,band_unit)
			write(band_unit,"(A,I12)") "Carrier occupation on the band structure at different time on band:",iband
			write(band_unit,"(*(1X,A12))")    "iband","ik","kx","ky","kz","E_nk",("psit","wsit",isnap=0,nsnap,dsnap)
			write(band_unit,"(6(1X,A12),*(1X,F12.2))") "index_band"," index_k","2pi/a0","2pi/a0", "2pi/a0", " eV" ,&
																						(dt*nstep*isnap*ry_to_fs,dt*nstep*isnap*ry_to_fs,isnap=0,nsnap,dsnap)
		do iband=1,nband
			write(band_unit,"(A,I12)") "#band=",iband
			do ik=1,nk
				ifre = (ik-1)*nband+iband
				write(band_unit,"(2(1X,I12),*(1X,F12.5))") iband,ik,(xk(ipol,ik),ipol=1,3),Enk(iband,ik)*RYTOEV,(psit(ifre,isnap),csit(ifre,isnap),isnap=0,nsnap,dsnap)
			enddo
			write(band_unit,*)
		enddo
		call close_file(band_file_,band_unit)
		write(stdout,"(A,A)") "Write Carrier occupation on the band structure at different time to file:",trim(band_file_)
	
	end subroutine plot_band_occupatin_withtime

 
	subroutine save_mskd(nsnap,mskd,mskd_file)
		implicit none
		integer,intent(in) :: nsnap
		real(kind=dp),intent(in) :: mskd(0:nsnap)
		character(len=*),intent(in) :: mskd_file
		
		character(len=maxlen) :: mskd_file_		
		integer :: mskd_unit
		
		mskd_unit = io_file_unit()
		mskd_file_ = trim(outdir)//trim(adjustl(mskd_file))
		call open_file(mskd_file_,mskd_unit)
		
		
		write(mskd_unit,"(A40)")  "mskd(isnap)"
		do isnap=0,nsnap
			write(mskd_unit,"(A5,F11.2,A4)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)"
      write(mskd_unit,"(7(1X,E12.5))") mskd(isnap)
    enddo  			
		
		call close_file(mskd_file_,mskd_unit)
		
	end subroutine save_mskd
 
	subroutine read_mskd(inode,icore,nsnap,mskd,mskd_file)
		implicit none
		integer,intent(in) :: nsnap,inode,icore
		real(kind=dp),intent(inout) :: mskd(0:nsnap)
		character(len=*),intent(in) :: mskd_file
		
		integer :: mskd_unit
		character(len=maxlen) :: mskd_file_,ctmp1,ctmp2
		
		real(kind=dp),allocatable :: mskd_(:)
		if(.not. allocated(mskd_)) allocate(mskd_(0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		mskd_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(mskd_file))
		
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

	subroutine plot_mskd(nsnap,mskd,mskd_file)
		implicit none
		integer,intent(in) :: nsnap
		real(kind=dp),intent(in) :: mskd(0:nsnap)
		character(len=*),intent(in) :: mskd_file
		
		integer :: mskd_unit

		character(len=maxlen) :: mskd_file_
		mskd_file_ = trim(outdir)//trim(adjustl(mskd_file))//".gnu"    
		
		mskd_unit = io_file_unit()
		call open_file(mskd_file_,mskd_unit)
		
		write(mskd_unit,"(2(1X,A12))") "time(fs) ","mskd"
		do isnap=0,nsnap
      write(mskd_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,mskd(isnap)
    enddo  			
		
		call close_file(mskd_file,mskd_unit)
		
	end subroutine plot_mskd
 

	
	subroutine save_mskds(nsnap,naver,mskds,mskds_file)
		implicit none
		integer,intent(in) :: nsnap,naver
		real(kind=dp),intent(in) :: mskds(0:nsnap,naver)
		character(len=*),intent(in) :: mskds_file
		
		character(len=maxlen) :: mskds_file_
		integer :: mskds_unit
		
		mskds_unit = io_file_unit()
		mskds_file_ = trim(outdir)//trim(adjustl(mskds_file))
		call open_file(mskds_file_,mskds_unit)
		
		
		write(mskds_unit,"(A6,I8,A40)") "naver=",naver, "mskds(isnap,iaver)"
		do isnap=0,nsnap
      write(mskds_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
			write(mskds_unit,"(7(1X,E12.5))") (mskds(isnap,iaver),iaver=1,naver)
    enddo  			
		
		call close_file(mskds_file_,mskds_unit)		
	
	end subroutine
	
	subroutine read_mskds(inode,icore,nsnap,naver,mskds,mskds_file)
		implicit none
		integer,intent(in) :: nsnap,naver,inode,icore
		real(kind=dp),intent(out) :: mskds(0:nsnap,naver)
		character(len=*),intent(in) :: mskds_file
		
		integer :: mskds_unit
		character(len=maxlen) :: mskds_file_,ctmp1,ctmp2
		
		real(kind=dp),allocatable:: mskds_(:,:)
		if(.not. allocated(mskds_)) allocate(mskds_(0:nsnap,naver))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		mskds_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(mskds_file))		
		
		
		mskds_unit = io_file_unit()
		call open_file(mskds_file_,mskds_unit)
		
		read(mskds_unit,*) 
		do isnap=0,nsnap
      read(mskds_unit,*)
			!write(mskds_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
			read(mskds_unit,"(7(1X,E12.5))") (mskds_(isnap,iaver),iaver=1,naver)
    enddo  			
		
		mskds = mskds_
		
		call close_file(mskds_file,mskds_unit)		
	
	end subroutine	

	subroutine plot_mskds(nnode,ncore,nsnap,naver,mskds,mskds_file)
		implicit none
		integer,intent(in) :: nsnap,naver,nnode,ncore
		real(kind=dp),intent(in) :: mskds(0:nsnap,naver*nnode*ncore)
		character(len=*),intent(in) :: mskds_file
		
		integer :: mskds_unit
		
		character(len=maxlen) :: mskds_file_ 
		mskds_file_ = trim(outdir)//trim(adjustl(mskds_file))//".gnu"    
		
		mskds_unit = io_file_unit()
		call open_file(mskds_file,mskds_unit)
		
		write(mskds_unit,"(2(1X,A12))") "time(fs) ","mskds"
		do isnap=0,nsnap
			write(mskds_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(mskds(isnap,iaver),iaver=1,naver*ncore*nnode)
    enddo  			
		
		call close_file(mskds_file,mskds_unit)		
	
	end subroutine plot_mskds
	

	
	subroutine save_ipr(nsnap,ipr,ipr_file)
		implicit none
		integer,intent(in) :: nsnap
		real(kind=dp),intent(in) :: ipr(0:nsnap)
		character(len=*),intent(in) :: ipr_file
		
		integer :: ipr_unit
		character(len=maxlen) :: ipr_file_
		ipr_file_ = trim(outdir)//trim(adjustl(ipr_file))
		ipr_unit = io_file_unit()
		call open_file(ipr_file_,ipr_unit)

		
		write(ipr_unit,"(A40)")  "ipr(isnap)"
		do isnap=0,nsnap
			write(ipr_unit,"(A5,F11.2,A5)") "time=",dt*nstep*isnap*ry_to_fs,"(fs)."
      write(ipr_unit,"(7(1X,E12.5))") ipr(isnap)
    enddo  			
		
		call close_file(ipr_file_,ipr_unit)	
	end subroutine save_ipr
	
	subroutine read_ipr(inode,icore,nsnap,ipr,ipr_file)
		implicit none
		integer,intent(in) :: nsnap,inode,icore
		real(kind=dp),intent(inout) :: ipr(0:nsnap)
		character(len=*),intent(in) :: ipr_file
		
		integer :: ipr_unit
		character(len=maxlen) :: ipr_file_,ctmp1,ctmp2
		
		real(kind=dp),allocatable :: ipr_(:)
		if(.not. allocated(ipr_)) allocate(ipr_(0:nsnap))
		
		write(ctmp1,*) inode
		write(ctmp2,*) icore
		ipr_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(ipr_file))
		
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

	subroutine plot_ipr(nsnap,ipr,ipr_file)
		implicit none
		integer,intent(in) :: nsnap
		real(kind=dp),intent(in) :: ipr(0:nsnap)
		character(len=*),intent(in) :: ipr_file
		
		integer :: ipr_unit

		character(len=maxlen) :: ipr_file_   
		ipr_file_ = trim(outdir)//trim(adjustl(ipr_file))//".gnu"    
		
		ipr_unit = io_file_unit()
		call open_file(ipr_file_,ipr_unit)
		
		write(ipr_unit,"(2(1X,A12))") "time(fs) ","ipr"
		do isnap=0,nsnap
      write(ipr_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,ipr(isnap)
    enddo  			
		
		call close_file(ipr_file_,ipr_unit)	
	end subroutine plot_ipr
	

	
end module saveinf