module saveinf
  use kinds,only : dp,dpc
  use surfacecom,only : dt,nstep,temp,nqv,l_ph_quantum,eps_acustic,indexk,indexq
  use io,only : io_file_unit,open_file,close_file,stdout
  use constants,only : ry_to_fs,maxlen,RYTOEV,ryd2meV,K_B_Ryd
  use elph2,only : wf,xqf,wf_sub
  use parameters,only : outdir
  implicit none
  integer :: iaver,isnap,ifre,iq,imode,inode,icore
  integer :: i
  
  character(len=maxlen) :: pes_e_file  = "pes_e.dat"
  character(len=maxlen) :: pes_h_file  = "pes_h.dat"
  character(len=maxlen) :: csit_e_file = "csit_e.dat"
  character(len=maxlen) :: csit_h_file = "csit_h.dat"
  character(len=maxlen) :: wsit_e_file = "wsit_e.dat"
  character(len=maxlen) :: wsit_h_file = "wsit_h.dat" 
  character(len=maxlen) :: psit_e_file = "psit_e.dat"
  character(len=maxlen) :: psit_h_file = "psit_h.dat"
  character(len=maxlen) :: phQ_file    = "phQsit.dat"  
  character(len=maxlen) :: phP_file    = "phPsit.dat"  
  character(len=maxlen) :: phU_file    = "phUsit.dat"  
  character(len=maxlen) :: phK_file    = "phKsit.dat"
  character(len=maxlen) :: phT_file    = "phLOtemp.dat"
  character(len=maxlen) :: band_e_file = "band_e"
  character(len=maxlen) :: band_h_file = "band_h"  
  character(len=maxlen) :: e_temp_file = "temp_e"
  character(len=maxlen) :: h_temp_file = "temp_h"
  
  character(len=12) :: ctmp1,ctmp2
  character(len=12),allocatable :: cphmode(:,:)
  character(len=12),allocatable :: cefree(:),cnk(:,:)
  
  contains 

  
  !Write average phQ infortmation to the file:   phQsit.dat
  subroutine save_phQ(nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    complex(kind=dpc),intent(in) :: phQsit(nmodes,nq,0:nsnap)
 
    integer :: phq_unit
    character(len=maxlen) :: phQ_file_   
    phQ_file_ = trim(outdir)//trim(adjustl(phQ_file))    
    phq_unit = io_file_unit()
    call open_file(phQ_file_,phq_unit)

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo

    write(phq_unit,"(A)") "Average of Normal mode coordinate for one core sample. ((phQ(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(phq_unit,"(A12,*(2(1X,A12)))") "time",(("Re[Q"//trim(adjustl(cphmode(imode,iq)))//"]",&
    "Im[Q"//trim(adjustl(cphmode(imode,iq)))//"]",imode=1,nmodes),iq=1,nq)
    write(phq_unit,"(A12,*(2(1X,A12)))") "fs  ",((" a.u. "," a.u. ",imode=1,nmodes),iq=1,nq)
    write(phq_unit,"(A12,*(2(1X,F12.5)))") "Omega(meV)",(((wf_sub(imode,iq)*ryd2meV,wf_sub(imode,iq)*ryd2meV),imode=1,nmodes),iq=1,nq)

    do isnap=0,nsnap
      write(phq_unit,"(F12.5,*(2(1X,E12.5)))") dt*nstep*isnap*ry_to_fs,((phQsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phQ_file_,phq_unit)
    
    write(stdout,"(/A,A)") "Save average phQ infortmation to the file: ",trim(phQ_file_)
  
  end subroutine save_phQ
  
  subroutine read_phQ(inode,icore,nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    complex(kind=dpc),intent(inout) :: phQsit(nmodes,nq,0:nsnap)
    
    integer :: phq_unit
    character(len=maxlen) :: phQ_file_
    character(len=maxlen) :: ctmp1,ctmp2
    complex(kind=dpc),allocatable :: phQsit_(:,:,:)
    
    if(.not. allocated(phQsit_)) allocate(phQsit_(nmodes,nq,0:nsnap))
    
    write(ctmp1,*) inode
    write(ctmp2,*) icore
    phQ_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phQ_file))
    
    phq_unit = io_file_unit()
    call open_file(phQ_file_,phq_unit)
    
    do i=1,4
      read(phq_unit,*)
    enddo
    do isnap=0,nsnap
      read(phq_unit,"(13X,*(2(1X,E12.5)))") ((phQsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    phQsit =phQsit + phQsit_
    
    call close_file(phQ_file_,phq_unit)
   
  end subroutine read_phQ  
  
  subroutine plot_phQ(nmodes,nq,nsnap,phQsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    complex(kind=dpc),intent(in) :: phQsit(nmodes,nq,0:nsnap)
    
    integer :: phq_unit
    character(len=maxlen) :: phQ_file_
    phQ_file_ = trim(outdir)//trim(adjustl(phQ_file)) 
    phq_unit = io_file_unit()
    call open_file(phQ_file_,phq_unit)

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo
    
    write(phq_unit,"(A)") "Average of Normal mode coordinate for one core sample. ((phQ(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(phq_unit,"(A12,*(2(1X,A12)))") "time",(("Re[Q"//trim(adjustl(cphmode(imode,iq)))//"]",&
    "Im[Q"//trim(adjustl(cphmode(imode,iq)))//"]",imode=1,nmodes),iq=1,nq)
    write(phq_unit,"(A12,*(2(1X,A12)))") "fs  ",((" a.u. "," a.u. ",imode=1,nmodes),iq=1,nq)
    write(phq_unit,"(A12,*(2(1X,F12.5)))") "Omega(meV)",(((wf_sub(imode,iq)*ryd2meV,wf_sub(imode,iq)*ryd2meV),imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
      write(phq_unit,"(F12.5,*(2(1X,E12.5)))") dt*nstep*isnap*ry_to_fs,((phQsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phQ_file_,phq_unit)
    
    write(stdout,"(A,A)") "Write average phQ infortmation to the file: ",trim(phQ_file_)
   
  end subroutine plot_phQ  

  

  !Write average phP infortmation to the file:   phPsit.dat.gnu  
  subroutine save_phP(nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    complex(kind=dp),intent(in) :: phPsit(nmodes,nq,0:nsnap)
    
    integer :: php_unit
    character(len=maxlen) :: phP_file_    
    php_unit = io_file_unit()
    phP_file_ = trim(outdir)//trim(adjustl(phP_file))
    call open_file(phP_file_,php_unit)

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo

    
    write(php_unit,"(A)") "Average of Normal mode verlocity for one core sample.((phP(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(php_unit,"(A12,*(2(1X,A12)))") "time",(("Re[P"//trim(adjustl(cphmode(imode,iq)))//"]",&
    "Im[P"//trim(adjustl(cphmode(imode,iq)))//"]",imode=1,nmodes),iq=1,nq)
    !write(php_unit,"(A12,*(2(1X,A12)))") "time",(("Re[phP]","Im[phP]",imode=1,nmodes),iq=1,nq)
    write(php_unit,"(A12,*(2(1X,A12)))") "fs  ",((" a.u. "," a.u. ",imode=1,nmodes),iq=1,nq)
    write(php_unit,"(A12,*(2(1X,F12.5)))") "Omega(meV)",(((wf_sub(imode,iq)*ryd2meV,wf_sub(imode,iq)*ryd2meV),imode=1,nmodes),iq=1,nq)    
    do isnap=0,nsnap
      write(php_unit,"(F12.5,*(2(1X,E12.5)))") dt*nstep*isnap*ry_to_fs,((phPsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phP_file_,php_unit)
   
    write(stdout,"(A,A)") "Save average phP infortmation to the file: ",trim(phP_file_)
    
  end subroutine save_phP  

  subroutine read_phP(inode,icore,nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    complex(kind=dpc),intent(inout) :: phPsit(nmodes,nq,0:nsnap)
    
    integer :: php_unit
    character(len=maxlen) :: phP_file_
    character(len=maxlen) :: ctmp1,ctmp2
    complex(kind=dpc),allocatable :: phPsit_(:,:,:)
    
    if(.not. allocated(phPsit_)) allocate(phPsit_(nmodes,nq,0:nsnap))
    
    write(ctmp1,*) inode
    write(ctmp2,*) icore
    phP_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phP_file))
    
    php_unit = io_file_unit()
    call open_file(phP_file_,php_unit)
    
    do i=1,4
      read(php_unit,*)
    enddo
    do isnap=0,nsnap
      read(php_unit,"(13X,*(2(1X,E12.5)))") ((phPsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    phPsit =phPsit + phPsit_
    
    call close_file(phP_file_,php_unit)
   
  end subroutine read_phP  
  
  subroutine plot_phP(nmodes,nq,nsnap,phPsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    complex(kind=dpc),intent(in) :: phPsit(nmodes,nq,0:nsnap)
    
    integer :: php_unit
    character(len=maxlen) :: phP_file_  
    phP_file_ = trim(outdir)//trim(adjustl(phP_file))     
    php_unit = io_file_unit()
    call open_file(phP_file_,php_unit)

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo
    
    write(php_unit,"(A)") "Average of Normal mode verlocity for one core sample.((phP(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(php_unit,"(A12,*(2(1X,A12)))") "time",(("Re[P"//trim(adjustl(cphmode(imode,iq)))//"]",&
    "Im[P"//trim(adjustl(cphmode(imode,iq)))//"]",imode=1,nmodes),iq=1,nq)    
    !write(php_unit,"(A12,*(2(1X,A12)))") "time",(("Re[phP]","Im[phP]",imode=1,nmodes),iq=1,nq)
    write(php_unit,"(A12,*(2(1X,A12)))") "fs  ",((" a.u. "," a.u. ",imode=1,nmodes),iq=1,nq)
    write(php_unit,"(A12,*(2(1X,F12.5)))") "Omega(meV)",(((wf_sub(imode,iq)*ryd2meV,wf_sub(imode,iq)*ryd2meV),imode=1,nmodes),iq=1,nq)    
    do isnap=0,nsnap
      write(php_unit,"(F12.5,*(2(1X,E12.5)))") dt*nstep*isnap*ry_to_fs,((phPsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
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

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo

    write(phK_unit,"(A)") "Average of Normal mode kinetic energy for one core sample.SUM_phK,((phK(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(phK_unit,"(*(1X,A12))") "time","SUM_phK",(("K"//trim(adjustl(cphmode(imode,iq))),imode=1,nmodes),iq=1,nq)      
    !write(phK_unit,"(*(1X,A12))") "time ","SUM_phK",(("phK(mode,q)",imode=1,nmodes),iq=1,nq)
    write(phK_unit,"(*(1X,A12))") "fs ","  meV  ",(("   meV     ",imode=1,nmodes),iq=1,nq)
    write(phK_unit,"(2(1X,A12),*(1X,F12.5))") "Omega(meV)","SUM_phK",((wf_sub(imode,iq)*ryd2meV,imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
        write(phK_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,SUM(phKsit(:,:,isnap))*ryd2meV,&
        ((phKsit(imode,iq,isnap)*ryd2meV,imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phK_file_,phK_unit)
   
    write(stdout,"(A,A)") "Save average phK infortmation to the file: ",trim(phK_file_)
   
  end subroutine save_phK

  subroutine read_phK(inode,icore,nmodes,nq,nsnap,phKsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phKsit(nmodes,nq,0:nsnap)
    
    integer :: phk_unit
    character(len=maxlen) :: phK_file_
    character(len=maxlen) :: ctmp1,ctmp2
    real(kind=dp),allocatable :: phKsit_(:,:,:)
    
    if(.not. allocated(phKsit_)) allocate(phKsit_(nmodes,nq,0:nsnap))
    phKsit_ = 0.0
    
    write(ctmp1,*) inode
    write(ctmp2,*) icore
    phK_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phK_file))
    
    phk_unit = io_file_unit()
    call open_file(phK_file_,phk_unit)
    
    do i=1,4
      read(phK_unit,*)
    enddo    
    do isnap=0,nsnap
      read(phk_unit,"(26X,*(1X,E12.5))") ((phKsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    phKsit =phKsit + phKsit_
    
    call close_file(phK_file_,phk_unit)
   
  end subroutine read_phK  
  
  subroutine plot_phK(nmodes,nq,nsnap,phKsit)  
    use dynamics,only : bolziman
    implicit none
    
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phKsit(nmodes,nq,0:nsnap)
    real(kind=dp),allocatable :: phKsit_ave(:,:,:)
     
    
    real(kind=dp) :: womiga
    integer :: phK_unit
    character(len=maxlen) :: phK_file_   
    
    phK_file_ = trim(outdir)//trim(adjustl(phK_file))     
    phK_unit = io_file_unit()
    call open_file(phK_file_,phK_unit)

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo

    allocate(phKsit_ave(1:nmodes,1:nq,0:nsnap))     
    phKsit_ave = 0.0
    
    if(.not. allocated(nqv)) allocate(nqv(nmodes,nq))
    nqv = 0.0

    do iq=1,nq
      do imode=1,nmodes
        womiga = wf_sub(imode,iq)
        if(womiga >= eps_acustic) then
          nqv(imode,iq) = bolziman(womiga,temp)+0.5
        endif
      enddo
    enddo    
    
    if(l_ph_quantum) then
      do isnap =0,nsnap
        phKsit_ave(:,:,isnap) = phKsit(:,:,isnap) - 0.5*nqv*wf_sub*ryd2meV
      enddo
    else
      phKsit_ave = phKsit - 0.5*temp*K_B_Ryd*ryd2meV
    endif
    
    if(l_ph_quantum) then
      write(phK_unit,"(A)") "Average of Normal mode kinetic energy - 0.5*nqv*hbar*wqv for all trajecotry.SUM_phK,((phK(imode,iq),imode=1,nmodes),iq=1,nq)"
    else
      write(phK_unit,"(A)") "Average of Normal mode kinetic energy - 0.5*KB*T for all trajecotry.SUM_phK,((phK(imode,iq),imode=1,nmodes),iq=1,nq)"
    endif
    write(phK_unit,"(*(1X,A12))") "time ","phK(mode,q)",(("K"//trim(adjustl(cphmode(imode,iq))),imode=1,nmodes),iq=1,nq)
    write(phK_unit,"(*(1X,A12))") "fs ","  meV  ",(("   meV     ",imode=1,nmodes),iq=1,nq)
    write(phK_unit,"(2(1X,A12),*(1X,F12.5))") "Omega(meV)","SUM_phK",((wf_sub(imode,iq)*ryd2meV,imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
      write(phK_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,SUM(phKsit(:,:,isnap)),&
      ((phKsit(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phK_file_,phK_unit)
   
    deallocate(phKsit_ave)
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

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo
    
    write(phU_unit,"(A)") "Average of Normal mode potential energy for one core sample.SUM_phU,((phU(imode,iq),imode=1,nmodes),iq=1,nq)"
    write(phU_unit,"(*(1X,A12))") "time ","SUM_phU",(("U"//trim(adjustl(cphmode(imode,iq))),imode=1,nmodes),iq=1,nq)
    write(phU_unit,"(*(1X,A12))") "fs ","  meV  ",(("   meV     ",imode=1,nmodes),iq=1,nq)
    write(phU_unit,"(2(1X,A12),*(1X,F12.5))") "Omega(meV)","SUM_phU",((wf_sub(imode,iq)*ryd2meV,imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
        write(phU_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,SUM(phUsit(:,:,isnap))*ryd2meV,&
        ((phUsit(imode,iq,isnap)*ryd2meV,imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phU_file_,phU_unit)
   
    write(stdout,"(A,A)") "Save average phU infortmation to the file: ",trim(phU_file_)
   
  end subroutine save_phU  

  subroutine read_phU(inode,icore,nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: inode,icore,nmodes,nq,nsnap
    real(kind=dp),intent(inout) :: phUsit(nmodes,nq,0:nsnap)
    
    integer :: phu_unit
    character(len=maxlen) :: phU_file_
    character(len=maxlen) :: ctmp1,ctmp2
    real(kind=dp),allocatable :: phUsit_(:,:,:)
    
    if(.not. allocated(phUsit_)) allocate(phUsit_(nmodes,nq,0:nsnap))
    phUsit_ = 0.0
    
    write(ctmp1,*) inode
    write(ctmp2,*) icore
    phU_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(phU_file))
    phu_unit = io_file_unit()
    call open_file(phU_file_,phu_unit)
    
    do i=1,4
      read(phu_unit,*)
    enddo  
    do isnap=0,nsnap
      read(phu_unit,"(26X,*(1X,E12.5))") ((phUsit_(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    phUsit =phUsit + phUsit_
    
    call close_file(phU_file_,phu_unit)
   
  end subroutine read_phU  

  subroutine plot_phU(nmodes,nq,nsnap,phUsit)  
    integer,intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phUsit(nmodes,nq,0:nsnap)
    
    integer :: phU_unit
    character(len=maxlen) :: phU_file_  
    real(kind=dp),allocatable :: phUsit_ave(:,:,:)
    
    
    
    allocate(phUsit_ave(1:nmodes,1:nq,0:nsnap))
    phUsit_ave = 0.0
    
    phU_file_ = trim(outdir)//trim(adjustl(phU_file))     
    phU_unit = io_file_unit()
    call open_file(phU_file_,phU_unit)

    if(.not. allocated(cphmode)) allocate(cphmode(nmodes,nq))
    do iq=1,nq
      do imode=1,nmodes
        write(ctmp1,*) imode
        write(ctmp2,*) iq
        cphmode(imode,iq) = "("//trim(adjustl(ctmp1))//","//trim(adjustl(ctmp2))//")"
      enddo
    enddo

    if(l_ph_quantum) then
      do isnap =0,nsnap
        phUsit_ave(:,:,isnap) = phUsit(:,:,isnap) - 0.5*nqv*wf_sub*ryd2meV
      enddo
    else
      phUsit_ave = phUsit_ave - 0.5*temp*K_B_Ryd*ryd2meV
    endif
    
    if(l_ph_quantum) then
      write(phU_unit,"(A)") "Average of Normal mode potential energy - 0.5*nqv*hbar*wqv for all trajecotry.SUM_phU,((phU(imode,iq),imode=1,nmodes),iq=1,nq)"
    else
      write(phU_unit,"(A)") "Average of Normal mode potential energy - 0.5*KB*T for all trajecotry.SUM_phU,((phU(imode,iq),imode=1,nmodes),iq=1,nq)"
    endif
    write(phU_unit,"(*(1X,A12))") "time ","phU(mode,q)",(("U"//trim(adjustl(cphmode(imode,iq))),imode=1,nmodes),iq=1,nq)
    write(phU_unit,"(*(1X,A12))") "fs ","  meV  ",(("   meV     ",imode=1,nmodes),iq=1,nq)
    write(phU_unit,"(2(1X,A12),*(1X,F12.5))") "Omega(meV)","SUM_phU",((wf_sub(imode,iq)*ryd2meV,imode=1,nmodes),iq=1,nq)
    do isnap=0,nsnap
        write(phU_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,SUM(phUsit_ave(:,:,isnap)),&
        ((phUsit_ave(imode,iq,isnap),imode=1,nmodes),iq=1,nq)
    enddo
    
    call close_file(phU_file_,phU_unit)
    deallocate(phUsit_ave)
   
    write(stdout,"(A,A)") "Write average phU infortmation to the file: ",trim(phU_file_)
    
   
  end subroutine plot_phU  



  subroutine plot_ph_temp(nmodes,nq,nsnap,phKsit,phUsit)
    implicit none
    integer , intent(in) :: nmodes,nq,nsnap
    real(kind=dp),intent(in) :: phKsit(nmodes,nq,0:nsnap),phUsit(nmodes,nq,0:nsnap)
    
    real(kind=dp), allocatable :: phTemp(:,:),phEsit(:,:)
    
    real(kind=dp) :: qx,qy,qz
    
    integer :: phT_unit
    character(len=maxlen) :: phT_file_

    allocate(phTemp(nq,0:nsnap))
    allocate(phEsit(nq,0:nsnap))    
    do imode=1,nmodes
      write(ctmp1,*) imode
      phT_file_ = trim(outdir)//"ph_temp_mode"//trim(adjustl(ctmp1))//".dat"  
    
    
      phT_unit = io_file_unit()
      call open_file(phT_file_,phT_unit)    

      phTemp = 0.0
      phEsit = 0.0
    
      if(l_ph_quantum) then
        do isnap=0,nsnap
          do iq=1,nq
            if(wf_sub(imode,iq)>0.0) then
              phEsit(iq,isnap) = (phKsit(imode,iq,isnap)+phUsit(imode,iq,isnap))/(wf_sub(imode,iq)*ryd2meV)
            else
              phEsit(iq,isnap) = 0.0
            endif
          enddo
        enddo
      
        do isnap=0,nsnap
          do iq=1,nq
            if(phEsit(iq,isnap)<=0.5) then
              phTemp(iq,isnap) = 0.0
            else
              phTemp(iq,isnap) = (wf_sub(nmodes,iq)/K_B_Ryd)/log((phEsit(iq,isnap)+0.5)/(phEsit(iq,isnap)-0.5))
            endif
          enddo
        enddo
      else
          phEsit = phKsit(imode,:,:)+phUsit(imode,:,:)
          phTemp = phEsit/ryd2meV/K_B_Ryd     
      endif
    
      if(l_ph_quantum) then
        write(phT_unit,"(A29,I12,A)") "The temperature of the imode=",imode," phonon use nqv."
      else
        write(phT_unit,"(A29,I12,A)") "The temperature of the imode=",imode," phonon use KbT."
      endif
    
      write(phT_unit,"(6(1X,A12))") "time","iq","qx","qy","qz","temperature"
      write(phT_unit,"(6(1X,A12))") " fs ","iq","b1","b2","b3"," K "
      do isnap=0,nsnap
        do iq=1,nq
          qx=xqf(1,indexq(iq))
          qy=xqf(2,indexq(iq))
          qz=xqf(3,indexq(iq))
          if(qx>0.5) qx=qx-1.0
          if(qy>0.5) qy=qy-1.0
          if(qz>0.5) qz=qz-1.0
          write(phT_unit,"(1X,F12.5,1X,I12,4(1X,F12.5))") dt*nstep*isnap*ry_to_fs,iq,qx,qy,qz,phTemp(iq,isnap)
          if(qx== 0.5 .or. qy==0.5 .or. qz==0.5) then
            if(qx==0.5) qx=qx-1.0
            if(qy==0.5) qy=qy-1.0
            if(qz==0.5) qz=qz-1.0
            write(phT_unit,"(1X,F12.5,1X,I12,4(1X,F12.5))") dt*nstep*isnap*ry_to_fs,iq,qx,qy,qz,phTemp(iq,isnap)
          endif
        enddo
      enddo
      
      call close_file(phT_file_,phT_unit)
    
    enddo
    
  end subroutine plot_ph_temp



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

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo

    write(wsit_unit,"(A)") "The average electron(hole) population on different adiabatic wave function for one core sample. "    
    write(wsit_unit,"(*(1X,A12))") "time ",("wsit"//trim(adjustl(cefree(ifre))),ifre=1,nfre)
    write(wsit_unit,"((1X,A12),*(1X,I12))") "fs  ",(ifre,ifre=1,nfre)
    do isnap=0,nsnap
      write(wsit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(wsit(ifre,isnap),ifre=1,nfre)  
    enddo    

    
    call close_file(wsit_file_,wsit_unit)

    write(stdout,"(A,A)")"Save the average electron(hole) wavefction population on &
              &different adiabatic wave function to file:",trim(wsit_file_)
  
    deallocate(cefree)
    
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
    
    read(wsit_unit,*)
    read(wsit_unit,*)
    read(wsit_unit,*)    
    do isnap=0,nsnap
      read(wsit_unit,"(13X,*(1X,E12.5))") (wsit_(ifre,isnap),ifre=1,nfre)   
    enddo
    
    wsit = wsit + wsit_
    
    call close_file(wsit_file,wsit_unit)
  
  end subroutine read_wsit  

  subroutine plot_wsit(nfre,nfre_sh,nsnap,naver,wsit,wsit_file)
    implicit none
    integer,intent(in) :: nfre,nfre_sh,nsnap,naver
    real(kind=dp),intent(in) :: wsit(nfre,0:nsnap)
    character(len=*),intent(in) :: wsit_file
    
    integer :: wsit_unit
    character(len=maxlen) :: wsit_file_  
    wsit_file_ = trim(outdir)//trim(adjustl(wsit_file))     
    wsit_unit = io_file_unit()
    call open_file(wsit_file_,wsit_unit)

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo

    write(wsit_unit,"(A)") "The average electron(hole) population on different adiabatic wave function for all trajecotry. "    
    write(wsit_unit,"(*(1X,A12))")          "time ",("wsit"//trim(adjustl(cefree(ifre))),ifre=1,nfre_sh)
    write(wsit_unit,"((1X,A12),*(1X,A12))") "fs  ",("|w|2",ifre=1,nfre_sh)
    do isnap=0,nsnap
      write(wsit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(wsit(ifre,isnap),ifre=1,nfre_sh)  
    enddo
    
    call close_file(wsit_file_,wsit_unit)
    write(stdout,"(A,A)")"Write the average electron(hole) wavefction population on &
              &different adiabatic wave function to file:",trim(wsit_file_)
  
    deallocate(cefree)
  end subroutine plot_wsit



  ! save the first trajecotry active PES and first trajecotry PES : pes_e.dat_f.gnu
  ! save the average active PES and PES for all trajecotry. : pes_e.dat_average.gnu
  subroutine save_pes(nfre,nsnap,naver,pes_one,pes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: pes_one(0:nfre,0:nsnap),pes(0:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file
    
    character(len=maxlen) :: pes_file_
    integer :: pes_unit
    pes_file_= trim(outdir)//trim(adjustl(pes_file))
    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo

    write(pes_unit,"(A,I8,A)") "Averager of naver=",naver, " pes(ifre,isnap),ifre=0,nfre)"
    write(pes_unit,"(*(1X,A12))") "time ","active_pes",("pes"//trim(adjustl(cefree(ifre))),ifre=1,nfre)
    write(pes_unit,"(*(1X,A12))") "fs",("eV",ifre=0,nfre)    
    do isnap=0,nsnap
      write(pes_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(pes(ifre,isnap)*RYTOEV,ifre=0,nfre)
    enddo

    write(stdout,"(A,A)") "Save average of active PES and PES for all trajecotry to the file:",trim(pes_file_)

    
    call close_file(pes_file_,pes_unit)

    pes_file_= trim(outdir)//trim(adjustl(pes_file))//".first.dat"
    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)

    write(pes_unit,"(A,I8,A)") "One trajecotry PES of iaver=",1, " pes_one(ifre,isnap),ifre=0,nfre)"
    write(pes_unit,"(*(1X,A12))") "time ","active_pes",("pesone(ifre)",ifre=1,nfre)
    write(pes_unit,"(*(1X,A12))") "fs",("eV",ifre=0,nfre)    
    do isnap=0,nsnap
      write(pes_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(pes_one(ifre,isnap)*RYTOEV,ifre=0,nfre)
    enddo

    call close_file(pes_file_,pes_unit)    
    write(stdout,"(A,A)") "Save the first trajecotry active PES and PES to the file:",trim(pes_file_)
    
    deallocate(cefree)
    
  end subroutine save_pes
  
  subroutine read_pes(inode,icore,nfre,nsnap,naver,pes_one,pes,pes_file)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver,inode,icore
    real(kind=dp),intent(inout) :: pes_one(0:nfre,0:nsnap),pes(0:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file

    
    integer :: pes_unit
    character(len=maxlen) :: pes_file_
    character(len=maxlen) :: ctmp1,ctmp2
    
    real(kind=dp),allocatable :: pes_(:,:)
    if(.not. allocated(pes_)) allocate(pes_(0:nfre,0:nsnap))
    
    write(ctmp1,*) inode
    write(ctmp2,*) icore
    pes_file_ = "./node"//trim(adjustl(ctmp1))//"/sample"//trim(adjustl(ctmp2))//"/"//trim(outdir)//trim(adjustl(pes_file))    
    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    
    read(pes_unit,*)
    read(pes_unit,*)
    read(pes_unit,*)
    do isnap=0,nsnap
      read(pes_unit,"(13X,*(1X,E12.5))") (pes_(ifre,isnap),ifre=0,nfre)
    enddo
    
    pes = pes + pes_
    
    call close_file(pes_file_,pes_unit)    
    
    pes_file_ = trim(adjustl(pes_file_))//".first.dat"
    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    
    read(pes_unit,*)
    read(pes_unit,*)
    read(pes_unit,*)
    do isnap=0,nsnap
      read(pes_unit,"(13X,*(1X,E12.5))") (pes_(ifre,isnap),ifre=0,nfre)
    enddo
    
    pes_one = pes_
    
    call close_file(pes_file_,pes_unit)        
    
    
  end subroutine read_pes
  
  subroutine plot_pes(nfre,nfre_sh,nsnap,pes_one,pes,wsit,pes_file)
    implicit none
    integer , intent(in) :: nfre,nfre_sh,nsnap
    real(kind=dp),intent(in) :: pes_one(0:nfre,0:nsnap),pes(0:nfre,0:nsnap),wsit(1:nfre,0:nsnap)
    character(len=*),intent(in) :: pes_file
    
    character(len=maxlen) :: pes_file_
    integer :: pes_unit
    pes_file_ = trim(outdir)//trim(adjustl(pes_file))//"_f.dat"

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    write(pes_unit,"(A)") "USing the first trajecotry as an example to Plotting potential Energy Surface(PES),and the active PES."
    write(pes_unit,"(*(1X,A12))") "time ","active_fpes",("pes"//trim(adjustl(cefree(ifre))),ifre=1,nfre_sh)
    write(pes_unit,"(*(1X,A12))") "fs  ","   eV     "  ,("  eV  ",ifre=1,nfre_sh)
    do isnap=0,nsnap
      write(pes_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(pes_one(ifre,isnap),ifre=0,nfre_sh)
    enddo
    
    call close_file(pes_file_,pes_unit)
    
    write(stdout,"(A,A)") "Write the first trajecotry active PES and PES to the file:",trim(pes_file_)

    pes_file_ = trim(outdir)//trim(adjustl(pes_file))//"_average.dat"

    pes_unit = io_file_unit()
    call open_file(pes_file_,pes_unit)
    write(pes_unit,"(A)") "Plotting averager of potential Energy Surface(PES),and the averager active PES and their weight."
    write(pes_unit,"(*(1X,A12))") "time ","ave_PES","wsit      ","avePES(a)","avePES(1)","avePES(nfre)"
    write(pes_unit,"(*(1X,A12))") " fs  ","   eV  ","w(ifre)**2","    eV   ","    eV   ","     eV     "

    do ifre =1 , nfre_sh
      !write(pes_unit,"(A,I12,1X,A12)") "#isurface=",ifre,"wsit"
      do isnap=0,nsnap
        if(ifre ==1 ) then
          write(pes_unit,"(6(1X,E12.5))") dt*nstep*isnap*ry_to_fs,pes(ifre,isnap),wsit(ifre,isnap),&
          pes(0,isnap),pes(1,isnap),pes(nfre,isnap)
        else
          write(pes_unit,"(3(1X,E12.5))") dt*nstep*isnap*ry_to_fs,pes(ifre,isnap),wsit(ifre,isnap)
        endif
      enddo
      write(pes_unit,*)
    enddo
    
    call close_file(pes_file_,pes_unit)
    
    write(stdout,"(A,A)") "Write average of active PES and PES for all trajecotry to the file:",trim(pes_file_)
    
    deallocate(cefree)
    
  end subroutine plot_pes
  


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

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo
    
    write(csit_unit,"(A)") "The average electron(hole) population on different diabatic wave function for one core sample. "
    write(csit_unit,"(*(1X,A12))") "time",("csit"//trim(adjustl(cefree(ifre))),ifre=1,nfre)
    write(csit_unit,"((1X,A12),*(1X,I12))") "fs ",(ifre,ifre=1,nfre)
    do isnap=0,nsnap
      write(csit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file_,csit_unit)
    write(stdout,"(A,A)")"Save the average electron(hole) wavefction population on &
              &different diabatic wave function to file:",trim(csit_file_)
  
    deallocate(cefree)
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

    read(csit_unit,*)
    read(csit_unit,*)
    read(csit_unit,*)    
    do isnap=0,nsnap
      read(csit_unit,"(13X,*(1X,E12.5))") (csit_(ifre,isnap),ifre=1,nfre)   
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
    csit_file_ = trim(outdir)//trim(adjustl(csit_file)) 
    csit_unit = io_file_unit()
    call open_file(csit_file_,csit_unit)

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo
    
    write(csit_unit,"(A)") "The electron(hole) population on different diabatic wave function. "
    write(csit_unit,"(*(1X,A12))") "time",("csit"//trim(adjustl(cefree(ifre))),ifre=1,nfre)
    write(csit_unit,"((1X,A12),*(1X,A12))") "fs ",("|c|2",ifre=1,nfre)
    do isnap=0,nsnap
      write(csit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(csit(ifre,isnap),ifre=1,nfre)   
    enddo
    
    call close_file(csit_file_,csit_unit)
    
    write(stdout,"(A,A)")"Write the average electron(hole) wavefction population on &
              &different diabatic wave function to file:",trim(csit_file_)
    
    deallocate(cefree)
  end subroutine plot_csit



  !The active adiabatic PES project to diabatic states. psit(ifre)=|P(ifre,isurface)|**2
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

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo

    write(psit_unit,"(A)") "The average active PES project on different diabatic wave function for one core sample. psit(ifre)=P(ifre,isurface)**2 "    
    write(psit_unit,"(*(1X,A12))") "time",("psit"//trim(adjustl(cefree(ifre))),ifre=1,nfre)
    write(psit_unit,"((1X,A12),*(1X,I12))") "fs ",(ifre,ifre=1,nfre)
    do isnap=0,nsnap
      write(psit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(psit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(psit_file_,psit_unit)

    write(stdout,"(A,A)")"Save the average electron(hole) active PES project to &
              &different diabatic wave function to file:",trim(psit_file_)
  
    deallocate(cefree)
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

    read(psit_unit,*)
    read(psit_unit,*)
    read(psit_unit,*)
    do isnap=0,nsnap
      read(psit_unit,"(13X,*(1X,E12.5))") (psit_(ifre,isnap),ifre=1,nfre)   
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
    psit_file_ = trim(outdir)//trim(adjustl(psit_file))         
    psit_unit = io_file_unit()
    call open_file(psit_file_,psit_unit)

    if(.not. allocated(cefree)) allocate(cefree(nfre))
    do ifre=1,nfre
      write(ctmp1,*) ifre
      cefree(ifre) = "("//trim(adjustl(ctmp1))//")"
    enddo
    
    write(psit_unit,"(A)") "The average active PES project on different diabatic wave function for all trajecotry. psit(ifre)=P(ifre,isurface)**2 "    
    write(psit_unit,"(*(1X,A12))") "time",("psit"//trim(adjustl(cefree(ifre))),ifre=1,nfre)
    write(psit_unit,"((1X,A12),*(1X,A12))") "fs ",("|P|2",ifre=1,nfre)
    do isnap=0,nsnap
      write(psit_unit,"(*(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(psit(ifre,isnap),ifre=1,nfre)  
    enddo
    
    call close_file(psit_file_,psit_unit)
  
    write(stdout,"(A,A)")"Write the average electron(hole) active PES project to &
              &different diabatic wave function to file:",trim(psit_file_)
  
    deallocate(cefree)
  end subroutine plot_psit  


  subroutine plot_band_occupatin_withtime(nband,nk,Enk,nsnap,psit,csit,dsnap,band_file)
    use elph2,only : xkf
    implicit none
    integer , intent(in) :: nband,nk,dsnap
    real(kind=dp),intent(in) :: Enk(nband,nk)
    integer , intent(in) :: nsnap
    real(kind=dp),intent(in) :: psit(1:nband*nk,0:nsnap),csit(1:nband*nk,0:nsnap)
    character(len=*),intent(in) :: band_file
    
    character(len=maxlen) :: band_file_
    integer :: band_unit
    integer :: ipol,iband,ik,ifre
    real(kind=dp) :: kx,ky,kz
    
    do iband=1,nband
      write(ctmp1,*) iband
      band_file_=trim(outdir)//trim(adjustl(band_file))//"_"//trim(adjustl(ctmp1))//".dat"  
      band_unit = io_file_unit()
      call open_file(band_file_,band_unit)
      write(band_unit,"(A32,I12,A)") "Carrier occupation on the iband=",iband," structure at different time"
      write(band_unit,"(*(1X,A12))")    "iband","ik","kx","ky","kz","E_nk",("psit","csit",isnap=0,nsnap,dsnap)
      write(band_unit,"(6(1X,A12),*(1X,F12.2))") "index_band"," index_k","b1","b2", "b3", " eV" ,&
                                          (dt*nstep*isnap*ry_to_fs,dt*nstep*isnap*ry_to_fs,isnap=0,nsnap,dsnap)
      do ik=1,nk
        ifre = (ik-1)*nband+iband
        kx = xkf(1,indexk(ik))
        ky = xkf(2,indexk(ik))
        kz = xkf(3,indexk(ik))
        if(kx > 0.5) kx=kx-1.0
        if(ky > 0.5) ky=ky-1.0
        if(kz > 0.5) kz=kz-1.0
        write(band_unit,"(2(1X,I12),*(1X,F12.5))") iband,ik,kx,ky,kz,Enk(iband,ik)*RYTOEV,(psit(ifre,isnap),csit(ifre,isnap),isnap=0,nsnap,dsnap)
        if(kx == 0.5 .or. ky==0.5 .or. kz==0.5) then
          if(kx==0.5) kx=kx-1.0
          if(ky==0.5) ky=ky-1.0
          if(kz==0.5) kz=kz-1.0
          write(band_unit,"(2(1X,I12),*(1X,F12.5))") iband,ik,kx,ky,kz,Enk(iband,ik)*RYTOEV,(psit(ifre,isnap),csit(ifre,isnap),isnap=0,nsnap,dsnap)
        endif
      enddo
    
      call close_file(band_file_,band_unit)
    enddo
    write(stdout,"(A,A)") "Write Carrier occupation on the band structure at different time to file:",trim(band_file_)
  
  end subroutine plot_band_occupatin_withtime

  
  subroutine plot_eh_temp(nband,nk,Enk,nsnap,psit,csit,dsnap,temp_file)
    use elph2,only : xkf
    implicit none
    integer , intent(in) :: nband,nk,dsnap
    real(kind=dp),intent(in) :: Enk(nband,nk)
    integer , intent(in) :: nsnap
    real(kind=dp),intent(in) :: psit(1:nband*nk,0:nsnap),csit(1:nband*nk,0:nsnap)
    character(len=*),intent(in) :: temp_file
    
    character(len=maxlen) :: temp_file_
    integer :: temp_unit
    integer :: ipol,iband,ik,ifre
    real(kind=dp) :: kx,ky,kz
    
    do iband=1,nband
      write(ctmp1,*) iband
      temp_file_=trim(outdir)//trim(adjustl(temp_file))//"_"//trim(adjustl(ctmp1))//".dat"  
      temp_unit = io_file_unit()
      call open_file(temp_file_,temp_unit)
      
      write(temp_unit,"(A34,I12,A)") "Carrier occupation on the iband = ",iband," structure at different time."
      write(temp_unit,"(9(A12,1X))") " time "," iband "," ik "," kx "," ky "," kz "," E_nk ","psit "," csit "
      write(temp_unit,"(9(A12,1X))") " fs "," iband "," ik "," b1 "," b2 ", " b3 "," eV "," psit "," csit " 
  
      do isnap =0,nsnap
        do ik=1,nk
          ifre = (ik-1)*nband+iband
          kx = xkf(1,indexk(ik))
          ky = xkf(2,indexk(ik))
          kz = xkf(3,indexk(ik))
          if(kx > 0.5) kx=kx-1.0
          if(ky > 0.5) ky=ky-1.0
          if(kz > 0.5) kz=kz-1.0
          write(temp_unit,"(F12.5,2(1X,I12),*(1X,F12.5))") &
          dt*nstep*isnap*ry_to_fs,iband,ik,kx,ky,kz,Enk(iband,ik)*RYTOEV,psit(ifre,isnap),csit(ifre,isnap)
          if(kx == 0.5 .or. ky==0.5 .or. kz==0.5) then
            if(kx==0.5) kx=kx-1.0
            if(ky==0.5) ky=ky-1.0
            if(kz==0.5) kz=kz-1.0
            write(temp_unit,"(F12.5,2(1X,I12),*(1X,F12.5))") &
            dt*nstep*isnap*ry_to_fs,iband,ik,kx,ky,kz,Enk(iband,ik)*RYTOEV,psit(ifre,isnap),csit(ifre,isnap)
          endif
        enddo
        write(temp_unit,*)
      enddo
      call close_file(temp_file_,temp_unit)
    enddo
    write(stdout,"(A,A)") "Write Carrier occupation on the band structure at different time to file:",trim(temp_file_)    
    
  end subroutine plot_eh_temp
  
  
end module saveinf