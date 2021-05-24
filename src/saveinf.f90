module saveinf
  use kinds,only : dp
  use surfacecom,only : dt,nstep
  use io,only : io_file_unit,open_file,close_file
  use constants,only : ry_to_fs
  implicit none
  integer :: iaver,isnap,ifre 
  
  character(len=9) :: pes_e_file="pes_e.dat"
  character(len=9) :: pes_h_file="pes_h.dat"
  
  contains 
  
  subroutine save_pes(nfre,nsnap,naver,pes,pes_filename)
    implicit none
    integer , intent(in) :: nfre,nsnap,naver
    real(kind=dp),intent(in) :: pes(0:nfre,1:nsnap,1:naver)
    character(len=*),intent(in) :: pes_filename
    
    integer :: pes_unit

    pes_unit = io_file_unit()
    call open_file(pes_filename,pes_unit)
    do iaver =1 , 1
      do isnap=1,nsnap
        write(pes_unit,"(5X,F11.3,7(1X,E12.5))") dt*nstep*isnap*ry_to_fs,(pes(ifre,isnap,iaver),ifre=0,nfre)
      enddo
    enddo
 
    call close_file(pes_filename,pes_unit)
    
  end subroutine save_pes
  
end module saveinf