module elphmat
use kinds,only : dp
use constants,only : maxlen
implicit none
  integer :: nk1,nk2,nk3,nq1,nq2,nq3
  integer :: nkf1,nkf2,nkf3,nqf1,nqf2,nqf3
  integer :: nmode,nu,nqftot,nkftot
  integer :: ibndmin,ibndmax
  real(kind=dp) :: ebndmin, ebndmax
  real(kind=dp),allocatable :: wf(:,:)         ! omega(nmode,nqftot)
  real(kind=dp),allocatable :: xqf(:,:)        ! xqf(3,nqftot)
  real(kind=dp),allocatable :: xkf(:,:)        !
  real(kind=dp),allocatable :: enkf(:,:)       ! En(ibndmin:ibndmax,nkftot)  
  real(kind=dp),allocatable :: epmatb(:,:,:,:,:)
  
  character(len=maxlen) :: ctmp1,ctmp2,ctmp3
  
  contains
  
  subroutine readepwout(filepwout)
    use io,only : io_file_unit,open_file,close_file,findkword,findkline
    implicit none
    character(len=*), intent(in) :: filepwout
    integer :: epw_unit
    character(len=maxlen) :: epw_name
    
    epw_unit=io_file_unit()
    epw_name=trim(adjustl(filepwout))
    call open_file(epw_name,epw_unit)
    
    call findkline(epw_unit,"     Using uniform q-mesh:",1,26)
    read(epw_unit,*) ctmp1,ctmp2,ctmp3,nqf1,nqf2,nqf3
    
    
  end subroutine readepwout

end module 