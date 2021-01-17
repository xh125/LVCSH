module environments
  use kinds, only: dp
  use global_version, only : version_number
  use io,only : stdout,stdout_name,io_file_unit
  use date_and_times,only : get_date_and_time
  use constants, only : constants_version_str1,constants_version_str2
  implicit none
  character(len=80) :: code_version
  integer :: max_threads
  real(kind=dp) :: time0
  real(kind=dp) :: t0
  character(len=80) :: mkl_version
  logical :: lsetthreads  ! set mkl_threads by yourself 
  integer :: mkl_threads
  
  contains
  
  subroutine environment_start( code )
    use mkl_service
    implicit none
    character(len=*),intent(in) :: code
    
    code_version = trim(adjustl(code)) // "v." // trim(adjustl(version_number))
    stdout = io_file_unit()
    stdout_name = "LVCSH.out"
    open(unit=stdout,file=stdout_name,status='REPLACE')
    !call io_date(cdate,ctime)
    call mkl_get_version_string(mkl_version )
    max_threads = mkl_get_max_threads()
    
    call openning_message( code_version)
    
  end subroutine environment_start

  subroutine openning_message( code_version)
    use kinds,only : print_kind_info
    implicit none
    
    character(len=*),intent(in) :: code_version
    character(len=9) :: cdate,ctime
    character(len=76) :: ctmp
    
    call get_date_and_time(cdate,ctime)
    
    write(stdout,"(1X,A77)") repeat("=",77)
    write(stdout,*) "LVCSH complied with using Intel Fortran and MKL:"
    write(stdout,*) trim(adjustl(mkl_version))
    write(stdout,"(1X,A77)") repeat("=",77)
    
    write(stdout,"(/,1X,A77)") repeat("=",77)
    write(stdout,'(1X,"Program ",A," startes on ",A9," at ",A9)') &
                  &trim(code_version),cdate,ctime
    write(ctmp  ,*) max_threads
    write(stdout,*) "By default, Intel MKL uses "//trim(adjustl(ctmp))//" threads"
    write(stdout,*) "where "//trim(adjustl(ctmp))//" is the number of physical cores on the system"
    write(stdout,"(1X,A77)") repeat("=",77)
    
    call print_kind_info (stdout)
    write(stdout,"(/,1X,A77)") repeat("=",77)
    write(stdout,*) constants_version_str1
    write(stdout,*) constants_version_str2
    write(stdout,"(1X,A77)") repeat("=",77)
    
    return
  end subroutine openning_message
  
  subroutine closing_message()
    implicit none
    character(len=9) :: cdate,ctime
    character(len=80) :: time_str
    
    call get_date_and_time(cdate,ctime)
    time_str = 'This run was terninated on : ' // ctime // ' ' // cdate
    
    write(stdout,*)
    write(stdout,"(3X,A60,/)") time_str
    write(stdout,"( '=',78('-'),'=' )" )
    
    return
  end subroutine closing_message
  
  subroutine compilation_info()
  !
  ! code borrowed by WanT - prints architecture / compilation details
  !
#if defined(__HAVE_CONFIG_INFO)
#include "configure.h"
  write( stdout, * )
  write( stdout, "( 2x,'        ARCH :',4x,a )") trim(adjustl(__conf_ARCH))
  write( stdout, "( 2x,'          CC :',4x,a )") trim(adjustl(__conf_CC))
  write( stdout, "( 2x,'         CPP :',4x,a )") trim(adjustl(__conf_CPP))
  write( stdout, "( 2x,'      MPIF90 :',4x,a )") trim(adjustl(__conf_MPIF90))
  write( stdout, "( 2x,'         F77 :',4x,a )") trim(adjustl(__conf_F77))
  write( stdout, "( 2x,'      DFLAGS :',4x,a )") trim(adjustl(__conf_DFLAGS))    
  write( stdout, "( 2x,'   BLAS_LIBS :',4x,a )") trim(adjustl(__conf_BLAS_LIBS))
  write( stdout, "( 2x,' LAPACK_LIBS :',4x,a )") trim(adjustl(__conf_LAPACK_LIBS))
  write( stdout, "( 2x,'    FFT_LIBS :',4x,a )") trim(adjustl(__conf_FFT_LIBS))
  write( stdout, "( 2x,'   MASS_LIBS :',4x,a )") trim(adjustl(__conf_MASS_LIBS))
#endif
  end subroutine compilation_info

  subroutine set_mkl_threads(mkl_threads)
    implicit none
    integer,intent(in) :: mkl_threads
    call mkl_set_num_threads(mkl_threads)
    write(stdout,*) repeat("=",57)
    write(stdout,"(5X,A,i3,A)") "Reset Intel MKL uses ",max_threads," threads." 
    write(stdout,*) 
  end subroutine set_mkl_threads

end module environments