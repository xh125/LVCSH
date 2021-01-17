module softinformation
  use io
  use mkl_service
  implicit none
  
  contains
  
  subroutine stdoutheader()
    use kinds
    use constants
    use parameters,only : ctmp,cdate,ctime
    
    implicit none
    
    integer           :: max_threads
    real(kind=dp)     :: time0  
    real(kind=dp)     :: t0
    
    time0   = io_time()
    stdout  = io_file_unit()
    stdout_name = "SCSH.out"
    open(unit=stdout,file=stdout_name,status='REPLACE')
    !open(unit=stdout,file="SCSH.out",status='REPLACE')
    call io_date(cdate,ctime)
    call mkl_get_version_string( ctmp )
    max_threads = mkl_get_max_threads()
    
    write(stdout,"(1X,A67)") repeat("=",67)
    write(stdout,*) "SCSH complied with using IMKL:"
    write(stdout,*) trim(adjustl(ctmp))
    write(stdout,"(1X,A67)") repeat("=",67)
    
    write(stdout,"(/,1X,A67)") repeat("=",67)
    write(stdout,*) 'SCSH :Execution started on ',cdate,' at ',ctime
    write(ctmp  ,*) max_threads
    write(stdout,*) "By default, Intel MKL uses "//trim(adjustl(ctmp))//" threads"
    write(stdout,*) "where "//trim(adjustl(ctmp))//" is the number of physical cores on the system"
    write(stdout,"(1X,A67)") repeat("=",67)
    
    call print_kind_info (stdout)
    
    write(stdout,"(/,1X,A67)") repeat("=",67)
    write(stdout,*) constants_version_str1
    write(stdout,*) constants_version_str2
    write(stdout,"(1X,A67)") repeat("=",67)
    
  end subroutine stdoutheader
  
  subroutine set_mkl_threads(mkl_threads)
    use parameters,only: ctmp
    implicit none
    integer,intent(in) :: mkl_threads
    call mkl_set_num_threads(mkl_threads)
    write(ctmp,*) mkl_threads
    write(stdout,*) repeat("=",57)
    write(stdout,*) "Reset Intel MKL uses "//trim(adjustl(ctmp))//" threads." 
    write(stdout,*)
  end subroutine set_mkl_threads
    
end module softinformation
