module environments
  use kinds, only: dp
  use global_version, only : version_number
  use io,only : stdout,stdout_name,io_file_unit
  use date_and_times,only : get_date_and_time
  use constants, only : maxlen,constants_version_str1,constants_version_str2
  implicit none
  character(len=maxlen) :: code_version
  integer :: max_threads
  real(kind=dp) :: time0
  real(kind=dp) :: t0
  character(len=maxlen) :: mkl_version
  logical :: lsetthreads  ! set mkl_threads by yourself 
  integer :: mkl_threads
  
  contains
  
  subroutine environment_start( code )
    use mkl_service
    implicit none
    character(len=*),intent(in) :: code
    
    code_version = trim(adjustl(code)) // " v." // trim(adjustl(version_number))
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
    character(len=maxlen) :: ctmp
    
    call get_date_and_time(cdate,ctime)
    
    ! Display the logo. Use https://www.asciiarts.net/ 
    write(stdout,'(10X,a)') "  _____  ____   ____   ______   ______   ____  ____  "
    write(stdout,'(10X,a)') " |_   _||_  _| |_  _|.' ___  |.' ____ \ |_   ||   _| "
    write(stdout,'(10X,a)') "   | |    \ \   / / / .'   \_|| (___ \_|  | |__| |   "
    write(stdout,'(10X,a)') "   | |   _ \ \ / /  | |        _.____`.   |  __  |   "
    write(stdout,'(10X,a)') "  _| |__/ | \ ' /   \ `.___.'\| \____) | _| |  | |_  "
    write(stdout,'(10X,a)') " |________|  \_/     `.____ .' \______.'|____||____| "    
    write(stdout,'(10X,a)') "                                                     "
    write(stdout,'(10X,a)') "                                                     "
    write(stdout,'(20X,a)') " Developed by XieHua at department of physic, USTC   "    
    write(stdout,'(20X,a)') "                        Email:xh125@mail.ustc.edu.cn "    
    write(stdout,'(a)') "                                                     "     
    
    !write(stdout,"(1X,A77)") repeat("=",77)
    write(stdout,"(1X,A)") "LVCSH complied with using Intel Fortran Complier and MKL:"
    !write(stdout,"(1X,A)") trim(adjustl(mkl_version))
    write(stdout,"(1X,A73)") mkl_version(1:73)
    write(stdout,"(1X,A)")   mkl_version(74:)
    !write(stdout,"(1X,A77)") repeat("=",77)
    write(ctmp  ,*) max_threads
    write(stdout,*) "By default, Intel MKL uses "//trim(adjustl(ctmp))//" threads"
    write(stdout,*) "where "//trim(adjustl(ctmp))//" is the number of physical cores on the system"


    call print_kind_info (stdout)
    write(stdout,"(/,1X,A77)") repeat("=",77)
    write(stdout,*) constants_version_str1
    write(stdout,*) constants_version_str2
    write(stdout,"(1X,A77)") repeat("=",77)
    
    !write(stdout,"(/,1X,A77)") repeat("=",77)
    write(stdout,'(/,1X,"Program ",A," startes on ",A9," at ",A9)') &
                  &trim(code_version),cdate,ctime
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
  

  subroutine set_mkl_threads(mkl_threads)
    implicit none
    integer,intent(in) :: mkl_threads
    if(mkl_threads <= max_threads) then
      call mkl_set_num_threads(mkl_threads)
    else
      write(stdout,"(5x,A)") "The mkl_threads is setting wrong!!"
      write(stdout,"(5x,A,I,A,I)"),"The max_threads is:",max_threads,"mkl_threads need <",max_threads
    endif
    write(stdout,*) repeat("=",77)
    write(stdout,"(5X,A,i3,A)") "Reset Intel MKL uses ",mkl_threads," threads." 
    write(stdout,*) 
  end subroutine set_mkl_threads

end module environments