module io
  !! Module to handle operations related to file input and output.
  use kinds     ,only : dp
  use constants ,only : maxlen
  implicit none
  
  integer, public, save           :: stdout,stdin
  character(len=maxlen)           :: stdout_name,stdin_name
  !! Unit on which stdout is written
  !integer, parameter, public      :: maxlen = 120  
  !! Max column width of input file
  integer                         :: ierr
  character(len=maxlen)           :: msg
  character(len=maxlen)           :: home_dir
  ! For parallel execution: I/O within an image
  ! These are set at startup by calling mp_world_start
  integer :: ionode_id= 0       ! index of the i/o node for this image
  Logical :: Lionode  = .True.  ! true if this processor is a i/o node
                                   ! for this image
                                   
  type timing_data
    !! Data about each stopwatch - for timing routines
    integer :: ncalls           
    !! Number of times stopwatch has been called
    real(kind=DP) :: ctime      
    !! Total time on stopwatch
    real(kind=DP) :: ptime
    !! Temporary record of time when watch is started     
    character(len=60) :: label      
    !! What is this stopwatch timing
  end type timing_data

  integer, parameter :: nmax = 100
  !! Maximum number of stopwatches
  type(timing_data) :: clocks(nmax)
  !! Data for the stopwatches 
  integer, save     :: nnames=0
  !! Number of active stopwatches
  
contains
  
  subroutine io_stopwatch(tag,mode)
  !=====================================
  !! Stopwatch to time parts of the code
  !=====================================  
  implicit none

    character(len=*), intent(in) :: tag
    !! Which stopwatch to act upon
    integer, intent(in)          :: mode
    !! Action  1=start 2=stop

    integer :: i
    real(kind=dp) :: t
    
    call cpu_time(t) !返回当前程序的处理机时间，单位为秒 参考P.740

    select case (mode)
      case (1)
        do i=1,nnames
          if (clocks(i)%label .eq. tag) then
            clocks(i)%ptime  = t
            clocks(i)%ncalls = clocks(i)%ncalls + 1
            return
          endif
        enddo
        nnames = nnames + 1
        if (nnames.gt.nmax) call io_error('Maximum number of calls to io_stopwatch exceeded')
        clocks(nnames)%label = tag
        clocks(nnames)%ctime = 0.0_dp
        clocks(nnames)%ptime = t
        clocks(nnames)%ncalls = 1
      case (2)
        do i=1,nnames
          if (clocks(i)%label .eq. tag) then
            clocks(i)%ctime = clocks(i)%ctime + t - clocks(i)%ptime
            return   !正确的读取了秒表信息，返回调用程序
          endif
        end do
        write(stdout,'(1x,3a)') 'WARNING: name = ',trim(tag),' not found in io_stopwatch' 
      case default
        write(stdout,*) ' Name = ',trim(tag),' mode = ',mode
        call io_error('Value of mode not recognised in io_stopwatch')
    end select
    
    return

  end subroutine io_stopwatch
  
  
  !=====================================
  subroutine io_print_timings()
  !=====================================
  !! Output timing information to stdout
  !=====================================

    implicit none

    integer :: i
    write(stdout,'(/1x,a)') '*===========================================================================*'
    write(stdout,'(1x,a)')  '|                             TIMING INFORMATION                            |'
    write(stdout,'(1x,a)')  '*===========================================================================*'
    write(stdout,'(1x,a)')  '|    Tag                                                Ncalls      Time (s)|'    
    write(stdout,'(1x,a)')  '|--------------------------------------------------:----------;--;----------|'    
    do i=1,nnames
      write(stdout,'(1x,"|",a50,":",i10,4x,f10.3,"|")') &
           clocks(i)%label,clocks(i)%ncalls,clocks(i)%ctime 
    enddo
    write(stdout,'(1x,a)')  '*---------------------------------------------------------------------------*'
    
    return

  end subroutine io_print_timings
  
  
  !========================================
  subroutine io_error ( error_msg )
  !========================================
  !! Abort the code giving an error message 
  !========================================

    implicit none
    character(len=*), intent(in) :: error_msg

    write(stdout,*)  'Exiting.......' 
    write(stdout, '(1x,a)') trim(error_msg)    
    close(stdout)    
    write(*, '(1x,a)') trim(error_msg)
    write(*,'(A)') "Error: examine the output/error file for details" 
    STOP
         
  end subroutine io_error
    
    
  !=======================================================
  !subroutine io_date(cdate, ctime)
  !=======================================================
  !                                                      
  !! Returns two strings containing the date and the time 
  !! in human-readable format. Uses a standard f90 call.
  !                                                    
  !=======================================================
  !  implicit none
  !  character (len=9), intent(out) :: cdate
  !  !! The date
  !  character (len=9), intent(out) :: ctime
  !  !! The time
  !
  !  character(len=3), dimension(12) :: months
  !  data months /'Jan','Feb','Mar','Apr','May','Jun',   &
  !       'Jul','Aug','Sep','Oct','Nov','Dec'/
  !  integer date_time(8)
  !  !
  !  call date_and_time(values=date_time)   !返回日期和时间，参考P.741
  !  !
  !  write (cdate,'(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
  !  write (ctime,'(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)
  !
  !end subroutine io_date
  

  function io_time()
  !===========================================================                                                        
  !! Returns elapsed CPU time in seconds since its first call.
  !! Uses standard f90 call                                                                                           
  !===========================================================
    use kinds, only : dp
    implicit none

    real(kind=dp) :: io_time

    ! t0 contains the time of the first call
    ! t1 contains the present time
    real(kind=dp) :: t0, t1
    logical :: first=.true.
    save first, t0
    !
    call cpu_time(t1)
    !
    if (first) then
      t0 = t1
      io_time = 0.0_dp
      first = .false.
    else
      io_time = t1 - t0
    endif
    return
  end function io_time


  function io_file_unit() !得到一个当前未使用的unit，用于打开文件
  !===========================================                                     
  !! Returns an unused unit number
  !! so we can later open a file on that unit.                                       
  !===========================================
  implicit none

    integer :: io_file_unit,unit_index
    logical :: file_open

    unit_index = 9
    file_open  = .true.
    do while ( file_open )
      unit_index = unit_index + 1
      inquire( unit=unit_index, OPENED = file_open ) !用于检查文件状态，参考P.536
    end do
    
    io_file_unit = unit_index

    return
  
  end function io_file_unit
  
  subroutine open_file(file_name,file_unit)
    implicit none
    
    character(len=*),intent(in) :: file_name
    integer,intent(in)          :: file_unit
    open(unit=file_unit, file=file_name,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem opening "'//trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif
  end subroutine open_file
  
  subroutine close_file(file_name,file_unit)
    implicit none
    
    integer,intent(in)          :: file_unit
    character(len=*),intent(in) :: file_name
    close(file_unit,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem close "'//trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif   
  end subroutine close_file
  
  function int_to_char(i)
  !-----------------------------------------------------------------------
  !
  ! ... converts an integer number of up to 12 figures
  ! ... into a left-justifed character variable 
  !  
  implicit none
  integer, intent(in) :: i
  character(len=12)   :: int_to_char
    CHARACTER :: c
    INTEGER   :: n, j, nc
    LOGICAL   :: neg
    !   
    nc = 6
    !
    IF( i < 0 ) then
       nc  = nc - 1
       n   = -i
       neg = .true.
    ELSE
       n   = i
       neg = .false.
    END IF
    !
    j = 1
    DO WHILE( j <= nc ) 
       int_to_char(j:j) = CHAR( MOD( n, 10 ) + ICHAR( '0' ) )
       n = n / 10
       IF( n == 0 ) EXIT
       j = j + 1
    END DO
    !
    IF( j <= nc ) THEN
       DO n = 1, j/2
          c = int_to_char( n : n )
          int_to_char( n : n ) = int_to_char( j-n+1 : j-n+1 )
          int_to_char( j-n+1 : j-n+1 ) = c
       END DO
       IF( j < nc ) int_to_char(j+1:nc) = ' '
    ELSE
       int_to_char(:) = '*'
    END IF
    !
    IF( neg ) THEN
       DO n = nc+1, 2, -1
          int_to_char(n:n) = int_to_char(n-1:n-1)
       END DO
       int_to_char(1:1) = '-'
    END IF
    !
    RETURN
    !
  END FUNCTION int_to_char  
  
  subroutine findkword(funit,kword)
    implicit none
      integer,intent(in):: funit
      character(len=*),intent(in) :: kword
      logical :: lfindkword 
      character(len=maxlen) :: ctmp
      
      lfindkword = .FAlSE.
      do while(.NOT. lfindkword)
        read(funit,*) ctmp
        if (trim(adjustl(ctmp))==trim(adjustl(kword))) then
          lfindkword = .TRUE.
          backspace(unit=funit)
        endif
      enddo  
  end subroutine findkword

  subroutine findkline(funit,kline,indexi,indexe)
    implicit none
      integer,intent(in):: funit
      character(len=*),intent(in) :: kline
      integer,intent(in)::indexi,indexe
      logical :: lfindkline
      character(len=maxlen) :: ctmp
      
      lfindkline = .FAlSE.
      do while(.NOT. lfindkline)
        read(funit,"(A)") ctmp
        if (ctmp(indexi:indexe)==kline) then
          lfindkline = .TRUE.
          backspace(unit=funit)
        endif
      enddo  
  end subroutine findkline
  
  
end module io
  
