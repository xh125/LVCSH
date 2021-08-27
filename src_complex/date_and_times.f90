module date_and_times
  implicit none
    
    character(len=3), dimension(12) :: months
    data months /'Jan','Feb','Mar','Apr','May','Jun',   &
         'Jul','Aug','Sep','Oct','Nov','Dec'/
    integer :: date_time(8)
  
  contains 
  
  subroutine get_date_and_time(cdate,ctime)
  !=======================================================
  !                                                      
  !! Returns two strings containing the date and the time 
  !! in human-readable format. Uses a standard f90 call.
  !                                                    
  !=======================================================
    implicit none
    character (len=9), intent(out) :: cdate
    !! The date
    character (len=9), intent(out) :: ctime
    !! The time

    !
    call date_and_time(values=date_time)   !返回日期和时间，参考P.741
    !
    write (cdate,'(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
    write (ctime,'(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)
  
  end subroutine get_date_and_time
  
end module date_and_times
