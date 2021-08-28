subroutine errore( calling_routine, message, ierr )
 !
 implicit none
 character(len=*), intent(in) :: calling_routine, message
 integer,          intent(in) :: ierr
 !
 ! ... the error message is written un the "*" unit
 !
 WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
 WRITE( UNIT = *, &
        FMT = '(5X,"from ",A," : error #",I10)' ) TRIM(calling_routine), ierr
 WRITE( UNIT = *, FMT = '(5X,A)' ) message
 WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
 !
 WRITE( *, '("     stopping ...")' )
 !
 STOP 2
 !
 RETURN
 ! 
end subroutine errore