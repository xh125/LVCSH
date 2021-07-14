module memory_report
	use kinds,only : dp
	use io,only    : stdout
	implicit none
	
	integer , parameter :: MB=1024*1024
	integer , parameter :: GB=1024*MB
	real(kind=dp), parameter :: complex_size=16_dp, real_size=8_dp, int_size=4_dp
	real(kind=dp) :: ram
	
	contains
	
	subroutine print_memory(ram_name,ram_)
		implicit none
		character(len=*),intent(in) :: ram_name
		real(kind=dp),intent(in) :: ram_
		
		write(stdout,1013) ram_name, ram_/MB


1013 format(/5X,'Dynamical RAM for ', A19, ': ', F10.2, ' MB')		
		
	end subroutine print_memory
	
end module memory_report