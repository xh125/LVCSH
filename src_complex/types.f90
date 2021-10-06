module types
	use kinds,only :dp
	implicit none
	!declare type of gmnvkq_n0 with gmnvkq>0.0
	type :: gmnvkq_n0
		integer :: m
		integer :: n
		integer :: v
		integer :: ik
		integer :: iq
		integer :: ikq
		real(kind=dp) :: g
	end type gmnvkq_n0	
	
end module types