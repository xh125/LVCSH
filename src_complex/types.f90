module types
	use kinds,only :dp,dpc
	implicit none
	!declare type of gmnvkq_n0 with gmnvkq>0.0
	type :: gmnvkq_n0
		integer :: iqv
		integer :: ink
		integer :: imkq
		complex(kind=dpc) :: g
	end type gmnvkq_n0	
	
end module types