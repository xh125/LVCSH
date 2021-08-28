module sortting
	use kinds,only : dp	,dpc
  use constants,only : cone,czero
	use io,only : stdout
	implicit none
	
	contains
	
	subroutine resort_eigen_energy_stat(nfre,ee,pp,ee_0,pp_0,P_lim)
		implicit none
		
		!=============================================================================================
		!min-cost or min-sum assignment problem
		!ref : 1 S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
    !ref : 2 G. Carpaneto, S. Martello, and P. Toth, Annals of Operations Research 13 (1988) 191.
		!ref : 3 E. Fabiano, T. W. Keal, and W. Thiel, Chem. Phys. 349 (2008) 334.
		!ref : 4 Temen, S. & Akimov, A. V. A Simple Solution to Trivial Crossings: A Stochastic State Tracking Approach. 
		!			!J Phys Chem Lett 12, 850-860, doi:10.1021/acs.jpclett.0c03428 (2021).
		!==============================================================================================
		
		integer,intent(in) :: nfre
		real(kind=dp),intent(inout) :: ee(nfre)
    complex(kind=dpc),intent(inout) :: pp(nfre,nfre)
		! energy and states at time t+δt
		real(kind=dp),intent(in) 		 :: ee_0(nfre)
    complex(kind=dpc),intent(in) :: pp_0(nfre,nfre)
		! energy and states at time t 
		real(kind=dp),intent(in)		:: P_lim
		! threshold to find the mixxing states
		
		!========================================================================
		! |pp_0(ifre)> = Σ_jfre |pp(jfre)><pp(jfre)|pp_0(ifre)>
		! |pp_0(ifre)> = Σ_jfre |pp(jfre)>a_ji(jfre,ifre)
		! |pp_0(ifre)> = Σ_jfre a_ji(jfre,ifre)|pp(jfre)>
		! a_ji(jfre,ifre) = <pp(jfre)|pp_0(ifre)> = S_ij(jfre,ifre)
		! a_ji(jfre,ifre) = SUM(CONJG(pp(:,jfre))*pp_0(:,ifre)) = S_ij(jfre,ifre)
		!========================================================================
		
		complex(kind=dpc),allocatable :: S_ij(:,:)
		! S_ij(t,t+δt)=  <pp(jfre)|pp_0(ifre)>   
		! equation(4) of S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		real(kind=dp),allocatable :: P_ij(:,:)
		! P_ij = |S_ij(t,t+δt)|^2
		! equation(5) of S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		real(kind=dp),allocatable :: P_i(:),P_ii(:)
		!P_i(:) = P_ij(:,ifre)
		!P_j(:) = P_ij(jfre,:)
		!P_ii(ifre) = P_ij(ifre,ifre)
		
		integer :: ifre,jfre,ifre_,jfre_
		integer :: cfre(1),maxfre
		real(kind=dp) :: maxoverlap
		integer :: nmix
		real(kind=dp) :: SUM_Pik
		! SUM_Pik = SUM(P_i(:))
		
		integer,allocatable :: i_mix(:),j_mix(:)
		
		allocate(S_ij(nfre,nfre))
		allocate(P_ij(nfre,nfre))
		allocate(P_i(nfre),P_ii(nfre))
		S_ij = czero
		P_ij = 0.0
		P_ii = 1.0
		P_i  = 0.0
		
		allocate(i_mix(nfre),j_mix(nfre))
		i_mix = 0
		j_mix = 0
		
		! ref : 1 S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		! calculate the overlap matrix
		do ifre =1 ,nfre 	 ! |pp_0(ifre)>
			do jfre = 1,nfre ! |pp(jfre)>
				S_ij(jfre,ifre) = SUM(CONJG(pp(:,jfre))*pp_0(:,ifre))
				! if states is a complex function
				!S_ij(jfre,ifre) = SUM(pp(:,jfre)*pp_0(:,ifre))
			enddo
		enddo
		
		 P_ij = S_ij*CONJG(S_ij)
		! equation(5) of S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		!P_ij = S_ij**2
		do ifre=1,nfre
			SUM_Pik = SUM(P_ij(:,ifre))
			!SUM_Pik = 1.0
			P_ij(:,ifre) = P_ij(:,ifre)/SUM_Pik
			P_ii(ifre)   = P_ij(ifre,ifre)
		enddo
		deallocate(S_ij)
		
		
		! not mixxing for different states
		! "trivial" unavoided crossing, "non-interacting-state" crossings		
		do ifre = 1, nfre ! |pp_0(ifre)>
			
			P_i(:) = P_ij(:,ifre)			
			
			cfre  		 = Maxloc(P_i)
			maxfre		 = cfre(1)
			maxoverlap = Maxval(P_i)
			
			! "trivial" unavoided crossing, "non-interacting-state" crossings
			if(maxfre /= ifre .and. maxoverlap > P_lim ) then ! exchange
				call exchange_states(nfre,ifre,maxfre,ee,pp,P_ij)
				!exchange the ee(:), pp(:,:) and P_ij
			endif	
			P_ii(ifre) = P_ij(ifre,ifre)
		enddo		

		!!for mixxing states, lets |pp(ifre)> must be with a largest or second larger population of |pp_0(ifre)>
		!do ifre = 1,nfre ! |pp_0(ifre)>
		!	if(P_ij(ifre,ifre)<= P_lim) then
		!		P_i(:) = P_ij(:,ifre)
		!		
		!		mix: do
		!			cfre  		 = Maxloc(P_i)
		!			maxfre		 = cfre(1)
		!			
		!			if(maxfre>=ifre .and. P_ij(maxfre,maxfre) < P_lim ) then
		!				call exchange_states(nfre,ifre,maxfre,ee,pp,P_ij)
		!				EXIT mix
		!			else
		!				P_i(maxfre) = 0.0
		!			endif
		!		enddo mix
		!		
		!	endif
		!	P_ii(ifre) = P_ij(ifre,ifre)
		!enddo
		!
		!
		!!mixxing for different states
		!do ifre = 1, nfre ! |pp_0(ifre)>
		!	if(P_ii(ifre)<= P_lim) then
		!		P_i(:) = P_ij(:,ifre)
		!		nmix   = 0
		!		maxoverlap = 0.0
		!		
		!		num_mix: do
		!			cfre  		 = Maxloc(P_i)
		!			maxfre		 = cfre(1)
		!			maxoverlap = Maxval(P_i)+maxoverlap
		!			P_i(maxfre)= 0.0
		!			nmix = nmix + 1
		!			j_mix(nmix) = maxfre
		!			if(maxoverlap > P_lim) EXIT
		!		enddo num_mix
		!		
		!		do jfre_=1,nmix
		!			jfre = j_mix(jfre_)
		!			if(P_ij(ifre,ifre)*P_ij(jfre,jfre) < P_ij(jfre,ifre)*P_ij(ifre,jfre) ) then
		!				call exchange_states(nfre,ifre,jfre,ee,pp,P_ij)
		!			endif
		!		enddo
		!		
		!	endif
		!	
		!enddo
		
		
		
		
		deallocate(P_ij,P_ii,P_i,i_mix,j_mix)
		
	end subroutine resort_eigen_energy_stat

	subroutine exchange_states(nfre,ifre,maxfre,ee,pp,P_ij)
		implicit none
		integer,intent(in) :: nfre,ifre,maxfre
		real(kind=dp),intent(inout) :: ee(nfre),P_ij(nfre,nfre)
    complex(kind=dpc),intent(inout) :: pp(nfre,nfre)

		! used as temp to resort the energy and states and overlap
		real(kind=dp) :: ee_tmp
		real(kind=dpc),allocatable :: pp_tmp(:)
		real(kind=dpc),allocatable :: Pj_tmp(:)

		allocate(pp_tmp(nfre))
		allocate(Pj_tmp(nfre))
		pp_tmp = 0.0
		Pj_tmp = 0.0	
		
		ee_tmp  = ee(ifre)
		pp_tmp  = pp(:,ifre)
		Pj_tmp  = P_ij(ifre,:)
		
		ee(ifre) = ee(maxfre)
		pp(:,ifre) = pp(:,maxfre)
		P_ij(ifre,:) = P_ij(maxfre,:)
		
		ee(maxfre)= ee_tmp
		pp(:,maxfre) = pp_tmp
		P_ij(maxfre,:)= Pj_tmp		
		
		deallocate(pp_tmp,Pj_tmp)
		
	end subroutine exchange_states

end module 