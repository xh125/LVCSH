module sortting
	use kinds,only : dp	
	implicit none
	
	contains
	
	subroutine resort_eigen_energy_stat(nfre,ee,pp,ee_0,pp_0,S_lim)
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
		real(kind=dp),intent(inout) :: ee(nfre),pp(nfre,nfre)
		! energy and states at time t+δt
		real(kind=dp),intent(in) 		:: ee_0(nfre),pp_0(nfre,nfre)
		! energy and states at time t 
		real(kind=dp),intent(in)		:: S_lim
		! threshold to find the mixxing states
		
		!=======================================================================
		! |pp_0(ifre)> = Σ_jfre |pp(jfre)><pp(jfre)|pp_0(ifre)>
		! |pp_0(ifre)> = Σ_jfre |pp(jfre)>a_ji(jfre,ifre)
		! |pp_0(ifre)> = Σ_jfre a_ji(jfre,ifre)|pp(jfre)>
		! a_ji(jfre,ifre) = <pp(jfre)|pp_0(ifre)> = S_ij(jfre,ifre)
		! a_ji(jfre,ifre) = SUM(CONJG(pp(:,jfre))*pp_0(:,ifre))
		!=======================================================================
		
		real(kind=dp),allocatable :: S_ij(:,:),S_ii(:)
		! S_ij(t,t+δt)=  <pp(jfre)|pp_0(ifre)>   
		! equation(4) of S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		real(kind=dp),allocatable :: P_ij(:,:)
		! P_ij = |S_ij(t,t+δt)|^2
		! equation(5) of S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		real(kind=dp),allocatable :: P_i(:)
		!P_i(:) = P_ij(:,ifre)
		
		! used as temp to resort the energy and states and overlap
		real(kind=dp) :: ee_tmp
		real(kind=dp),allocatable :: pp_tmp(:)
		real(kind=dp),allocatable :: S_tmp(:)
		!real(kind=dp) :: flad,c_tmp
		
		integer :: ifre,jfre
		integer :: cfre(1),maxfre
		real(kind=dp) :: maxoverlap!,sumc2
		!integer :: nmax
		real(kind=dp) :: SUM_Sik
		
		allocate(S_ij(nfre,nfre),S_ii(nfre))
		allocate(P_i(nfre))
		S_ij = 0.0
		S_ii = 1.0
		P_i  = 0.0
		
		allocate(pp_tmp(nfre))
		allocate(S_tmp(nfre))
		pp_tmp = 0.0
		S_tmp  = 0.0	
		
		! ref : 1 S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		! calculate the overlap matrix
		do ifre =1 ,nfre 	 ! |pp_0(ifre)>
			do jfre = 1,nfre ! |pp(jfre)>
				!S_ij(jfre,ifre) = SUM(CONJG(pp(:,jfre))*pp_0(:,ifre))
				! if states is a complex function
				S_ij(jfre,ifre) = SUM(pp(:,jfre)*pp_0(:,ifre))
			enddo
		enddo
		
		! S_ij = S_ij*CONJG(S_ij)
		! equation(5) of S. Fernandez-Alberti et al., J Chem Phys 137 (2012) 014512.
		S_ij = S_ij**2
		do ifre=1,nfre
			SUM_Sik = SUM(S_ij(:,ifre))
			S_ij(:,ifre) = S_ij(:,ifre)/SUM_Sik
			S_ii(ifre) = S_ij(ifre,ifre)
		enddo
		
		! not mixxing for different states
		! "trivial" unavoided crossing, "non-interacting-state" crossings		
		do ifre = 1, nfre ! |pp_0(ifre)>
			
			P_i(:) = S_ij(:,ifre)			
			
			cfre  		 = Maxloc(P_i)
			maxfre		 = cfre(1)
			maxoverlap = Maxval(P_i)
			
			! "trivial" unavoided crossing, "non-interacting-state" crossings
			if(maxfre /= ifre .and. maxoverlap > S_lim ) then ! exchange
				!exchange the ee(:), pp(:,:) and S_ij
				ee_tmp = ee(ifre)
				pp_tmp = pp(:,ifre)
				S_tmp  = S_ij(ifre,:)
				ee(ifre) = ee(maxfre)
				pp(:,ifre) = pp(:,maxfre)
				S_ij(ifre,:) = S_ij(maxfre,:)
				ee(maxfre)= ee_tmp
				pp(:,maxfre) = pp_tmp
				S_ij(maxfre,:)= S_tmp
			endif	
			S_ii(ifre) = S_ij(ifre,ifre)
		enddo		
		
		!do ifre=1, nfre ! pp_0(ifre)
		!	
		!	pij2 = 0.0
		!	do jfre=1,nfre !pp(jfre)
		!		pij2(jfre) = SUM(pp_0(:,ifre)*pp(:,jfre))
		!		!pp_0(ifre) = SUM_jfre(pij2(jfre)*pp(jfre))
		!	enddo
		!	cfre = Maxloc(pij2**2)
		!	maxfre = cfre(1)
		!	maxsv= Maxval(pij2**2)
		!	nmax = 0
		!	do jfre=1,nfre
		!		if(ABS(pij2(jfre)) == maxsv) then
		!			nmax=nmax+1
		!			if(pij2(jfre) > pij2(maxfre)) maxfre = jfre
		!		endif
		!	enddo		
		!	
		!	e_tmp= ee(ifre)
		!	p_tmp= pp(:,ifre)
		!	ee(ifre) = ee(maxfre)
		!	pp(:,ifre) = pp(:,maxfre)
		!	ee(maxfre)= e_tmp
		!	pp(:,maxfre) = p_tmp
		!	if(nmax>1) then
		!		write(stdout,"(A,I8,A,I8,A)")  "pp_0(",ifre,"), project to ",nmax ,"pp(jfre) as same maxval."
		!	endif
		!enddo
		
		deallocate(S_ij,S_ii,P_i,pp_tmp,S_tmp)
		
	end subroutine

end module 