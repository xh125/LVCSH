module hamiltonian
  use kinds,only : dp,dpc
	use types
  use parameters
  use elph2, only  : nbndfst,nk=>nktotf,gmnvkq,nq => nqtotf,wf,&
                     ibndmin,ibndmax
  use modes, only  : nmodes
  use readepw,only : E_nk

  use surfacecom,only     : lelecsh,lholesh,ieband_min,ieband_max,&
                            ihband_min,ihband_max
  
  implicit none
  
  integer :: neband,nhband,nefre,nhfre,nphfre
  real(kind=dp),allocatable :: H0_e(:,:),H_e(:,:),Enk_e(:,:)
  real(kind=dp),allocatable :: H0_e_nk(:,:,:,:),H_e_nk(:,:,:,:)
  real(kind=dp),allocatable :: H0_h(:,:),H_h(:,:),Enk_h(:,:)
  real(kind=dp),allocatable :: H0_h_nk(:,:,:,:),H_h_nk(:,:,:,:)

	real(kind=dp),allocatable :: H_e_eq(:,:),E_e_eq(:),P_e_eq(:,:) !for position of equilibrium
	real(kind=dp),allocatable :: H_h_eq(:,:),E_h_eq(:),P_h_eq(:,:) !for position of equilibrium
  
  real(kind=dp),allocatable :: gmnvkq_e(:,:,:,:,:)
  real(kind=dp),allocatable :: gmnvkq_h(:,:,:,:,:)
	
  logical,allocatable :: lgmnvkq_e(:,:,:,:,:)
  logical,allocatable :: lgmnvkq_h(:,:,:,:,:)	
	integer :: ngfre_e,ngfre_h
	type(gmnvkq_n0),allocatable :: gmnvkq_n0_e(:),gmnvkq_n0_h(:)
	
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable :: E_e(:),P_e(:,:),P_e_nk(:,:,:),&
                               E0_e(:),P0_e(:,:),P0_e_nk(:,:,:)
  real(kind=dp),allocatable :: E_h(:),P_h(:,:),P_h_nk(:,:,:),&
                               E0_h(:),P0_h(:,:),P0_h_nk(:,:,:)
  !哈密顿量的本征值与本征矢   
  integer :: ierr
  contains
  
  subroutine allocate_hamiltonian(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
    implicit none
    logical,intent(in) :: lelecsh,lholesh
    integer,intent(in) :: ieband_min,ieband_max,ihband_min,ihband_max
    if(lelecsh) then
      neband = ieband_max - ieband_min + 1
      nefre   = neband * nk
			allocate(Enk_e(neband,nk))
      allocate(gmnvkq_e(neband,neband,nmodes,nk,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating gmnvkq_e',1)
      gmnvkq_e = 0.0
      allocate(lgmnvkq_e(neband,neband,nmodes,nk,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating lgmnvkq_e',1)
      lgmnvkq_e = .false.			
      allocate(H0_e(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H0_e',1)    
      H0_e = 0.0
      allocate(H0_e_nk(neband,nk,neband,nk),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H0_e_nk',1)
      H0_e_nk = 0.0
      allocate(H_e(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H_e',1)
      H_e = 0.0
      allocate(H_e_nk(neband,nk,neband,nk),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H_e_nk',1)
      H_e_nk = 0.0
      allocate(E_e_eq(nefre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating E_e_eq',1)
      E_e_eq = 0.0
      allocate(P_e_eq(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating P_e_eq',1)
      P_e_eq = 0.0      
      allocate(H_e_eq(nefre,nefre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H_e_eq',1)
      H_e_eq = 0.0      			
			
    endif
    
    if(lholesh) then
      nhband = ihband_max - ihband_min + 1
      nhfre   = nhband * nk
			allocate(Enk_h(nhband,nk))
      allocate(gmnvkq_h(nhband,nhband,nmodes,nk,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating gmnvkq_h',1)
      gmnvkq_h = 0.0
      allocate(lgmnvkq_h(nhband,nhband,nmodes,nk,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating lgmnvkq_h',1)
      lgmnvkq_h = .false.            
      allocate(H0_h(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H0_h',1)    
      H0_h = 0.0
      allocate(H0_h_nk(nhband,nk,nhband,nk),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H0_h_nk',1)
      H0_h_nk = 0.0
      allocate(H_h(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H_h',1)
      H_h = 0.0
      allocate(H_h_nk(nhband,nk,nhband,nk),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating H_h_nk',1)
      H_h_nk = 0.0   
      allocate(E_h_eq(nhfre),stat=ierr)
      if(ierr /=0) call errore('Energy of hamiltonian for position of equilibrium ','Error allocating E_h_eq',1)
      E_h_eq = 0.0
      allocate(P_h_eq(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating P_h_eq',1)
      P_h_eq = 0.0      
      allocate(H_h_eq(nhfre,nhfre),stat=ierr)
      if(ierr /=0) call errore('hamiltonian for position of equilibrium ','Error allocating E_h_eq',1)
      H_h_eq = 0.0
			
    endif
    
  end subroutine allocate_hamiltonian
  
  subroutine set_H0_nk(nk,nband,Enk,H0_nk,iehband_min,gmnvkq_eh,lit,ngfre)
    use elph2, only : etf,gmnvkq
    implicit none
    integer , intent(in) :: nk
    integer , intent(in) :: nband
    integer , intent(in) :: iehband_min
    real(kind=dp), intent(out) :: Enk(nband,nk),H0_nk(nband,nk,nband,nk)
    real(kind=dp), intent(out) :: gmnvkq_eh(nband,nband,nmodes,nk,nq)
    real(kind=dp), intent(in) :: lit
		integer, intent(out) :: ngfre
		
    integer :: ik,iband,jband,nfre
    integer :: m,n,nu,iq
		
    nfre = nband*nk
    
    H0_nk = 0.0
    do ik=1,nk
      do iband=1,nband
				Enk(iband,ik) = etf(iband+iehband_min-1,ik)
        H0_nk(iband,ik,iband,ik)=etf(iband+iehband_min-1,ik)
      enddo
    enddo  
      
    iband = iehband_min
    jband = iehband_min+nband-1
    gmnvkq_eh = gmnvkq(iband:jband,iband:jband,:,:,:)    
		
		ngfre= 0
		do iq=1,nq
			do ik=1,nk
				do nu=1,nmodes
					do n=iband,jband
						do m=iband,jband
							if(ABS(gmnvkq(m,n,nu,ik,iq))>lit) ngfre = ngfre+1
						enddo
					enddo
				enddo
			enddo
		enddo
		
  end subroutine set_H0_nk
  

  !ref: 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !     (3.20) (6.4)
  !ref: 1 F. Giustino, Reviews of Modern Physics 89 (2017) 015003.
  !     (1)  
  subroutine set_H_nk(nband,nk,nmodes,nq,ph_Q,wf,gmnvkq_eh,H0_nk,H_nk)
    use epwcom,only  : kqmap
    implicit none
    integer , intent(in) :: nband,nk,nmodes,nq
    real(kind=dp),intent(in) :: ph_Q(nmodes,nq),wf(nmodes,nq)
    real(kind=dp),intent(in) :: H0_nk(nband,nk,nband,nk)
    real(kind=dp),intent(in) :: gmnvkq_eh(nband,nband,nmodes,nk,nq)
    real(kind=dp),intent(out):: H_nk(nband,nk,nband,nk)
    
    integer :: iq,nu,iband,jband
    integer :: ik,ikq    
    
    H_nk = H0_nk
    
    do iq=1,nq
      do ik=1,nk
        ikq=kqmap(ik,iq)
				do nu=1,nmodes					
          do iband=1,nband !|iband,ik>
            do jband=1,nband !|jband,ikq>
              H_nk(jband,ikq,iband,ik) = H_nk(jband,ikq,iband,ik)+&
              gmnvkq_eh(iband,jband,nu,ik,iq)*ph_Q(nu,iq)*sqrt(2.0*wf(nu,iq)/nq)
            enddo
          enddo
        enddo
      enddo
    enddo
        
  end subroutine set_H_nk

	subroutine get_gmnvkq_n0(lit,nk,nband,nq,nmode,gmnvkq,lgmnvkq,ngfre,g_n0)
		use epwcom,only  : kqmap
		use types
		implicit none
		integer ,intent(in) :: nk,nband,nq,nmode,ngfre
		real(kind=dp),intent(in) :: lit
		real(kind=dp),intent(in) :: gmnvkq(nband,nband,nmode,nk,nq)
		logical,intent(out) :: lgmnvkq(nband,nband,nmode,nk,nq)
		type(gmnvkq_n0),intent(out) :: g_n0(ngfre)
		
		integer :: ik,iband,jband,iq,imode,igfre
		
		igfre=0
		lgmnvkq = .false.
		do iq=1,nq
			do ik=1,nk
				do imode=1,nmode
					do iband=1,nband
						do jband=1,nband
							if(abs(gmnvkq(jband,iband,imode,ik,iq))>lit) then
								igfre=igfre+1
								g_n0(igfre)%m=jband
								g_n0(igfre)%n=iband
								g_n0(igfre)%v=imode
								g_n0(igfre)%ik=ik
								g_n0(igfre)%iq=iq
								g_n0(igfre)%ikq=kqmap(ik,iq)
								g_n0(igfre)%g= gmnvkq(jband,iband,imode,ik,iq)
								lgmnvkq(jband,iband,imode,ik,iq) = .true.
							endif
						enddo
					enddo
				enddo
			enddo
		enddo
		
		
		
	end subroutine get_gmnvkq_n0



  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!  
  subroutine calculate_eigen_energy_state(nfre,H,ee,pp)
    use f95_precision
    use lapack95
    implicit none
    integer ,intent(in) :: nfre
    real(kind=dp),intent(in) :: H(nfre,nfre)
    real(kind=dp),intent(out):: ee(nfre),pp(nfre,nfre) 
		
    pp = H
    call syev(pp,ee,'V','U')	
    !!USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !!On exit, hh array is overwritten
    !call heev(pp,ee,'V','U')    !P1143 MKL
    
  end subroutine calculate_eigen_energy_state
  
	subroutine resort_eigen_energy_stat(nfre,ee,pp,ee_eq,pp_eq)
		implicit none
		integer,intent(in) :: nfre
		real(kind=dp),intent(inout) :: ee(nfre),pp(nfre,nfre)
		real(kind=dp),intent(in) :: ee_eq(nfre),pp_eq(nfre,nfre)
		
		real(kind=dp),allocatable :: p_tmp(:),pdotp(:)
		real(kind=dp) :: e_tmp,flad
		
		integer :: ifre,jfre,cfre(1)
		real(kind=dp) :: maxsv
		integer :: nmax
		
		if(.not. allocated(p_tmp)) allocate(p_tmp(nfre))
		if(.not. allocated(pdotp)) allocate(pdotp(nfre))	
		
		do ifre =1 ,nfre
				pdotp = 0.0
				do jfre=ifre,nfre
					pdotp(jfre) = SUM(pp(:,jfre)*pp_eq(:,ifre))
				enddo
				cfre = Maxloc(ABS(pdotp))
				maxsv= Maxval(ABS(pdotp))
				nmax = 0
				do jfre=ifre,nfre
					if(ABS(pdotp(jfre)) == maxsv) then
						nmax=nmax+1
						if(pdotp(jfre) > pdotp(cfre(1))) cfre(1) = jfre
					endif
				enddo
				
				e_tmp= ee(ifre)
				p_tmp= pp(:,ifre)
				ee(ifre) = ee(cfre(1))
				pp(:,ifre) = pp(:,cfre(1))
				ee(cfre(1))= e_tmp
				pp(:,cfre(1)) = p_tmp
		enddo
		
	end subroutine
	
end module hamiltonian