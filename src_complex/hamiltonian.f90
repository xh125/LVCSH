module hamiltonian
  use kinds,only : dp,dpc
	use io, only : msg,io_error,stdout
	use types
  use parameters
  use elph2, only  : nbndfst,nk=>nktotf,gmnvkq,nq => nqtotf,wf
  use modes, only  : nmodes
	use memory_report,only : MB,GB,complex_size, real_size,int_size,ram,print_memory
  use readepw,only : E_nk
	use constants,only : ryd2mev,czero,cone
  use surfacecom,only     : lelecsh,lholesh,ieband_min,ieband_max,&
                            ihband_min,ihband_max
  
  implicit none
  
  integer :: neband,nhband,nefre,nhfre,nphfre
	complex(kind=dp),allocatable :: H0_e(:,:),H_e(:,:),H0_e_nk(:,:,:,:),H_e_nk(:,:,:,:)
  complex(kind=dp),allocatable :: H0_h(:,:),H_h(:,:),H0_h_nk(:,:,:,:),H_h_nk(:,:,:,:)
	real(kind=dp),allocatable    :: Enk_e(:,:),Enk_h(:,:)

  
  complex(kind=dp),allocatable :: gmnvkq_e(:,:,:,:,:)
  complex(kind=dp),allocatable :: gmnvkq_h(:,:,:,:,:)
	
	
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable    :: E_e(:),E0_e(:)
	complex(kind=dp),allocatable :: P_e(:,:),P_e_nk(:,:,:),&
                                  P0_e(:,:),P0_e_nk(:,:,:)
  real(kind=dp),allocatable    :: E_h(:),E0_h(:)
	complex(kind=dp),allocatable :: P_h(:,:),P_h_nk(:,:,:),&
                                  P0_h(:,:),P0_h_nk(:,:,:)
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
			allocate(Enk_e(neband,nk),stat=ierr,errmsg=msg)
			if(ierr /= 0) call io_error(msg)
      allocate(gmnvkq_e(neband,neband,nmodes,nk,nq),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating gmnvkq_e',1)
				call io_error(msg)
			endif
			gmnvkq_e = 0.0
			ram = complex_size*neband*neband*nmodes*nk*nq
			call print_memory("gmnvkq_e",ram)
      allocate(H0_e(nefre,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H0_e',1)    
				call io_error(msg)
			endif
			H0_e = 0.0
			ram = complex_size*nefre*nefre
			call print_memory("H0_e",ram)
      allocate(H0_e_nk(neband,nk,neband,nk),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H0_e_nk',1)
				call io_error(msg)
			endif
			H0_e_nk = 0.0
			call print_memory("H0_e_nk",ram)
      allocate(H_e(nefre,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H_e',1)
				call io_error(msg)
			endif
			H_e = 0.0
			call print_memory("H_e",ram)
      allocate(H_e_nk(neband,nk,neband,nk),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H_e_nk',1)
				call io_error(msg)
			endif
			H_e_nk = 0.0
			call print_memory("H_e_nk",ram)
			
    endif
    
    if(lholesh) then
      nhband = ihband_max - ihband_min + 1
      nhfre   = nhband * nk
			allocate(Enk_h(nhband,nk),stat=ierr,errmsg=msg)
			if(ierr /= 0) call io_error(msg)
      allocate(gmnvkq_h(nhband,nhband,nmodes,nk,nq),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating gmnvkq_h',1)
				call io_error(msg)
		  endif
			gmnvkq_h = 0.0
			ram = complex_size*nhband*nhband*nmodes*nk*nq
			call print_memory("gmnvkq_h",ram)
      allocate(H0_h(nhfre,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H0_h',1)    
				call io_error(msg)
			endif
			H0_h = 0.0
			ram = complex_size*nhfre*nhfre
			call print_memory("H0_h",ram)			
      allocate(H0_h_nk(nhband,nk,nhband,nk),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H0_h_nk',1)
				call io_error(msg)
			endif
			H0_h_nk = 0.0
			call print_memory("H0_h_nk",ram)	
      allocate(H_h(nhfre,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H_h',1)
				call io_error(msg)
			endif
			H_h = 0.0
			ram = complex_size*nhfre*nhfre
			call print_memory("H_h",ram)
      allocate(H_h_nk(nhband,nk,nhband,nk),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H_h_nk',1)
				call io_error(msg)
			endif
			H_h_nk = 0.0
			call print_memory("H_h_nk",ram)				

    endif
    
  end subroutine allocate_hamiltonian
  
  subroutine set_H0_nk(nk,ibandmin,ibandmax,Enk,H0_nk,gmnvkq_eh)
    use elph2, only : etf,epmatq
    implicit none
    integer , intent(in) :: nk
    integer , intent(in) :: ibandmin,ibandmax
    real(kind=dp), intent(out) :: Enk(ibandmax-ibandmin+1,nk)
		complex(kind=dp), intent(out) :: H0_nk(ibandmax-ibandmin+1,nk,ibandmax-ibandmin+1,nk)
    complex(kind=dp), intent(out) :: gmnvkq_eh(ibandmax-ibandmin+1,ibandmax-ibandmin+1,nmodes,nk,nq)
		
		integer :: nband
    integer :: ik,iband
    integer :: m,n,nu,iq
		
		nband = ibandmax-ibandmin+1
    
    H0_nk = czero
    do ik=1,nk
      do iband=1,nband
				Enk(iband,ik) = etf(iband+ibandmin-1,ik)
        H0_nk(iband,ik,iband,ik)=etf(iband+ibandmin-1,ik)*cone
      enddo
    enddo  
      
    gmnvkq_eh = epmatq(ibandmin:ibandmax,ibandmin:ibandmax,:,:,:)    
		
  end subroutine set_H0_nk
  

  !ref: 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !     (3.20) (6.4)
  !ref: 1 F. Giustino, Reviews of Modern Physics 89 (2017) 015003.
  !     (1)  
  subroutine set_H_nk(nband,nk,nmodes,nq,ph_Q,gmnvkq_eh,H0_nk,H_nk)
    use epwcom,only  : kqmap
    implicit none
    integer , intent(in) :: nband,nk,nmodes,nq
    real(kind=dp),intent(in) :: ph_Q(nmodes,nq)
    complex(kind=dp),intent(in) :: H0_nk(nband,nk,nband,nk)
    complex(kind=dp),intent(in) :: gmnvkq_eh(nband,nband,nmodes,nk,nq)
    complex(kind=dp),intent(out):: H_nk(nband,nk,nband,nk)
    
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
              gmnvkq_eh(iband,jband,nu,ik,iq)*ph_Q(nu,iq)
            enddo
          enddo
        enddo
      enddo
    enddo
        
  end subroutine set_H_nk


  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!  
  subroutine calculate_eigen_energy_state(nfre,H,ee,pp)
    use f95_precision
    use lapack95
    implicit none
    integer ,intent(in) :: nfre
    complex(kind=dp),intent(in)  :: H(nfre,nfre)
    real(kind=dp),intent(out)    :: ee(nfre)
		complex(kind=dp),intent(out) :: pp(nfre,nfre) 
		
    pp = H
		call heevd(pp,ee,'V','U')
    !call syev(pp,ee,'V','U')	
    !!USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !!On exit, hh array is overwritten
    !call heev(pp,ee,'V','U')    !P1143 MKL
    
  end subroutine calculate_eigen_energy_state
  
end module hamiltonian