module hamiltonian
  use kinds,only : dp,dpc
	use io,   only : msg,io_error,stdout
	use types
  use parameters
  use elph2, only  : nbndfst,epmatq,wf,wf_sub,nqtotf
  use modes, only  : nmodes
  use readepw,only : E_nk
	use constants,only : ryd2mev,cone,czero
  use surfacecom,only: lelecsh,lholesh,ieband_min,ieband_max,&
                       ihband_min,ihband_max,nq_sub
	use memory_report,only : MB,GB,complex_size, real_size,int_size,ram,print_memory  
  
  implicit none
  
  integer :: neband,nhband,nefre,nhfre,nphfre
  
  real(kind=dp),   allocatable :: Enk_e(:,:),Enk_h(:,:)
  complex(kind=dp),allocatable :: H0_e(:,:),H_e(:,:)
  complex(kind=dp),allocatable :: H0_e_nk(:,:,:,:),H_e_nk(:,:,:,:)
  complex(kind=dp),allocatable :: H0_h(:,:),H_h(:,:)
  complex(kind=dp),allocatable :: H0_h_nk(:,:,:,:),H_h_nk(:,:,:,:)
  
  complex(kind=dp),allocatable :: gmnvkq_e(:,:,:,:,:)
  complex(kind=dp),allocatable :: gmnvkq_h(:,:,:,:,:)
	
  type(gmnvkq_n0),allocatable :: gijk_e(:),gijk_h(:)
  
  integer :: gnzfree_e,gnzfree_h
  
	integer :: ngfre_e,ngfre_h
	
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable    :: E_e(:),E0_e(:),E_h(:),E0_h(:)
  complex(kind=dp),allocatable :: P_e(:,:),P_e_nk(:,:,:),P0_e(:,:),P0_e_nk(:,:,:),&
                                  P_h(:,:),P_h_nk(:,:,:),P0_h(:,:),P0_h_nk(:,:,:)
  !哈密顿量的本征值与本征矢   
  integer :: ierr
  contains
  
  subroutine allocate_hamiltonian(lelecsh,lholesh,nk_sub,ieband_min,ieband_max,ihband_min,ihband_max)
    implicit none
    logical,intent(in) :: lelecsh,lholesh
    integer,intent(in) :: nk_sub,ieband_min,ieband_max,ihband_min,ihband_max
    if(lelecsh) then
      neband = ieband_max - ieband_min + 1
      nefre   = neband * nk_sub
			allocate(Enk_e(neband,nk_sub),stat=ierr,errmsg=msg)
			if(ierr /= 0) call io_error(msg)
      allocate(gmnvkq_e(neband,neband,nmodes,nk_sub,nq_sub),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating gmnvkq_e',1)
				call io_error(msg)
			endif
			gmnvkq_e = 0.0
			ram = complex_size*neband*neband*nmodes*nk_sub*nq_sub
			call print_memory("gmnvkq_e",ram)
      allocate(H0_e(nefre,nefre),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H0_e',1)    
				call io_error(msg)
			endif
			H0_e = 0.0
			ram = complex_size*nefre*nefre
			call print_memory("H0_e",ram)
      allocate(H0_e_nk(neband,nk_sub,neband,nk_sub),stat=ierr,errmsg=msg)
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
      allocate(H_e_nk(neband,nk_sub,neband,nk_sub),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H_e_nk',1)
				call io_error(msg)
			endif
			H_e_nk = 0.0
			call print_memory("H_e_nk",ram)
			
    endif
    
    if(lholesh) then
      nhband = ihband_max - ihband_min + 1
      nhfre   = nhband * nk_sub
			allocate(Enk_h(nhband,nk_sub),stat=ierr,errmsg=msg)
			if(ierr /= 0) call io_error(msg)
      allocate(gmnvkq_h(nhband,nhband,nmodes,nk_sub,nq_sub),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating gmnvkq_h',1)
				call io_error(msg)
		  endif
			gmnvkq_h = 0.0
			ram = complex_size*nhband*nhband*nmodes*nk_sub*nq_sub
			call print_memory("gmnvkq_h",ram)
      allocate(H0_h(nhfre,nhfre),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H0_h',1)    
				call io_error(msg)
			endif
			H0_h = 0.0
			ram = complex_size*nhfre*nhfre
			call print_memory("H0_h",ram)			
      allocate(H0_h_nk(nhband,nk_sub,nhband,nk_sub),stat=ierr,errmsg=msg)
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
      allocate(H_h_nk(nhband,nk_sub,nhband,nk_sub),stat=ierr,errmsg=msg)
      if(ierr /=0) then
				call errore('hamiltonian','Error allocating H_h_nk',1)
				call io_error(msg)
			endif
			H_h_nk = 0.0
			call print_memory("H_h_nk",ram)				

    endif
    
  end subroutine allocate_hamiltonian
  
  subroutine set_H0_nk(nk_sub,indexk,nq_sub,indexq,ibandmin,ibandmax,Enk,H0_nk,gmnvkq_eh,lit_g,gnzfree)
    use elph2, only : etf_sub,gmnvkq
    use epwcom,only : kqmap_sub
    implicit none
    integer , intent(in) :: nk_sub,nq_sub
    integer , intent(in) :: indexk(nk_sub),indexq(nq_sub),ibandmin,ibandmax
    real(kind=dp), intent(out) :: Enk(ibandmax-ibandmin+1,nk_sub)
		complex(kind=dp), intent(out) :: H0_nk(ibandmax-ibandmin+1,nk_sub,ibandmax-ibandmin+1,nk_sub)
    complex(kind=dp), intent(out) :: gmnvkq_eh(ibandmax-ibandmin+1,ibandmax-ibandmin+1,nmodes,nk_sub,nq_sub)
		integer , intent(out) :: gnzfree
    real(kind=dp) , intent(in) :: lit_g
    
		integer :: nband
    integer :: ik,iband,jband
    integer :: m,n,nu,iq,ikq
		
		nband = ibandmax-ibandmin+1
    
    H0_nk = 0.0
    do ik=1,nk_sub
      do iband=1,nband
				Enk(iband,ik) = etf_sub(iband+ibandmin-1,ik)
        H0_nk(iband,ik,iband,ik)=etf_sub(iband+ibandmin-1,ik)*cone
      enddo
    enddo  
    
    do iq=1,nq_sub
      do ik=1,nk_sub
        gmnvkq_eh(:,:,:,ik,iq) = epmatq(ibandmin:ibandmax,ibandmin:ibandmax,:,indexk(ik),indexq(iq))    
      enddo
    enddo
    
    gnzfree = 0
    do iq=1,nq_sub
      do ik=1,nk_sub
        ikq = kqmap_sub(ik,iq)
        if(ikq /= 0) then
          do nu=1,nmodes
            do iband=1,nband
              do jband=1,nband
                if(ABS(gmnvkq_eh(jband,iband,nu,ik,iq)) > sqrt(2.0*wf_sub(nu,iq)/nqtotf)*lit_g) gnzfree = gnzfree +1
              enddo
            enddo
          enddo
        endif
      enddo
    enddo
    
    
  end subroutine set_H0_nk
  
  subroutine set_gijk(nband,nmodes,nk_sub,nq_sub,gmnvkq,lit_g,gnzfree,gijk)
    use epwcom,only : kqmap_sub
    implicit none
    integer, intent(in) :: nband,nmodes,nk_sub,nq_sub,gnzfree
    complex(kind=dpc),intent(in) :: gmnvkq(nband,nband,nmodes,nk_sub,nq_sub)
    real(kind=dp) , intent(in)   :: lit_g
    type(gmnvkq_n0),intent(out)  :: gijk(gnzfree)
    
    integer :: iq,ik,ikq,nu,iband,jband,ink,iqv,imkq
    integer :: ig
    
    ig = 0
    do iq=1,nq_sub
      do ik=1,nk_sub
        ikq = kqmap_sub(ik,iq)
        if(ikq /= 0) then
          do nu=1,nmodes
            do iband=1,nband
              do jband=1,nband
                if(ABS(gmnvkq(jband,iband,nu,ik,iq)) > sqrt(2.0*wf_sub(nu,iq)/nqtotf)*lit_g) then
                  ig = ig +1
                  gijk(ig)% iqv = (iq-1)*nmodes + nu
                  gijk(ig)% ink = (ik-1)*nband  + iband
                  gijk(ig)% imkq= (ikq-1)*nband + jband
                  gijk(ig)% g   = gmnvkq(jband,iband,nu,ik,iq)
                endif  
              enddo
            enddo
          enddo
        endif
      enddo
    enddo
    
    
    
  end subroutine

  !ref: 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !     (3.20) (6.4)
  !ref: 1 F. Giustino, Reviews of Modern Physics 89 (2017) 015003.
  !     (1)  
  subroutine set_H_nk(nband,nk,nmodes,nq,ph_Q,gmnvkq_eh,H0_nk,H_nk)
    use epwcom,only  : kqmap_sub
    implicit none
    integer , intent(in) :: nband,nk,nmodes,nq
    complex(kind=dpc),intent(in) :: ph_Q(nmodes,nq)
    complex(kind=dpc),intent(in) :: H0_nk(nband,nk,nband,nk)
    complex(kind=dpc),intent(in) :: gmnvkq_eh(nband,nband,nmodes,nk,nq)
    complex(kind=dpc),intent(out):: H_nk(nband,nk,nband,nk)
    
    integer :: iq,nu,iband,jband
    integer :: ik,ikq    
    
    H_nk = H0_nk
    
    do iq=1,nq
      do ik=1,nk
        ikq=kqmap_sub(ik,iq)
        if(ikq/=0) then
          do nu=1,nmodes					
            do iband=1,nband !|iband,ik>
              do jband=1,nband !|jband,ikq>
                H_nk(jband,ikq,iband,ik) = H_nk(jband,ikq,iband,ik)+&
                gmnvkq_eh(iband,jband,nu,ik,iq)*ph_Q(nu,iq)
              enddo
            enddo
          enddo
        endif
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
    complex(kind=dpc),intent(in)  :: H(nfre,nfre)
    real(kind=dp),intent(out)     :: ee(nfre)
    complex(kind=dpc),intent(out) :: pp(nfre,nfre) 
		
    pp = H
    call heevd(pp,ee,'V','U')
    !call syev(pp,ee,'V','U')	
    !!USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !!On exit, hh array is overwritten
    !call heev(pp,ee,'V','U')    !P1143 MKL
    
  end subroutine calculate_eigen_energy_state
  
end module hamiltonian