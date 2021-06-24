module hamiltonian
  use kinds,only : dp,dpc
  use parameters
  use elph2, only  : nbndfst,nk=>nktotf,gmnvkq,nq => nqtotf,wf,&
                     ibndmin,ibndmax
  use modes, only  : nmodes
  use readepw,only : E_nk

  use surfacecom,only     : lelecsh,lholesh,ieband_min,ieband_max,&
                            ihband_min,ihband_max
  
  implicit none
  
  integer :: neband,nhband,nefre,nhfre,nphfre
  real(kind=dp),allocatable :: H0_e(:,:),H_e(:,:)
  real(kind=dp),allocatable :: H0_e_nk(:,:,:,:),H_e_nk(:,:,:,:)
  real(kind=dp),allocatable :: H0_h(:,:),H_h(:,:)
  real(kind=dp),allocatable :: H0_h_nk(:,:,:,:),H_h_nk(:,:,:,:)
  
  real(kind=dp),allocatable :: gmnvkq_e(:,:,:,:,:)
  real(kind=dp),allocatable :: gmnvkq_h(:,:,:,:,:)
  
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable :: E_e(:),P_e(:,:),P_e_nk(:,:,:),&
                               E0_e(:),P0_e(:,:),P0_e_nk(:,:,:)
  real(kind=dp),allocatable :: E_h(:),P_h(:,:),P_h_nk(:,:,:),&
                               E0_h(:),P0_h(:,:),P0_h_nk(:,:,:)
  real(kind=dp),allocatable :: E_h_(:)
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
      allocate(gmnvkq_e(neband,neband,nmodes,nk,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating gmnvkq_e',1)
      gmnvkq_e = 0.0
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
      
    endif
    
    if(lholesh) then
      nhband = ihband_max - ihband_min + 1
      nhfre   = nhband * nk
      allocate(gmnvkq_h(nhband,nhband,nmodes,nk,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating gmnvkq_h',1)
      gmnvkq_h = 0.0      
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
    endif
    
  end subroutine allocate_hamiltonian
  
  subroutine set_H0_nk(nk,nband,H0_nk,iehband_min,gmnvkq_eh)
    use elph2, only : etf,gmnvkq
    implicit none
    integer , intent(in) :: nk
    integer , intent(in) :: nband
    integer , intent(in) :: iehband_min
    real(kind=dp), intent(out) :: H0_nk(nband,nk,nband,nk)
    real(kind=dp), intent(out) :: gmnvkq_eh(nband,nband,nmodes,nk,nq)
    
    integer :: ik,iband,jband,nfre
    
    nfre = nband*nk
    
    H0_nk = 0.0
    do ik=1,nk
      do iband=1,nband
        H0_nk(iband,ik,iband,ik)=etf(iband+iehband_min-1,ik)
      enddo
    enddo  
      
    iband = iehband_min
    jband = iehband_min+nband-1
    gmnvkq_eh = gmnvkq(iband:jband,iband:jband,:,:,:)    
    
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
              ph_Q(nu,iq)*sqrt(2.0*wf(nu,iq)/nq)*gmnvkq_eh(iband,jband,nu,ik,iq)
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
    real(kind=dp),intent(in) :: H(nfre,nfre)
    real(kind=dp),intent(out):: ee(nfre),pp(nfre,nfre) 
		
    pp = H
    call syev(pp,ee,'V','U')	
    !!USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !!On exit, hh array is overwritten
    !call heev(pp,ee,'V','U')    !P1143 MKL
    
  end subroutine calculate_eigen_energy_state
  
end module hamiltonian