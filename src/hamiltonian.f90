module hamiltonian
  use kinds,only : dp,dpc
  use parameters
  use elph2, only  : nbndfst,nk=>nktotf,epcq,nq => nqtotf,wf,&
                     ibndmin,ibndmax
  use modes, only  : nmodes
  use readepw,only : E_nk
  use epwcom,only  : kqmap
  use surfacecom,only     : lelecsh,lholesh,ieband_min,ieband_max,&
                            ihband_min,ihband_max
  implicit none
  real(kind=dp),allocatable :: H0_e(:,:),H_e(:,:)
  real(kind=dp),allocatable :: H0_e_nk(:,:,:,:),H_e_nk(:,:,:,:)
  real(kind=dp),allocatable :: H0_h(:,:),H_h(:,:)
  real(kind=dp),allocatable :: H0_h_nk(:,:,:,:),H_h_nk(:,:,:,:)
  
  real(kind=dp),allocatable     :: H0_nk(:,:,:,:),H_nk(:,:,:,:)
  real(kind=dp),allocatable     :: H0(:,:),H(:,:)
  
  real(kind=dp),allocatable :: epcq_e(:,:,:,:,:)
  real(kind=dp),allocatable :: epcq_h(:,:,:,:,:)
  
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable :: E_e(:),P_e(:,:),P_e_nk(:,:,:),&
                               E0_e(:),P0_e(:,:),P0_e_nk(:,:,:)
  real(kind=dp),allocatable :: E_h(:),P_h(:,:),P_h_nk(:,:,:),&
                               E0_h(:),P0_h(:,:),P0_h_nk(:,:,:)
  real(kind=dp),allocatable :: E(:),P(:,:),E0(:),P0(:,:)!,P_nk(:,:,:)
  !哈密顿量的本征值与本征矢   
  integer :: neband,nhband,nefre,nhfre,nphfre
  integer :: ierr
  contains
  
  subroutine set_H0_nk_eh(nk,nband,H0_nk,iehband_min,epcq_eh)
    use elph2, only : etf,epcq
    implicit none
    integer , intent(in) :: nk
    integer , intent(in) :: nband
    integer , intent(in) :: iehband_min
    real(kind=dp), intent(out) :: H0_nk(nband,nk,nband,nk)
    real(kind=dp), intent(out) :: epcq_eh(nband,nband,nk,nmodes,nq)
    
    integer :: ik,iband,jband,nfre
    
    nfre = nband*nk
    
    H0_nk = 0.0
    do ik=1,nk
      do iband=1,nband
        H0_nk(iband,ik,iband,ik)=etf(iband+iehband_min-1,2*ik-1) !E_nk(iband,ik)
      enddo
    enddo  
      
    iband = iehband_min-ibndmin+1
    jband = iehband_min-ibndmin+nband
    epcq_eh = epcq(iband:jband,iband:jband,:,:,:)    
    
    
  end subroutine set_H0_nk_eh
  
  subroutine set_H0_nk()
    use elph2, only : etf
    implicit none
    integer :: ik,iband,jband
    if(lelecsh) then
      neband = ieband_max - ieband_min + 1
      nefre   = neband * nk
      allocate(epcq_e(neband,neband,nk,nmodes,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating epcq_e',1)
      epcq_e = 0.0
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
      
      call set_H0_nk_eh(nk,neband,H0_e_nk,ieband_min,epcq_e)
      !do ik=1,nk
      !  do iband=1,neband
      !    H0_e_nk(iband,ik,iband,ik)=etf(iband+ieband_min-1,2*ik-1) !E_nk(iband,ik)
      !  enddo
      !enddo
      !H0_e = reshape(H0_e_nk,(/ nefre,nefre /))
      
      !iband = ieband_min-ibndmin+1
      !jband = ieband_max-ibndmin+1
      !epcq_e = epcq(iband:jband,iband:jband,:,:,:)
      
    endif
    
    if(lholesh) then
      nhband = ihband_max - ihband_min + 1
      nhfre   = nhband * nk
      allocate(epcq_h(nhband,nhband,nk,nmodes,nq),stat=ierr)
      if(ierr /=0) call errore('hamiltonian','Error allocating epcq_h',1)
      epcq_h = 0.0      
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
      
      call set_H0_nk_eh(nk,nhband,H0_h_nk,ihband_min,epcq_h)
      !do ik=1,nk
      !  do iband=1,nhband
      !    H0_e_nk(iband,ik,iband,ik)=etf(iband+ihband_min-1,2*ik-1) !E_nk(iband,ik)
      !  enddo
      !enddo
      !H0_h = reshape(H0_h_nk,(/ nhfre,nhfre /)) 
      !
      !iband  = ihband_min-ibndmin+1
      !jband  = ihband_max-ibndmin+1
      !epcq_h = epcq(iband:jband,iband:jband,:,:,:)
      
    endif
 
  end subroutine set_H0_nk
  

  !ref: 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !     (3.20) (6.4)
  !ref: 1 F. Giustino, Reviews of Modern Physics 89 (2017) 015003.
  !     (1)  
  subroutine set_H_nk(ph_Q)
    use kinds,only: dp
    
    implicit none
    real(kind=dp),intent(in) :: ph_Q(nmodes,nq)
    
    
    if(lelecsh) then
      call set_H_nk_eh(neband,nk,nmodes,nq,ph_Q,wf,epcq_e,H0_e_nk,H_e_nk)
      H_e = reshape(H_e_nk,(/ neband*nk,neband*nk /))
    endif
    
    if(lholesh) then
      call set_H_nk_eh(nhband,nk,nmodes,nq,ph_Q,wf,epcq_h,H0_h_nk,H_h_nk)
      H_h = reshape(H_h_nk,(/ nhband*nk,nhband*nk /))
    endif
      
    
    
    !H_nk = H0_nk
    !
    !do iq=1,nq
    !  do nu=1,nmodes
    !    do ik=1,nk
    !      ikq=kqmap(ik,iq)
    !      do iband=1,nband !|iband,ik>
    !        do jband=1,nband !|jband,ikq>
    !          H_nk(jband,ikq,iband,ik) = H_nk(jband,ikq,iband,ik)+&
    !          ph_Q(nu,iq)*sqrt(2.0*wf(nu,iq)/nq)*epcq(iband,jband,ik,nu,iq)
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !enddo
    !
    !H = reshape(H_nk,(/ nband*nk,nband*nk /))
    
  end subroutine set_H_nk

  subroutine set_H_nk_eh(nband,nk,nmodes,nq,ph_Q,wf,epcq_eh,H0_nk,H_nk)
    implicit none
    integer , intent(in) :: nband,nk,nmodes,nq
    real(kind=dp),intent(in) :: ph_Q(nmodes,nq),wf(nmodes,nq)
    real(kind=dp),intent(in) :: H0_nk(nband,nk,nband,nk)
    real(kind=dp),intent(in) :: epcq_eh(nband,nband,nk,nmodes,nq)
    real(kind=dp),intent(out):: H_nk(nband,nk,nband,nk)
    
    integer :: iq,nu,iband,jband
    integer :: ik,ikq    
    
    H_nk = H0_nk
    
    do iq=1,nq
      do nu=1,nmodes
        do ik=1,nk
          ikq=kqmap(ik,iq)
          do iband=1,nband !|iband,ik>
            do jband=1,nband !|jband,ikq>
              H_nk(jband,ikq,iband,ik) = H_nk(jband,ikq,iband,ik)+&
              ph_Q(nu,iq)*sqrt(2.0*wf(nu,iq)/nq)*epcq_eh(iband,jband,ik,nu,iq)
            enddo
          enddo
        enddo
      enddo
    enddo
    
    !H = reshape(H_nk,(/ nband*nk,nband*nk /))    
  end subroutine set_H_nk_eh


  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!  
  subroutine calculate_eigen_energy_state(nfre,H,ee,pp)
    !use kinds,only : dp
    !use mkl95_precision
    !use mkl95_lapack
    use f95_precision
    use lapack95
    implicit none
    integer ,intent(in) :: nfre
    real(kind=dp),intent(in) :: H(nfre,nfre)
    real(kind=dp),intent(out):: ee(nfre),pp(nfre,nfre) 
    
    !real(kind=dp):: pp(nband*nk,nband*nk)
    pp = H
    ! pp=reshape(H_nk,(/ nband*nk,nband*nk/))
    
    call syev(pp,ee,'V','U')
    !!USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !!On exit, hh array is overwritten
    !call heev(pp,ee,'V','U')    !P1143 MKL
    
    !p_nk = reshape(pp,(/ nband,nk,nband*nk/))
    
  end subroutine calculate_eigen_energy_state
  
end module hamiltonian