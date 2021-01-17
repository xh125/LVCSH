module hamiltonian
  use kinds,only : dp,dpc
  use parameters
  use elph2, only  : nbndfst,nktotf,epcq,nqtotf,wf
  use modes, only  : nmodes
  use readepw,only : Enk
  implicit none
  real(kind=dp),allocatable     :: H0_e(:,:,:,:)
  real(kind=dp),allocatable     :: H0(:,:)
  real(kind=dp),allocatable     :: Hep(:,:,:),HH(:,:),H_tmp(:,:)
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable     :: E(:),P(:,:),E0(:),P0(:,:)
  !哈密顿量的本征值与本征矢   
  integer :: nefre,nphfre
  integer :: ierr
  contains
  
  subroutine set_H0()
    implicit none
    integer :: ik,iband,ibasis
    nefre = nbndfst * nktotf
    nphfre = nmodes * nqtotf
    
    allocate(H0(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('hamiltonian','Error allocating H0',1)    
    H0=0.0
    allocate(H0_e(nbndfst,nktotf,nbndfst,nktotf),stat=ierr)
    if(ierr /=0) call errore('hamiltonian','Error allocating H0_e',1)
    H0_e = 0.0
    !reshape( H0,(/ nbndfst,nktotf,nbndfst,nktotf /))
    allocate(HH(nefre,nefre),stat=ierr)
    if(ierr /=0) call errore('hamiltonian','Error allocating HH',1)        
    
    do ik=1,nktotf
      do iband=1,nbndfst
        H0_e(iband,ik,iband,ik)=Enk(iband,ik)
      enddo
    enddo
    H0=reshape(H0_e,(/ nefre,nefre /))
    
    
  end subroutine set_H0
  
  

  
end module hamiltonian