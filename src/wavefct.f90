module wavefct
!
!! The variables needed to compute the band structure
!
use kinds,only : dp
implicit none

  integer :: nband
  real(kind=dp),allocatable :: et(:,:)
  !! eigenvalues of the hamiltonian  
  real(kind=dp),allocatable :: wg(:,:)
  !! the weight of each k point and band
  
  
end module wavefct