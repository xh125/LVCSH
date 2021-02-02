module kpoint_grid
  use kinds, only : dp
  implicit none
  integer :: nrot
  integer :: npk
  integer :: k1, k2, k3
  integer :: nk1, nk2, nk3
  integer :: t_rev(48)
  integer :: s(3,3,48)
  logical :: time_reversal
  logical :: skip_equivalence
  real(kind=dp) :: bg(3,3)
  !! bg(:,i) are the reciprocal lattice vectors, b_i,
  !! in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba    
  integer :: nks
  !! number of k points
  real(kind=dp), allocatable :: xk(:,:)
  !! coordinates of k points
  real(kind=dp), allocatable :: wk(:)
  !! weight of k points
  
  contains 
  
  
  
end module 
