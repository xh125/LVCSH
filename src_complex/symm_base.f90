module symm_base
  use kinds , only : dp
  implicit none
  
  INTEGER :: s(3,3,48)
  !! symmetry matrices, in crystal axis
  real(kind=dp) :: sr(3,3,48)
  !! symmetry matrices , in cartesian axis
  character(len=45) :: sname(48)
  !! name of the symmetries  
  INTEGER :: t_rev(48) = 0 
  !! time reversal flag, for noncolinear magnetism  
  LOGICAL :: time_reversal = .TRUE.
  !! if .TRUE. the system has time reversal symmetry
  INTEGER :: invs(48)
  !! index of inverse operation: S^{-1}_i=S(invs(i))
  integer :: nrot
  !! number of bravais lattice symmetries
  integer :: spacegroup
  !! space group index, as read from input
  integer :: nsym,isym
  !! total number of crystal symmetries
  
  character(len=11) :: group_name
  !! name of the group
  
  
  
  
  
end module 