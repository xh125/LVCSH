module kpoints
  use kinds,only : dp
  implicit none
    !integer :: nkstot = 0, nk1 = 0, nk2 = 0, nk3 = 0, k1 = 0, k2 = 0, k3 = 0
    !real(kind=dp), allocatable :: xk(:,:) , wk(:)
    
    integer :: npk
    !! max number of k-points
    integer  :: k1
    !! the offset from the origin, direction 1
    integer  :: k2
    !! the offset from the origin, direction 2
    integer  :: k3
    !! the offset from the origin, direction 3
    integer  :: nk1
    !! the special-point grid, direction 1
    integer  :: nk2
    !! the special-point grid, direction 2
    integer  :: nk3
    !! the special-point grid, direction 3
    integer  :: t_rev(48)
    !! time reversal flag, for noncolinear magnetism
    logical  :: time_reversal
    !! if .TRUE. the system has time reversal symmetry
    logical  :: skip_equivalence
    !! if .TRUE. skip check of k-points equivalence
    integer :: nirkps
    !! number of irreducible k points
    
    
    real(kind=dp),allocatable :: xk(:,:),wk(:)
    !! coordinates of k points
    !! weight of k points    
    real(kind=dp),allocatable :: xkg(:,:),wkk(:)
    integer :: nktot
    
    !!!!!
    ! nkr=nk1*nk2*nk3
    ! allocate(xkg(3,nkr),wkk(nkr))
    
  
end module kpoints