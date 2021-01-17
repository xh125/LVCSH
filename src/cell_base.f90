module cell_base
  use kinds,only      : dp
  use constants,only  : pi,maxlen
  use io,only         : stdout
  
  implicit none
  save
  !  ibrav: index of the bravais lattice  (see latgen.f90)
  !----------------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111)
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !-----------------------------------------------------------------------------  
  integer :: ibrav
  !  celldm: old-style parameters of the simulation cell (se latgen.f90)
  real(kind=dp) :: celldm(6) = (/ 0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp /)
  !  traditional crystallographic cell parameters (alpha=cosbc and so on)  
  real(kind=dp) :: a, b, c, cosab, cosac, cosbc
  real(kind=dp) :: a1(3),a2(3),a3(3)
  ! format of input cell parameters:
  ! 'alat','bohr','angstrom'
  character(len=maxlen) :: cell_units
  !  alat: lattice parameter - often used to scale quantities, or
  !  in combination to other parameters/constants to define new units  
  real(kind=dp) :: alat = 0.0_dp
  ! omega: volume of the simulation cell
  real(kind=dp) :: omega= 0.0_dp
  ! tpiba: 2 PI/alat, tpiba2=tpiba^2
  real(kind=dp) :: tpiba = 0.0_dp , tpiba2 = 0.0_dp
  !  direct and reciprocal lattice primitive vectors
  !  at(:,i) are the lattice vectors of the simulation cell, a_i,
  !          in alat units: a_i(:) = at(:,i)/alat
  !  bg(:,i) are the reciprocal lattice vectors, b_i,
  !          in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba  
  real(kind=dp) :: at(3,3) = reshape( (/ 0.0_dp /), (/ 3,3 /), (/ 0.0_dp /) )
  real(kind=dp) :: bg(3,3) = reshape( (/ 0.0_dp /), (/ 3,3 /), (/ 0.0_dp /) )
  
  ! atomic coordinate referred to the crystal axes
  !real(kind=dp),allocatable :: xau(:,:)
  !real(kind=dp),allocatable :: tau(:,:)
  !integer :: nat, ntyp
  
  
  contains
  
end module cell_base    