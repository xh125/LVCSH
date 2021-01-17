module klist
  use kinds,      only : dp
  use parameters, only : npk
  !  integer,parameter :: npk = 40000 ! max number of k-points in pw.x calculation
  implicit none
  character(len=32) :: smearing
  !! smearing type
  real(kind=dp),allocatable :: xk(:,:) !xk(3,npk)
  !! coordinates of k points  
  real(kind=dp),allocatable :: wk(:) !wk(npk)
  !! weight of k points  
  real(kind=dp) :: xqq(3)
  !! coordinates of q point (used in the ACFDT part)  
  REAL(DP) :: degauss
  !! smearing parameter
  REAL(DP) :: nelec
  !! number of electrons
  REAL(DP) :: nelup=0.0_dp
  !! number of spin-up electrons (if two_fermi_energies=t)
  REAL(DP) :: neldw=0.0_dp
  !! number of spin-dw electrons (if two_fermi_energies=t)
  REAL(DP) :: tot_magnetization
  !! nelup-neldw >= 0 (negative value means unspecified)
  REAL(DP) :: tot_charge
  !! total charge
  REAL(DP) :: qnorm= 0.0_dp
  !! |q|, used in phonon+US calculations only
  INTEGER, ALLOCATABLE :: igk_k(:,:)
  !! index of G corresponding to a given index of k+G
  INTEGER, ALLOCATABLE :: ngk(:)
  !! number of plane waves for each k point
  !
  INTEGER :: nks
  !! number of k points in this pool
  INTEGER :: nkstot
  !! total number of k points in scf calculation
  INTEGER :: ngauss
  !! type of smearing technique
  LOGICAL :: lgauss
  !! if .TRUE.: use gaussian broadening
  LOGICAL :: ltetra
  !! if .TRUE.: use tetrahedra
  LOGICAL :: lxkcry=.FALSE.
  !! if .TRUE.:k-pnts in cryst. basis accepted in input
  LOGICAL :: two_fermi_energies
  !! if .TRUE.: nelup and neldw set ef_up and ef_dw separately
  !  
  
  contains
  
end module klist

MODULE lsda_mod
  !
  !! It contains the variables needed for the LSDA calculation.
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx, npk
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: lsda
  !! true if lsda is active
  REAL(DP) :: magtot
  !! total magnetization
  REAL(DP) :: absmag
  !! total absolute magnetization
  REAL(DP) :: starting_magnetization(ntypx)
  !! the magnetization used to start with
  INTEGER :: nspin
  !! number of spin polarization: 2 if lsda, 1 other
  INTEGER :: current_spin
  !! spin of the current kpoint
  INTEGER :: isk(npk)
  !! for each k-point: 1=spin up, 2=spin down
  !
END MODULE lsda_mod