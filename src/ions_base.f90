module ions_base
  use kinds, only: dp
  
  implicit none
    !     nsp       = number of species
    !     na(is)    = number of atoms of species is
    !     nax       = max number of atoms of a given species
    !     nat       = total number of atoms of all species
    !     iat        i-th atom
    
    integer,parameter :: ntypx = 10
    INTEGER :: nsp     = 0
    INTEGER :: na(ntypx) = 0
    INTEGER :: nax     = 0
    INTEGER :: nat     = 0
    !! total number of atoms of all species    
    integer :: ntyp    = 0 
    
    integer :: iat

    !     zv(is)    = (pseudo-)atomic charge
    !     amass(is) = mass of ions, in atomic mass units
    !     rcmax(is) = Ewald radius (for ion-ion interactions)

    REAL(DP) :: zv(ntypx)    = 0.0_DP
    REAL(DP) :: amass(ntypx) = 0.0_DP
    REAL(DP) :: rcmax(ntypx) = 0.0_DP

    !     ityp( i ) = the type of i-th atom in stdin
    !     iatm( i ) = the name of the i-th atom  in stdin
    !     atm( j )  = name of the type of the j-th atomic specie
    !     tau( 1:3, i ) = position of the i-th atom

    INTEGER,  ALLOCATABLE         :: ityp(:)
    character(len=3), allocatable :: iatm(:)
    real(kind=dp),allocatable     :: iamass(:)    ! iamass(nat)
    REAL(DP), ALLOCATABLE         :: tau(:,:)     !  initial positions read from stdin (in bohr)
    real(dp), allocatable         :: xau(:,:)     !  
    REAL(DP), ALLOCATABLE         :: vel(:,:)     !  initial velocities read from stdin (in bohr)
    CHARACTER(LEN=3)              :: atm( ntypx )
    CHARACTER(LEN=80)             :: tau_format   ! format of input atomic positions:
                                          ! 'alat','crystal','bohr','angstrom'  
  
  
end module 