module constants
  use kinds , only : dp,dpc
  !
  ! ... The constants needed everywhere
  !  
  implicit none
  save
  integer,parameter,public        :: maxlen = 256
  character(len=5),parameter      :: crash_file = 'CRASH'
  !
  !...Mathmatical constants
  !
  real(kind=dp),parameter,public  ::             &
    sqrt2 = dsqrt(2.0d0)                        ,&
    sqrt3 = dsqrt(3.0d0)                        ,&
    sqrt5 = dsqrt(5.0d0)                        ,&
    sqrt7 = dsqrt(7.0d0)
  real(kind=dp),parameter,public  :: pi  = 3.141592653589793238462643383279_dp 
  real(kind=dp),parameter,public  :: tpi = 2.0*pi
  real(kind=dp),parameter,public  :: fpi = 4.0*pi
  real(kind=dp),parameter,public  :: sqrtpi = 1.77245385090551602729_dp
  !real(kind=dp),parameter,public  :: sqrtpi = dsqrt(pi)
  real(kind=dp),parameter,public  :: sqrtpim1 = 1.0_dp/sqrtpi
  
  complex(kind=dpc),parameter,public   :: eye=(0.0d0,1.0d0)
  complex(kind=dpc), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
  !! i as a complex variable
  complex(kind=dpc), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  !! 0 as a complex variable
  complex(kind=dpc), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)
  !! 1 as a complex variable
  COMPLEX(KIND = DP), PARAMETER :: ci    = (0.0E0_DP, 1.0E0_DP)
  COMPLEX(KIND = DP), PARAMETER :: cone  = (1.0E0_DP, 0.0E0_DP)
  COMPLEX(KIND = DP), PARAMETER :: czero = (0.0E0_DP, 0.0E0_DP)  
	
	
  !~~ NUMERICAL CONVERGENCE CONSTANTS ~~!
  real(kind=dp), parameter, public    :: eps2  = 1.0e-2_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps5  = 1.0e-5_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps6  = 1.0e-6_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps7  = 1.0e-7_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps8  = 1.0e-8_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp
  !! numerical convergence constant  
  

  !~~ PHYSICAL CONSTANTS ~~!
  !
  ! Values of the fundamental constants taken from 
  ! http://physics.nist.gov/cuu/Constants/index.html
  ! ##### CODATA 2010 ##### !
  ! #warning "SCSH INFO: Using CODATA 2010 constant values"
    !au2amu 原子单位质量(me)转换为相对原子质量(C12/12)
    !the atomic mass of a carbon-12 atom is about 1.998467052 × 10−26 kg
    !one twelfth of the mass of an unbound neutral atom of carbon-12 
    !in its nuclear and electronic ground state and at rest has a value 
    !of 1.660539040(20)×10−27 kg,
    !electron mass 9.10938291(40)×10−31 kg  

  !
  ! ... Physical constants, SI (NIST 2018)
  !     http://physics.nist.gov/constants
  !
  
  real(kind=dp), parameter, public :: H_PLANCK_SI= 6.62607015E-34_DP      ! J*s  !J Hz-1
  real(kind=dp), parameter, public :: hbar_SI=1.054571726e-34_dp          ! J*s =H_PLANCK_SI/tpi
  !! hbar               ->  $$\hbar$$
  real(kind=dp), parameter, public :: k_B_SI=1.3806488e-23_dp             ! J/K
  real(kind=dp), parameter, public :: K_BOLTZMANN_SI = 1.380649E-23_DP    ! J * K^-1
  !! Boltzman Constant  ->  $$k_B$$   
  real(kind=dp), parameter, public :: Avo_con_SI=6.022140857e23_dp        !mol-1
  real(kind=dp), parameter, public :: elem_charge_SI=1.602176565e-19_dp   ! C
  !real(kind=dp), parameter, public :: sqrt_elem_charge_SI=elem_charge_SI**2
  !! elemental charge   ->  e**2
  real(kind=dp), parameter, public :: electronvolt_SI = 1.602176565e-19_dp! J
  real(kind=dp), parameter, public :: elec_mass_SI    = 9.10938291e-31_dp      ! kg
  real(kind=dp), parameter, public :: ELECTRONMASS_SI = 9.1093837015E-31_DP! Kg
  real(kind=dp), parameter, public :: Hartree_SI      = 4.3597447222071E-18_DP   ! J
  real(kind=dp), parameter, public :: Rydbebg_SI      = Hartree_SI/2.0_dp       ! J
  REAL(kind=dp), PARAMETER, public :: BOHR_RADIUS_SI  = 0.529177210903E-10_DP  ! m  
  real(kind=dp), parameter, public :: amu_mass_SI     = 1.660539040e-27_dp     ! kg
  REAL(kind=DP), PARAMETER, public :: AMU_SI          = 1.66053906660E-27_DP  ! Kg
  REAL(kind=DP), PARAMETER, public :: C_SI            = 2.99792458E+8_DP    ! m sec^-1  
  !! electron mass      ->  $$m_e$$

  !
  ! ... Physical constants, atomic units:
  ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
  ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)


  real(kind=dp),parameter :: e2 = 2.0  !the square of the electron charge
  !
  ! ... Unit conversion factors: energy and masses
  !
  REAL(DP), PARAMETER :: AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
  REAL(DP), PARAMETER :: RYTOEV           = AUTOEV / 2.0_DP
  real(dp), parameter :: Ry2J             = RYTOEV *electronvolt_SI
  REAL(DP), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI
  REAL(DP), PARAMETER :: AMU_RY           = AMU_AU / 2.0_DP
  REAL(KIND = DP), PARAMETER :: ryd2ev   = rytoev                     ! 13.6056981
  REAL(KIND = DP), PARAMETER :: ryd2mev  = ryd2ev * 1.0E3_DP          ! 13605.6981
  real(kind = dp), parameter :: ryd2V    = ryd2ev*sqrt(e2)
  real(kind = dp), parameter :: Ryd2V_m  = Ryd2V/BOHR_RADIUS_SI
  real(kind = dp), parameter :: K_B_Ryd  = k_B_SI/ Ry2J
  !! 
  ! ... Unit conversion factors: atomic unit of time, in s and ps
  !
  REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI 
  ! AU_SEC = T/tpi   (rad/s)^-1
  REAL(DP), PARAMETER :: AU_PS            = AU_SEC * 1.0E+12_DP
  !ev2thz
  real(dp), parameter :: ev_sec          = H_PLANCK_SI/tpi/(electronvolt_SI)
  real(dp), parameter :: ev_ps           = ev_sec * 1.0E+12_DP
  real(dp), parameter :: ev_ThZ          = 1.0/ev_ps/tpi
  real(dp), parameter :: mev_THZ         = (0.001*electronvolt_SI)/H_PLANCK_SI/1.0E+12_DP





  real(kind=dp), parameter, public :: bohr_magn_SI=927.400968e-26_dp      ! J/T
  !! Bohr magneton      ->  $$\mu_B$$
  real(kind=dp), parameter, public :: eps0_SI=8.854187817e-12_dp          !F / m
  !! Vacuum Dielectric Constant ->  $$\epsilon_0$$
  real(kind=dp), parameter, public :: speedlight_SI=299792458.0_dp        ! m / s
  !! Speed of light     ->  $$c$$
  real(kind=dp), parameter, public :: eV_au=3.674932379e-2_dp             ! (see table of Conv. Factors)  
  !! Electron Volt in atomic units
  real(kind=dp), parameter, public :: bohr_angstrom_internal=0.52917721092_dp
  !! Bohr to Anstrom Conversion factor
  real(kind=dp), parameter, public :: bohr = bohr_angstrom_internal
  !! 4*pi*eps0
  real(kind=dp), parameter, public :: fopieps0 = 4*pi*eps0_SI
  !real(kind=dp), parameter, public :: THz2womiga = 1.0e12_dp                     !f=  rad/s
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% SET PHYSICAL CONSTANTS                              %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/ATOMIC_UNIT       %!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/ATOMIC_MASS_UNIT  %!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/PHYSICAL_CONSTANT %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  real(kind=dp),parameter,public  ::             &
    au2cm = 2.194887656d5                       ,&
    au2ev = 2.7211d1                            ,&
    au2mev= 2.7211d4                            ,&
    au2j  = 4.35974417d-18                      
  real(kind=dp),parameter,public  ::             &
    au2s  = 2.418884326505d-17                  ,&
    au2fs = 2.418884326505d-2                   ,&
    au2ps = 2.418884326505d-5                  
  real(kind=dp),parameter,public  ::  au2amu= 5.485798701848d-4
  real(kind=dp),parameter,public  ::  au2k  = 3.1577464d5 
  real(kind=dp),parameter,public  ::  au2ang= 5.291772108d-1
  
  
  !!!!!! COMPATIBIILITY
  !
  REAL(DP), PARAMETER :: BOHR_RADIUS_CM = BOHR_RADIUS_SI * 100.0_DP
  REAL(DP), PARAMETER :: BOHR_RADIUS_ANGS = BOHR_RADIUS_CM * 1.0E8_DP
  REAL(DP), PARAMETER :: ANGSTROM_AU = 1.0_DP/BOHR_RADIUS_ANGS
  REAL(DP), PARAMETER :: AU_TERAHERTZ  = AU_PS
  REAL(DP), PARAMETER :: AU_TO_OHMCMM1 = 46000.0_DP ! (ohm cm)^-1
  REAL(DP), PARAMETER :: RY_TO_THZ = 1.0_DP / AU_TERAHERTZ / FPI
  REAL(DP), PARAMETER :: RY_TO_GHZ = RY_TO_THZ*1000.0_DP
  REAL(DP), PARAMETER :: RY_TO_CMM1= 1.E+10_DP * RY_TO_THZ / C_SI
  real(dp), parameter :: Ry_TO_fs  = 2.0 * au2fs
  !
  REAL(DP), PARAMETER :: AVOGADRO = 6.02214076D+23  
  
  
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: constants_version_str1 = "-> Using CODATA 2010 constant values"
  character(len=75), parameter, public :: constants_version_str2 = "   (http://physics.nist.gov/cuu/Constants/index.html)"

end module constants
  