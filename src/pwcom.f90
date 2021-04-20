module klist
  !
  !! This module contains the variables related to the k-points.
  !
  use kinds,      only : dp
  use parameters, only : npk
  !  integer,parameter :: npk = 40000 ! max number of k-points in pw.x calculation
  implicit none
  character(len=32) :: smearing
  !! smearing type
  real(kind=dp),allocatable :: xk(:,:) !xk(3,npk)
  !! coordinates of k points in units 2pi/alat
  real(kind=dp),allocatable :: xkg(:,:) !xkg(3,npk)
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

module gvec
  use kinds,only : dp
  implicit none
  real(kind=dp) :: ecutwfc = 0.0
  real(kind=dp) :: ecutrho = 0.0
  logical :: dft_is_hybrid
  real(kind=dp) :: ecutfock
  logical :: lscf
  
  !   values for costant cut-off computations

  REAL(DP) :: ecfixed=0.0_DP     ! value of the constant cut-off
  REAL(DP) :: qcutz = 0.0_DP     ! height of the penalty function (above ecfix)
  REAL(DP) :: q2sigma=0.0_DP     ! spread of the penalty function around ecfix
  ! augmented cut-off for k-point calculation                                   

  ! ... G vectors with |G|^2 < 4*ecutwfc, cut-off for wavefunctions
  ! ... ("smooth" grid). Gamma tricks and units as for the "dense" grid
  !
  INTEGER :: ngms = 0  ! local  number of smooth vectors (on this processor)
  INTEGER :: ngms_g=0  ! global number of smooth vectors (summed on procs) 
                       ! in serial execution this is equal to ngms
  INTEGER :: ngm  = 0  ! local  number of G vectors (on this processor)
                       ! with gamma tricks, only vectors in G>
  INTEGER :: ngm_g= 0  ! global number of G vectors (summed on all procs)
                       ! in serial execution, ngm_g = ngm
  INTEGER :: ngl = 0   ! number of G-vector shells
  INTEGER :: ngmx = 0  ! local number of G vectors, maximum across all procs  
  INTEGER :: ngsx = 0  ! local number of smooth vectors, max across procs

  REAL(DP) :: ecuts = 0.0_DP   ! energy cut-off = 4*ecutwfc
  REAL(DP) :: gcutms= 0.0_DP   ! ecuts/(2 pi/a)^2, cut-off for |G|^2

  REAL(DP) :: dual = 0.0_DP    ! ecutrho=dual*ecutwfc
  LOGICAL  :: doublegrid = .FALSE. ! true if smooth and dense grid differ
                                   ! doublegrid = (dual > 4)
  
end module gvec

module funct
  use kinds,only : dp
  implicit none
  character(len=25) :: dft = 'not set'
  !
  ! PRIVATE variables defining the DFT functional
  !
  real(kind=dp) :: exx_fraction
  integer :: iexch, icorr, igcx, igcc, inlc, imeta, imetac
end module funct

module pw_control_flags
  use kinds,only : dp
  implicit none
  !=--------------------------------------------------------------------------=!
  !
  ! ... this module contains all basic variables that controls the
  ! ... pw.x execution flow
  !----------------------------------------------
  !
  INTEGER :: nstepe   = 1
                            !  parameters to control how many electronic steps
                            !  between ions move  
  !
  ! ... pw self-consistency
  !
  INTEGER, PUBLIC :: &
    ngm0,             &! used in mix_rho
    niter,            &! the maximum number of iteration
    nmix,             &! the number of iteration kept in the history
    imix               ! the type of mixing (0=plain,1=TF,2=local-TF)
  character(len=9) :: mixing_style 

  INTEGER,  PUBLIC :: &
    n_scf_steps        ! number of scf iterations to reach convergence
  REAL(DP), PUBLIC :: &
    mixing_beta,      &! the mixing parameter
    tr2,              &! the convergence threshold for potential
    scf_error=0.0      ! actual convergence reached

  LOGICAL, PUBLIC :: &
    conv_elec          ! if .TRUE. electron convergence has been reached
  ! next 3 variables used for EXX calculations
  LOGICAL, PUBLIC :: &
    adapt_thr       ! if .TRUE. an adaptive convergence threshold is used
                       ! for the scf cycle in an EXX calculation.
  REAL(DP), PUBLIC  :: &
    tr2_init,         &! initial value of tr2 for adaptive thresholds
    tr2_multi          ! the dexx multiplier for adaptive thresholds
                       ! tr2 = tr2_multi * dexx after each V_exx update 
  LOGICAL, PUBLIC :: scf_must_converge
  
  !
  ! ... Several variables controlling the run ( used mainly in PW calculations )
  !
  ! ... logical flags controlling the execution
  !
  LOGICAL, PUBLIC :: &
    lscf    =.FALSE., &! if .TRUE. the calc. is selfconsistent
    lbfgs   =.FALSE., &! if .TRUE. the calc. is a relaxation based on BFGS
    lmd     =.FALSE., &! if .TRUE. the calc. is a dynamics
    lwf     =.FALSE., &! if .TRUE. the calc. is with wannier functions
    !=================================================================
    !exx_wf related 
    lwfnscf =.FALSE., &
    lwfpbe0nscf=.FALSE.,&
    !=================================================================
    lbands  =.FALSE., &! if .TRUE. the calc. is band structure
    lconstrain=.FALSE.,&! if .TRUE. the calc. is constraint
    llondon =.FALSE., & ! if .TRUE. compute Grimme D2 dispersion corrections
    ldftd3 =.FALSE., & ! if .TRUE. compute Grimme D3 dispersion corrections
    ts_vdw  =.FALSE., & ! as above for Tkatchenko-Scheffler disp.corrections
    lxdm    =.FALSE., & ! if .TRUE. compute XDM dispersion corrections
    lensemb =.FALSE., &! if .TRUE. compute ensemble energies
    restart =.FALSE.   ! if .TRUE. restart from results of a preceding run
  
  !
  ! ... ionic dynamics
  !
  INTEGER, PUBLIC :: &
    nstep = 1,       &! number of ionic steps
    istep = 0          ! current ionic step
  LOGICAL, PUBLIC :: &
    conv_ions          ! if .TRUE. ionic convergence has been reached
  REAL(DP), PUBLIC  :: &
    upscale            ! maximum reduction of convergence threshold  

  ! ... Real-sapce algorithms
  !
  LOGICAL,          PUBLIC :: tqr=.FALSE. ! if true the Q are in real space

  !LOGICAL,          PUBLIC :: real_space=.false. ! beta functions in real space
  !
  ! ... Augmetation charge and beta smoothing
  !
  LOGICAL,          PUBLIC :: tq_smoothing=.FALSE. ! if true the Q are smoothed 
  LOGICAL,          PUBLIC :: tbeta_smoothing=.FALSE. ! if true the betas are smoothed 
  ! ... External Forces on Ions
  !
  LOGICAL,          PUBLIC :: textfor = .FALSE.


  LOGICAL,          PUBLIC :: treinit_gvecs = .FALSE.
  
  ! ... printout control
  !
  INTEGER, PUBLIC :: &
    io_level = 1       ! variable controlling the amount of I/O to file
  INTEGER, PUBLIC :: & ! variable controlling the amount of I/O to output
    iverbosity = 0     ! -1 minimal, 0 low, 1 medium, 2 high, 3 debug

    
end module pw_control_flags

MODULE vlocal
  !
  !! The variables needed for the local potential in reciprocal space.
  !
  USE kinds,       ONLY : DP
  USE parameters,  ONLY : ntypx
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: strf(:,:)
  !! the structure factor
  REAL(DP), ALLOCATABLE :: vloc(:,:)
  !! the local potential for each atom type
  REAL(DP) :: starting_charge(ntypx)
  !! the atomic charge used to start with
  !
END MODULE vlocal

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

module noncolin_module
  use kinds,only : dp
  implicit none
  LOGICAL :: &
    noncolin, &           !  true if noncollinear magnetism is allowed
    lsign=.FALSE.         !  if true use the sign feature to calculate
                            !  rhoup and rhodw
end module noncolin_module

MODULE spin_orb
  !
  !! Variables needed for calculations with spin-orbit
  !
  USE kinds,       ONLY : DP
  !USE upf_params,  ONLY : lmaxx, lqmax
  !! FIXME: rot_ylm could be dynamically allocated
  !
  SAVE
  !
  LOGICAL :: lspinorb
  !! if .TRUE. this is a spin-orbit calculation
  LOGICAL :: lforcet
  !! if .TRUE. apply Force Theorem to calculate MAE 
  LOGICAL :: starting_spin_angle
  !! if .TRUE. the initial wavefunctions are spin-angle functions. 
  LOGICAL :: domag
  !! if .TRUE. magnetization is computed
  !COMPLEX (DP) :: rot_ylm(lqmax,lqmax)
  !! transform real spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:)
  !! function needed to account for spinors.
  !
END MODULE spin_orb

MODULE relax
  !
  !! The variables used to control ionic relaxations
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(DP) :: epse = 0.0_dp
  !! threshold on total energy
  REAL(DP) :: epsf
  !! threshold on forces
  REAL(DP) :: epsp
  !! threshold on pressure
  REAL(DP) :: starting_scf_threshold
  !! self-explanatory
  !
END MODULE relax

MODULE wvfct
  !
  !! The variables needed to compute the band structure
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  INTEGER ::  npwx
  !! maximum number of PW for wavefunctions
  INTEGER ::  nbndx
  !! max number of bands use in iterative diag
  INTEGER ::  nbnd
  !! number of bands
  INTEGER ::  npw
  !! the number of plane waves
  INTEGER ::  current_k
  !! the index of k-point under consideration
  REAL(DP), ALLOCATABLE :: et(:,:)
  !! eigenvalues of the hamiltonian
  REAL(DP), ALLOCATABLE :: wg(:,:)
  !! the weight of each k point and band
  REAL(DP), ALLOCATABLE :: g2kin(:)
  !! kinetic energy
  INTEGER, ALLOCATABLE :: btype(:,:) 
  !! one if the corresponding state has to be
  !! converged to full accuracy, zero otherwise
  !
END MODULE wvfct

MODULE ener
  !
  !! The variables needed to compute the energies
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(DP) :: etot
  !! the total Kohn-Sham energy of the solid
  REAL(DP) :: hwf_energy
  !! this is the Harris-Weinert-Foulkes energy
  REAL(DP) :: eband
  !! the band energy
  REAL(DP) :: deband
  !! scf correction to have variational energy
  REAL(DP) :: ehart
  !! the Hartree energy
  REAL(DP) :: etxc
  !! the exchange and correlation energy
  REAL(DP) :: vtxc
  !! another exchange-correlation energy
  REAL(DP) :: etxcc
  !! the nlcc exchange and correlation
  REAL(DP) :: ewld
  !! the ewald energy
  REAL(DP) :: elondon
  !! the semi-empirical dispersion energy
  REAL(DP) :: edftd3
  !! the Grimme-d3 dispersion energy
  REAL(DP) :: exdm
  !! the XDM dispersion energy
  REAL(DP) :: demet
  !! variational correction ("-TS") for metals
  REAL(DP) :: epaw
  !! sum of one-center paw contributions
  REAL(DP) :: ef
  !! the Fermi energy
  REAL(DP) :: ef_up
  !! the Fermi energy up (if two_fermi_energies=.TRUE.)
  REAL(DP) :: ef_dw
  !! the Fermi energy down (if two_fermi_energies=.TRUE.)
  !
END MODULE ener

MODULE fixed_occ
  !
  !! The quantities needed in calculations with fixed occupations.
  !
  USE kinds,      ONLY : DP
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: f_inp(:,:)
  !! the occupations for each spin
  LOGICAL :: tfixed_occ
  !! if .TRUE. the occupations are fixed.
  LOGICAL :: one_atom_occupations
  !! if .TRUE. the occupations are decided according to the projections of the
  !! wavefunctions on the initial atomic  wavefunctions (to be used only
  !! for an isolated atom).
  !
END MODULE fixed_occ

MODULE paw_variables
    !
    USE kinds,      ONLY : DP
    !
    IMPLICIT NONE
    PUBLIC
    SAVE

    !!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Control flags: !!!!

    ! Set to true after initialization, to prevent double allocs:
    LOGICAL              :: paw_is_init = .false.
    ! Analogous to okvan in  "uspp_param" (Modules/uspp.f90)
    LOGICAL :: &
         okpaw = .FALSE.          ! if .TRUE. at least one pseudo is PAW

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Pseudopotential data: !!!!

    ! There is (almost) no pseudopotential data here, it is all stored in the upf type.
    ! See files pseudo_types.f90 and read_uspp.f90

    ! Constant to be added to etot to get all-electron energy 
    REAL(DP) :: total_core_energy = 0._dp
    ! true if all the pseudopotentials are PAW
    LOGICAL  :: only_paw

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Initialization data: !!!!

    INTEGER,PARAMETER    :: lm_fact = 3   ! To converge E_xc integrate up to LM = lm_fact * lm_max
    INTEGER,PARAMETER    :: lm_fact_x = 3 ! As above, for gradient corrected functionals
    INTEGER,PARAMETER    :: xlm = 2       ! Additional factor to add to have a good grad.corr.
    INTEGER,PARAMETER    :: radial_grad_style = 0 ! = 0 or 1, algorithm to use for d/dr

    TYPE paw_radial_integrator
        ! the following variables are used to integrate radial sampling
        INTEGER          :: lmax        ! max l component that can be integrated correctly
        INTEGER          :: ladd        ! additional l max that have been added for grad.corr.
        INTEGER          :: lm_max      ! as above, but +1 and squared
        INTEGER          :: nx          ! number of integration directions
        REAL(DP),POINTER :: ww(:)       ! integration weights (one per direction)
        REAL(DP),POINTER :: ylm(:,:)    ! Y_lm(nx,lm_max)
        REAL(DP),POINTER :: wwylm(:,:)  ! ww(nx) * Y_lm(nx,lm_max)
        ! additional variables for gradient correction
        REAL(DP),POINTER :: dylmt(:,:),&! |d(ylm)/dtheta|**2
                            dylmp(:,:)  ! |d(ylm)/dphi|**2
        REAL(DP),POINTER :: cos_phi(:)  ! cos(phi)  
        REAL(DP),POINTER :: sin_phi(:)  ! sin(phi)  
        REAL(DP),POINTER :: cos_th(:)  ! cos(theta)  (for divergence)
        REAL(DP),POINTER :: sin_th(:)  ! sin(theta)  (for divergence)
        REAL(DP),POINTER :: cotg_th(:)  ! cos(theta)/sin(theta)  (for divergence)
    END TYPE
    TYPE(paw_radial_integrator), ALLOCATABLE :: &
        rad(:) ! information to integrate different atomic species

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! self-consistent variables: !!!!

    ! This type contains some useful data that has to be passed to all
    ! functions, but cannot stay in global variables for parallel:
    TYPE paw_info
        INTEGER :: a ! atom index
        INTEGER :: t ! atom type index = itype(a)
        INTEGER :: m ! atom mesh = g(t)%mesh
        INTEGER :: b ! number of beta functions = upf(t)%nbeta
        INTEGER :: l ! max angular index l+1 -> (l+1)**2 is max
                     ! lm index, it is used to allocate rho
        INTEGER :: ae ! tells if we are doing all-electron (1) or pseudo (2)
    END TYPE

    ! Analogous to deeq in "uspp_param" (Modules/uspp.f90)
    REAL(DP), ALLOCATABLE :: &
         ddd_paw(:,:,:)  ! D: D^1_{ij} - \tilde{D}^1_{ij} (only Hxc part)

    REAL(DP), ALLOCATABLE ::  vs_rad(:,:,:)

 END MODULE paw_variables

module tsvdw_module
  use kinds,only : dp
  ! PUBLIC variables 
  !
  LOGICAL, PUBLIC :: vdw_isolated    ! isolated system control
  REAL(DP), PUBLIC:: vdw_econv_thr   ! energy convergence threshold for periodic systems
  REAL(DP), PUBLIC :: EtsvdW                                   !the TS-vdW energy
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: UtsvdW        !the TS-vdW wavefunction forces (dispersion potential)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: FtsvdW      !the TS-vdW ionic forces (-dE/dR)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: HtsvdW      !the TS-vdW cell forces (dE/dh)
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: VefftsvdW     !the TS-vdW effective Hirshfeld volume
end module tsvdw_module

MODULE extfield
  !! The quantities needed in calculations with external field.
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: tefield
  !! if .TRUE. a finite electric field is added to the
  !! local potential
  LOGICAL :: dipfield
  !! if .TRUE. the dipole field is subtracted
  ! TB
  LOGICAL :: relaxz
  !! relax in z direction
  LOGICAL :: block
  !! add potential barrier
  LOGICAL :: gate
  !! if .TRUE. and system is charged, charge is represented
  !! with charged plate (gate)
  !
  INTEGER :: edir
  !! direction of the field
  REAL(DP) :: emaxpos
  !! position of the maximum of the field (0<emaxpos<1)
  REAL(DP) :: eopreg
  !! amplitude of the inverse region (0<eopreg<1)
  REAL(DP) :: eamp
  !! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
  REAL(DP) :: etotefield
  !! energy correction due to the field
  REAL(DP) :: el_dipole
  !! electronic_dipole, used when dipole correction is on
  REAL(DP) :: ion_dipole
  !! ionic_dipole, used when dipole correction is on
  REAL(DP) :: tot_dipole
  !! total dipole, used when dipole correction is on
  ! TB
  REAL(DP) :: zgate
  !! position of charged plate
  REAL(DP) :: block_1
  !! ...blocking potential
  REAL(DP) :: block_2
  !! ...blocking potential
  REAL(DP) :: block_height
  !! ...blocking potential
  REAL(DP) :: etotgatefield
  !! energy correction due to the gate
  REAL(DP), ALLOCATABLE :: forcefield(:,:)
  !
  REAL(DP), ALLOCATABLE :: forcegate(:,:)
  !! TB gate forces
  !
end module extfield

MODULE pwcom
  !
  USE klist
  USE lsda_mod
  USE vlocal
  USE wvfct
  USE ener
  use fixed_occ
  USE relax
  USE spin_orb
  use fixed_occ
  !
END MODULE pwcom