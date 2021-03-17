module phdisp
  use kinds,only :dp,dpc
  implicit none
    INTEGER :: nq1, nq2, nq3  ! number of q-points in each direction
    INTEGER :: nqs            ! number of q points to be calculated
    integer :: iq
    REAL(DP), ALLOCATABLE :: x_q(:,:), & ! coordinates of the q points
                          wq(:) ! for plot
  
    REAL(DP), ALLOCATABLE :: omega_disp(:,:),omega_thz(:,:),omega_cm(:,:)
    real(kind=dp),allocatable :: ph_wqv(:,:)    !rad/s
    real(kind=dp),allocatable :: phQ(:,:),phP(:,:) !
    real(kind=dp),allocatable :: ph_nqv(:,:)  !1.0 / (EXP(hbar * ph_wqv(i, qp) / (K_BOLTZMANN_SI * T)) - 1)
    real(kind=dp),allocatable :: ph_lqv(:,:)  ! phonon filed amplitude Q_qv=sqrt(hbar/(2*ph_wqv)) * ph_lqv
    real(kind=dp),allocatable :: ph_pqv(:,:)  ! d_lqv/dt    P_qv=d(Q_qv)/dt=sqrt(hbar/(2*ph_wqv)) * d(ph_lqv)/dt
    real(kind=dp),allocatable :: ph_l(:,:),ph_l0(:,:)
    real(kind=dp),allocatable :: ph_p(:,:),ph_p0(:,:)
    
    
    complex(dpc),allocatable :: evecter_disp(:,:,:)
    
    complex(dpc),allocatable :: phid(:,:,:,:,:)
    !!! Dsas'a'(q)
    
    complex(dpc),allocatable :: disdyn(:,:,:,:)
    
  contains
  
  subroutine ph_configuration(nqtotf,nmodes,wf,T)
    use constants,only : hbar_SI,tpi,mev_Thz,K_BOLTZMANN_SI
    implicit none
    integer,intent(in) :: nqtotf
    integer,intent(in) :: nmodes
    real(kind=dp),intent(in) :: wf(nmodes,nqtotf)  ! in meV
    real(kind=dp),intent(in) :: T
    
    allocate(ph_wqv(nmodes,nqtotf))
    allocate(ph_nqv(nmodes,nqtotf),ph_lqv(nmodes,nqtotf),ph_pqv(nmodes,nqtotf))
    allocate(ph_l(nmodes,nqtotf),ph_p(nmodes,nqtotf))
    ph_wqv   = 0.0
    ph_nqv = 0.0
    ph_lqv = 0.0
    ph_pqv = 0.0
    
    ph_wqv = wf * mev_Thz * (1.0E12) * tpi
    ! correct frequency for phonons in SI (rad/s)
    ph_nqv = 1.0 / (EXP(hbar_SI * ph_wqv / (K_BOLTZMANN_SI * T)) - 1.0)
    ph_nqv(1:3,1) = 0.0 ! gamma qpoint (1-3) branch 沿xyz三个方向的平移

    ph_lqv = dsqrt(1.0+2.0*ph_nqv) !ph_Q = sqrt(hbar/(2*ph_wqv))*ph_lqv
    ph_pqv = ph_wqv*ph_lqv
    
  end subroutine ph_configuration
  
end module phdisp

module modes
  use kinds, only : dp
  implicit none
  INTEGER :: nirr, nmodes,irr
  integer :: imode,nu
  ! number of irreducible representations contained in the dynamical matrix
  ! number of modes
  INTEGER, ALLOCATABLE, TARGET :: npert(:) !3 * nat )
  ! the number of perturbations per IR
  INTEGER :: npertx
  ! max number of perturbations per IR
  COMPLEX (DP), POINTER :: &
       u(:,:),                     &!  3 * nat, 3 * nat),
       t(:,:,:,:),                 &! npertx, npertx, 48,3 * nat),
       tmq(:,:,:)                   ! npertx, npertx, 3 * nat)
  ! the transformation modes patterns
  ! the mode for deltarho
  ! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa

  CHARACTER(15), ALLOCATABLE :: name_rap_mode(:) ! symmetry type of each mode
  INTEGER, ALLOCATABLE :: num_rap_mode(:)  ! number of the representation for
                                           ! each mode
  !
end module 


MODULE efield_mod
  USE kinds, ONLY :  DP
  !
  ! ... the variables for the electric field perturbation
  !
  SAVE
  !
  REAL (DP) :: epsil (3, 3)
  REAL (DP), ALLOCATABLE :: &
       zstareu(:,:,:),       &! 3, 3, nat),
       zstarue(:,:,:)         ! 3, nat, 3)
  real(kind=dp),allocatable :: zeu(:,:,:) ! (3,3,nat)
  ! the dielectric constant
  ! the effective charges Z(E,Us) (E=scf,Us=bare)
  ! the effective charges Z(Us,E) (Us=scf,E=bare)
  COMPLEX (DP), ALLOCATABLE :: &
       zstareu0(:,:),        &! 3, 3 * nat),
       zstarue0(:,:),        &! 3 * nat, 3)
       zstarue0_rec(:,:)      ! 3 * nat, 3)
  ! the effective charges
  !
END MODULE efield_mod