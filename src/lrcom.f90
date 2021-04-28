!----------------------------------------------------------------------------
!
! ... Common variables for LR_Modules routines
!
MODULE lr_symm_base
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq=0, irotmq
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  LOGICAL :: minus_q, & ! if .TRUE. there is the symmetry sending q<->-q
             invsymq    ! if .TRUE. the small group of q has inversion
  logical :: l_nsymq_le_1! if nsymq <=1 l_nsymq_le_1 = .true.
  !
END MODULE lr_symm_base