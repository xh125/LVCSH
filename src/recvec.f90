!=----------------------------------------------------------------------------=!
MODULE gvect
!=----------------------------------------------------------------------------=!

     ! ... variables describing the reciprocal lattice vectors
     ! ... G vectors with |G|^2 < ecutrho, cut-off for charge density
     ! ... With gamma tricks, G-vectors are divided into two half-spheres,
     ! ... G> and G<, containing G and -G (G=0 is in G>)
     ! ... This is referred to as the "dense" (or "hard", or "thick") grid

     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     INTEGER :: ngm  = 0  ! local  number of G vectors (on this processor)
                          ! with gamma tricks, only vectors in G>
     INTEGER :: ngm_g= 0  ! global number of G vectors (summed on all procs)
                          ! in serial execution, ngm_g = ngm
     INTEGER :: ngl = 0   ! number of G-vector shells
     INTEGER :: ngmx = 0  ! local number of G vectors, maximum across all procs

     REAL(DP) :: ecutrho = 0.0_DP ! energy cut-off for charge density 
     REAL(DP) :: gcutm = 0.0_DP   ! ecutrho/(2 pi/a)^2, cut-off for |G|^2

     INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
                           ! Needed in parallel execution: gstart=2 for the
                           ! proc that holds G=0, gstart=1 for all others

     !     G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
     !
     REAL(DP), ALLOCATABLE, TARGET :: gg(:) 

     !     gl(i) = i-th shell of G^2 (in units of tpiba2)
     !     igtongl(n) = shell index for n-th G-vector
     !
     REAL(DP), POINTER, PROTECTED            :: gl(:)
     INTEGER, ALLOCATABLE, TARGET, PROTECTED :: igtongl(:)
     !
     !     G-vectors cartesian components ( in units tpiba =(2pi/a)  )
     !
     REAL(DP), ALLOCATABLE, TARGET :: g(:,:) 

     !     mill = miller index of G vectors (local to each processor)
     !            G(:) = mill(1)*bg(:,1)+mill(2)*bg(:,2)+mill(3)*bg(:,3) 
     !            where bg are the reciprocal lattice basis vectors 
     !
     INTEGER, ALLOCATABLE, TARGET :: mill(:,:)
     
     !     ig_l2g  = converts a local G-vector index into the global index
     !               ("l2g" means local to global): ig_l2g(i) = index of i-th
     !               local G-vector in the global array of G-vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: ig_l2g(:)
     !
     !     mill_g  = miller index of all G vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_g(:,:)
     !
     ! the phases e^{-iG*tau_s} used to calculate structure factors
     !
     COMPLEX(DP), ALLOCATABLE :: eigts1(:,:), eigts2(:,:), eigts3(:,:)
     !
end module gvect

!=----------------------------------------------------------------------------=!
   MODULE gvecs
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     ! ... G vectors with |G|^2 < 4*ecutwfc, cut-off for wavefunctions
     ! ... ("smooth" grid). Gamma tricks and units as for the "dense" grid
     !
     INTEGER :: ngms = 0  ! local  number of smooth vectors (on this processor)
     INTEGER :: ngms_g=0  ! global number of smooth vectors (summed on procs) 
                          ! in serial execution this is equal to ngms
     INTEGER :: ngsx = 0  ! local number of smooth vectors, max across procs

     REAL(DP) :: ecuts = 0.0_DP   ! energy cut-off = 4*ecutwfc
     REAL(DP) :: gcutms= 0.0_DP   ! ecuts/(2 pi/a)^2, cut-off for |G|^2

     REAL(DP) :: dual = 0.0_DP    ! ecutrho=dual*ecutwfc
     LOGICAL  :: doublegrid = .FALSE. ! true if smooth and dense grid differ
                                      ! doublegrid = (dual > 4)
!=----------------------------------------------------------------------------=!
   END MODULE gvecs
!=----------------------------------------------------------------------------=!