module surfacecom
  use kinds,only : dp,dpc
  implicit none
  integer :: iaver
  integer :: isnap,istep
  integer :: iesurface,ihsurface
  integer :: ierr
  
  !method of surface hopping
  character(len=8) :: MethodSH
  ! FSSH, SC-FSSH, CC-FSSH
  ! FSSH    ref:1 J. C. Tully, J. Chem. Phys. 93 (1990) 1061.
  ! SC-FSSH ref:1 L. Wang, and O. V. Prezhdo, Journal of Physical Chemistry Letters 5 (2014) 713.
  ! CC-FSSH ref:
  ! phonons normal mode coordinate,and phonons P
  real(kind=dp),allocatable :: phQ(:,:),phP(:,:),phQ0(:,:),phP0(:,:)
  real(kind=dp),allocatable :: e(:),p(:,:),p_nk(:,:,:),d(:,:,:,:),ge(:),gh(:)
  real(kind=dp),allocatable :: e0(:),p0(:,:),d0(:,:,:,:),ge1(:),gh1(:)  
  real(kind=dp),allocatable :: pes(:,:,:),inf(:,:,:),csit(:,:),wsit(:,:),&
                               psit(:,:),xsit(:,:),ksit(:,:) 
  real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  complex(kind=dpc),allocatable :: celec_nk(:,:),w_e(:),w0_e(:)
  complex(kind=dpc),allocatable :: chole_nk(:,:),w_h(:),w0_h(:)
  real(kind=dp) :: minde_e,minde_h,sumg0_e,sumg0_h,sumg1_e,sumg1_h
  contains
  
  
end module surfacecom