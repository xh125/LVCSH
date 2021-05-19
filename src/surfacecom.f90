module surfacecom
  use kinds,only : dp,dpc
  implicit none
  integer :: iaver
  integer :: isnap,istep
  integer :: iesurface,ihsurface,iesurface_j,ihsurface_j,&
             esurface_type,hsurface_type
  integer :: iesurface_,ihsurface_
  integer         :: naver
  integer         :: nstep
  integer         :: nsnap
  real(kind=dp)   :: dt
  real(kind=dp)   :: temp
  real(kind=dp)   :: gamma    ! gamma is the friction coefficient,dimention is 1/t(ps-1)  THZ  
  
  logical :: lelecsh
  logical :: lholesh
  logical :: lehpairsh
  
  integer :: ieband_min,ieband_max,ihband_min,ihband_max
  
  !method of surface hopping
  character(len=8) :: MethodSH
  ! FSSH, SC-FSSH, CC-FSSH
  ! FSSH    ref:1 J. C. Tully, J. Chem. Phys. 93 (1990) 1061.
  ! SC-FSSH ref:1 L. Wang, and O. V. Prezhdo, Journal of Physical Chemistry Letters 5 (2014) 713.
  ! CC-FSSH ref:
  
  ! phonons normal mode coordinate,and phonons P
  real(kind=dp),allocatable :: phQ(:,:),phP(:,:),phQ0(:,:),phP0(:,:)
  real(kind=dp),allocatable :: dEa_dQ(:,:),dEa_dQ_e(:,:),dEa_dQ_h(:,:)
  real(kind=dp),allocatable :: dEa2_dQ2(:,:),dEa2_dQ2_e(:,:),dEa2_dQ2_h(:,:)
  
  real(kind=dp),allocatable :: phQsit(:,:,:),phPsit(:,:,:),phKsit(:,:,:),phUsit(:,:,:)
  ! phonons normal mode 
  ! ref : <固体物理> (3-44) (3-45)
  real(kind=dp),allocatable :: phU(:,:),phK(:,:)
  real(kind=dp) :: SUM_phU,SUM_phK,SUM_phE

  !
  real(kind=dp),allocatable :: d_e(:,:,:,:) ,d_h(:,:,:,:) ,g_e(:) ,g_h(:)
  real(kind=dp),allocatable :: d0_e(:,:,:,:),d0_h(:,:,:,:),g1_e(:),g1_h(:)  
  real(kind=dp),allocatable :: pes_e(:,:,:),csit_e(:,:),wsit_e(:,:),&
                               psit_e(:,:),& 
                               pes_h(:,:,:),csit_h(:,:),wsit_h(:,:),&
                               psit_h(:,:) 
                               
  real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  
  complex(kind=dpc),allocatable :: c_e(:),c_e_nk(:,:),w_e(:),w0_e(:)
  complex(kind=dpc),allocatable :: c_h(:),c_h_nk(:,:),w_h(:),w0_h(:)
  
  real(kind=dp) :: sumg0_e,sumg0_h,sumg1_e,sumg1_h
  contains
  
  
end module surfacecom