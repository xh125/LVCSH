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
  integer         :: pre_nstep ! only phonons dynamica nstep before non-diabatic dynamica.
  real(kind=dp)   :: pre_dt
  real(kind=dp)   :: dt
  real(kind=dp)   :: temp
  real(kind=dp)   :: gamma    ! gamma is the friction coefficient,dimention is 1/t(ps-1)  THZ  
  real(kind=dp)   :: ld_fric  ! 2pi*gamma_qv/w_qv
  !calculate averager energy from simulation at different ld_gamma and reference temperature.
  !ref: 1 D. M. F. M. Germana Paterlini, Chemical Physics 236 (1998) 243. Table 1
  logical         :: l_gamma_energy
  real(kind=dp)   :: gamma_min,gamma_max
  real(kind=dp)   :: ld_fric_min,ld_fric_max
  integer         :: n_gamma
  real(kind=dp) 	:: lit_gmnvkq  !in unit of mev
	
  logical         :: l_ph_quantum
  logical :: lfeedback
  
	!ref:1	Bedard-Hearn, M. J., Larsen, R. E. & Schwartz, B. J. Mean-field dynamics with stochastic decoherence (MF-SD):
	!				a new algorithm for nonadiabatic mixed quantum/classical molecular-dynamics simulations with nuclear-induced decoherence.
	!				J Chem Phys 123, 234106, doi:10.1063/1.2131056 (2005).
	!ref:2 1	Granucci, G. & Persico, M. Critical appraisal of the fewest switches algorithm for surface hopping. 
	!				J Chem Phys 126, 134114, doi:10.1063/1.2715585 (2007).
	!ref:3 1	Zhu, C., Nangia, S., Jasper, A. W. & Truhlar, D. G.
	!				Coherent switching with decay of mixing: an improved treatment of electronic coherence for non-Born-Oppenheimer trajectories. 
	!				J Chem Phys 121, 7658-7670, doi:10.1063/1.1793991 (2004).
	!ref:4 1	Qiu, J., Bai, X. & Wang, L. Crossing Classified and Corrected Fewest Switches Surface Hopping.
	!					The Journal of Physical Chemistry Letters 9, 4319-4325, doi:10.1021/acs.jpclett.8b01902 (2018).
	logical :: ldecoherence
	real(kind=dp) :: cdecoherence  !0.1 Hartree
	
	
	
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
  
  ! phonons normal mode Langevin dynamica friction coefficient
  real(kind=dp),allocatable :: ld_gamma(:,:)
  
  ! phonons normal mode coordinate,and phonons P
  real(kind=dp),allocatable :: phQ(:,:),phP(:,:),phQ0(:,:),phP0(:,:)
  real(kind=dp),allocatable :: dEa_dQ(:,:),dEa_dQ_e(:,:),dEa_dQ_h(:,:)
  real(kind=dp),allocatable :: dEa2_dQ2(:,:),dEa2_dQ2_e(:,:),dEa2_dQ2_h(:,:)
  
  real(kind=dp),allocatable :: phQsit(:,:,:),phPsit(:,:,:),phKsit(:,:,:),phUsit(:,:,:)
  ! phonons normal mode 
  ! ref : <固体物理> (3-44) (3-45)
  real(kind=dp),allocatable :: phU(:,:),phK(:,:)
  real(kind=dp) :: SUM_phU,SUM_phK,SUM_phE
	real(kind=dp) :: SUM_phU0,SUM_phK0,SUM_phE0
  real(kind=dp) :: E_ph_CA_sum,E_ph_QA_sum
  ! averager crystal energy and temperature T,in classical and quantum . 
  
  real(kind=dp),allocatable :: d_e(:,:,:,:) ,d_h(:,:,:,:) ,g_e(:) ,g_h(:)
  real(kind=dp),allocatable :: d0_e(:,:,:,:),d0_h(:,:,:,:),g1_e(:),g1_h(:)  
  real(kind=dp),allocatable :: pes_e(:,:,:),csit_e(:,:),wsit_e(:,:),&
                               psit_e(:,:),mskds_e(:,:),mskds_sum_e(:,:),mskd_e(:),ipr_e(:),& 
                               pes_h(:,:,:),csit_h(:,:),wsit_h(:,:),&
                               psit_h(:,:),mskds_h(:,:),mskds_sum_h(:,:),mskd_h(:),ipr_h(:),&
															 apes_e(:,:),apes_sum_e(:,:),apes_h(:,:),apes_sum_h(:,:)
                               
  real(kind=dp),allocatable :: msd(:),ipr(:),msds(:,:)
  
  complex(kind=dpc),allocatable :: c_e(:),c_e_nk(:,:),w_e(:),w0_e(:)
  complex(kind=dpc),allocatable :: c_h(:),c_h_nk(:,:),w_h(:),w0_h(:)
  
  real(kind=dp) :: sumg0_e,sumg0_h,sumg1_e,sumg1_h
  
  
  
end module surfacecom