module parameters
  !! This module contains parameters to control the actions of SCSH.
  !! Also routines to read the parameters and write them out again.
  use kinds,only : dp
  use constants,only    : maxlen
  use environments,only : mkl_threads,lsetthreads
  use lasercom,only     : llaser,efield,efield_cart,w_laser,fwhm
  use klist,only        : nelec
  use surfacecom,only   : methodsh,lfeedback,naver,nstep,nsnap,dt,pre_nstep,pre_dt,&
                          temp,gamma,ld_fric,l_ph_quantum,lit_gmnvkq,lit_ephonon,&
                          l_gamma_energy,gamma_min,gamma_max,ld_fric_min,ld_fric_max,n_gamma,&
                          lelecsh,lholesh,lehpairsh,ldecoherence,Cdecoherence,&
                          ieband_min,ieband_max,ihband_min,ihband_max,nefre_sh,nhfre_sh
	use modes,only : nmodes
  implicit none
  
  !integer,parameter :: npk = 40000 ! max number of k-points in pw.x calculation
  !integer,parameter :: ntypx = 10  ! max number of different types of atom
  character(len=maxlen) :: scfoutname,phoutname,fildyn,epwoutname

  real(kind=dp)   :: init_kx,init_ky,init_kz  !激发后初始的电子(空穴)的k坐标
  integer         :: init_hband,init_eband    !激发后初始的电子(空穴)所处的能带
  integer         :: init_ikx,init_iky,init_ikz,init_ik 
  real(kind=dp)   :: mix_thr
	! threshold to find the mixxing states
	character(len=maxlen) :: inputfilename = "LVCSH.in"
	character(len=maxlen) :: calculation
	character(len=maxlen) :: verbosity
	character(len=maxlen) :: outdir
	integer :: nnode,ncore,naver_sum,savedsnap
  logical :: lreadscfout,lreadphout,lreadfildyn,lsortpes
  
  namelist / shinput / &
           calculation,verbosity,outdir,ldecoherence,Cdecoherence,lit_gmnvkq,lit_ephonon,&
					 lreadscfout,lreadphout,scfoutname,phoutname,lreadfildyn,fildyn,&
					 epwoutname,methodsh,lfeedback,naver,nstep,nsnap,mix_thr,&
           pre_nstep,pre_dt,gamma,ld_fric,l_ph_quantum,dt,temp,&
           l_gamma_energy,gamma_min,gamma_max,ld_fric_min,ld_fric_max,n_gamma,&
           init_kx,init_ky,init_kz,init_hband,init_eband,&
           llaser,efield,efield_cart,w_laser,fwhm,&
           lsetthreads,mkl_threads,lelecsh,lholesh,lehpairsh,&
           nelec,ieband_min,ieband_max,ihband_min,ihband_max,nefre_sh,nhfre_sh,&
					 nnode,ncore,savedsnap,nmodes,lsortpes

  contains
  

  
end module parameters   