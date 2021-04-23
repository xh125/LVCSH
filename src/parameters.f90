module parameters
  !! This module contains parameters to control the actions of SCSH.
  !! Also routines to read the parameters and write them out again.
  use kinds,only : dp
  !use control
  use constants,only    : maxlen
  use environments,only : mkl_threads,lsetthreads
  use lasercom,only     : llaser,efield,efield_cart,w_laser,fwhm
  use surfacecom,only   : methodsh
  implicit none
  
  integer,parameter :: npk = 40000 ! max number of k-points in pw.x calculation
  integer,parameter :: ntypx = 10  ! max number of different types of atom
  character(len=maxlen) :: scfoutname,phoutname,fildyn,epwoutname
  integer         :: naver
  integer         :: nstep
  integer         :: nsnap
  real(kind=dp)   :: dt
  real(kind=dp)   :: temp
  real(kind=dp)   :: gamma    ! gamma is the friction coefficient,dimention is 1/t(ps-1)  THZ
  real(kind=dp)   :: init_kx,init_ky,init_kz  !激发后初始的电子(空穴)的k坐标
  integer         :: init_cband,init_vband    !激发后初始的电子(空穴)所处的能带
  integer         :: init_ikx,init_iky,init_ikz,init_ik 
  character(len=maxlen) :: inputfilename = "LVCSH.in"
  logical :: lreadscfout,lreadphout,lreadfildyn
  
  
  namelist / shinput / &
           lreadscfout,lreadphout, scfoutname,phoutname,lreadfildyn,fildyn,epwoutname,&
           naver,nstep,nsnap,gamma,dt,temp,&
           init_kx,init_ky,init_kz,init_cband,init_vband,&
           llaser,efield,efield_cart,w_laser,fwhm,&
           lsetthreads,mkl_threads,methodsh

  contains
  

  
end module parameters   
  