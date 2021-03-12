module parameters
  !! This module contains parameters to control the actions of SCSH.
  !! Also routines to read the parameters and write them out again.
  use kinds,only : dp
  !use control
  use constants,only    : maxlen
  use environments,only : mkl_threads,lsetthreads
  implicit none
  
  integer,parameter :: npk = 40000 ! max number of k-points in pw.x calculation
  integer,parameter :: ntypx = 10  ! max number of different types of atom
  character(len=maxlen) :: scfoutname,phoutname,fildyn,epwoutname
  integer         :: naver
  integer         :: nstep
  integer         :: nsnap
  real(kind=dp)   :: dt
  real(kind=dp)   :: temp
  real(kind=dp)   :: gamma
  real(kind=dp)   :: init_kx,init_ky,init_kz  !激发后初始的电子(空穴)的k坐标
  integer(kind=dp):: init_band                !激发后初始的电子所处的能带
  integer         :: init_ikx,init_iky,init_ikz,init_ik 
  
  real(kind=dp)   :: w_laser     ! the laser energy in eV.
  real(kind=dp)   :: fwhm, fwhm_T            
  ! The laser shape is assumed to be a Gaussian f(t)=exp(-t^2/2T^2)
  !T is related to the full width at half-maximum as fwhm = 2sqrt(2ln2)T
  real(kind=dp)   :: Efield_x,Efield_y,Efield_z
  ! The laser Electric field strength along the xyz direction
  
  
  namelist / shinput / &
           scfoutname,phoutname,fildyn,epwoutname,naver,nstep,nsnap,gamma,&
           dt,temp,init_kx,init_ky,init_kz,init_band, w_laser,fwhm,fwhm_T,&
           Ex_laser,Ey_laser,Ez_laser,&
           lsetthreads,mkl_threads

  contains
  

  
end module parameters   
  