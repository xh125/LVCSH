!----------------------------------------------------------------------------!
! Copyright (C) 2020 WeiZhuang's group                                       !
! This file is distributed under the terms of the                            !
! GNU General Public License. See the file `License'                         !
! in the root directory of the present distribution,                         !
! or http://www.gnu.org/copyleft/gpl.txt .                                   !
!                                                                            !
!----------------------------------------------------------------------------!

program lvcsh
  !=====================================================================================================!
  != this program is used to simulate carrier transport with the self-consistent surface hopping method=!
  !=====================================================================================================!
  != of materials. We use a liner vibronic coupling model with electron-phonon couplings and           =!
  !=====================================================================================================!
  !=  system interaction with environment by system-bath interactions                                  =!     
  !=====================================================================================================!
  != Last updata 2020-12.30 Version:0.1.1                                                              =!   
  != Developed by XieHua at department of physic, USTC;xh125@mail.ustc.edu.cn                          =!
  !=====================================================================================================!
  !! Author: HuaXie
  !! Version: v0.1.1
  !! License: GNU
  !!=================================================================================================  
  use mkl_service
  use omp_lib
  use constants,only      : maxlen
  use environments,only   : environment_start,mkl_threads,&
                            set_mkl_threads,lsetthreads
  use readinput,only      : get_inputfile
  use readscf,only        : readpwscf_out
  use readphout,only      : readph_out
  use readepw,only        : readepwout
  use parameters, only    : lreadscfout,scfoutname,lreadphout,phoutname,epwoutname,temp,&
                            nsnap,nstep,dt,inputfilename,init_ik,llaser
  use hamiltonian,only    : H_nk,set_H_nk,set_H0_nk,calculate_eigen_energy_state
  use randoms,only        : init_random_seed
  use initialsh,only      : init_normalmode_coordinate_velocity,init_dynamical_variable
  use surfacehopping,only : iaver,isnap,istep,naver,phQ,phP,phQ0,phP0,e,p,e0,p0,d,d0,g,g1,&
                            allocatesh,celec_nk,chole_nk,w_e,w_h,w0_e,w0_h,&
                            calculate_nonadiabatic_coupling,convert_diabatic_adiabatic
  use elph2,only          : wf,nqtotf,nktotf,nbndfst
  use disp,only           : ph_configuration
  use modes,only          : nmodes
  use io      ,only       : stdout,io_time,time1
  use dynamics,only       : rk4_nuclei,rk4_electron_diabatic,lelec,lhole
  !use dynamica
  implicit none
  !===============!
  != preparation =!
  !===============!
  
  call environment_start( 'LVCSH' )
  call get_inputfile(inputfilename)
  if(lreadscfout) call readpwscf_out(scfoutname)
  if(lreadphout) call readph_out(phoutname)
  call readepwout(epwoutname)
  call set_H0_nk()
  !call ph_configuration(nqtotf,nmodes,wf,temp)
  call init_random_seed()
  if(lsetthreads) call set_mkl_threads(mkl_threads)
  call allocatesh(nmodes)
  
  
  !==========================!
  != loop over realizations =!
  !==========================!  
  do iaver=1,naver
    write(stdout,'(a,I4,a)') '###### iaver=',iaver,' ######'    
    !==================!
    != initialization =!
    !==================! 
    !!得到简正坐标的初始位置phQ和速度phP
    call init_normalmode_coordinate_velocity(wf,temp,phQ,phP)
    !!得到初始电子和空穴的状态
    call init_dynamical_variable(phQ,llaser,celec_nk,chole_nk,e,p,w_e,w_h)
    call calculate_nonadiabatic_coupling(nmodes,e,p,d)
    phQ0=phQ; phP0=phP; e0=e; p0=p; d0=d; w0_e=w_e; w0_h=w_h
    

    !=======================!
    != loop over snapshots =!
    !=======================!

    do isnap=1,nsnap
      do istep=1,nstep
        
        time1   = io_time()  
        !==============================!
        != update phQ,phP,c,e,p,d,w,g =!
        !==============================!
        
        !use rk4 to calculate the dynamical of phonon normal modes
        call rk4_nuclei(p0,phQ,phP,dt)
        call rk4_electron_diabatic(phQ0,celec_nk,dt,lelec)
        call rk4_electron_diabatic(phQ0,chole_nk,dt,lhole)
        call set_H_nk(phQ,H_nk)
        call calculate_eigen_energy_state(nktotf,nbndfst,H_nk,e,p)
        call calculate_nonadiabatic_coupling(nmodes,e,p,d)
        call convert_diabatic_adiabatic(p,celec_nk,w_e)
        call convert_diabatic_adiabatic(p,chole_nk,w_h)
        call calculate_hopping_probability(w0_e,phP0,d0,dt,g,g1)
!        !rk4方法计算电子空穴在透热表象下的演化
!        if(lelecsh) then
!          call dia_syH(nbasis,HH_e,E_e,P_e)
!          call calculate_nonadiabatic_coupling(nbasis,nfreem,E_e,P_e,Hep,d_e,dE_dQ_e)
!          call convert_diabatic_adiabatic(nbasis,P_e,C_e,W_e)
!          call calculate_hopping_probability(nbasis,nfreem,w0_e,phP0,d0_e,dt,g_e,g1_e,isurface_e)
!          call correct_hopping_probability(MSH,nbasis,g_e,g1_e,E0_e,w0_e,w_e,P0_e,P_e,isurface_e)
!          call nonadiabatic_transition(MSH,nbasis,nfreem,E0_e,P0_e,P_e,d_e,isurface_e,g_e,w0_e,phP,c_e)
!        endif
!        if(lholesh) then
!          call dia_syH(nbasis,HH_h,E_h,P_h)
!          call calculate_nonadiabatic_coupling(nbasis,nfreem,E_h,P_h,Hep,d_h,dE_dQ_h)
!          d_h     = -d_h
!          dE_dQ_h = -dE_dQ_h
!          call convert_diabatic_adiabatic(nbasis,P_h,C_h,W_h)
!          call calculate_hopping_probability(nbasis,nfreem,w0_h,phP0,d0_h,dt,g_h,g1_h,isurface_h)
!          call correct_hopping_probability(MSH,nbasis,g_h,g1_h,E0_h,w0_h,w_h,P0_h,P_h,isurface_h)
!          call nonadiabatic_transition(MSH,nbasis,nfreem,E0_h,P0_h,P_h,d_h,isurface_h,g_h,w0_h,phP,c_e)          
!        endif
!        
!        !===================!
!        != add bath effect =!
!        !===================!        
!        call add_bath_effect(nbasis,nfreem,d0_e,p0_e,d0_h,p0_h,womiga,dt,phQ,phP)
!
!        !============================!
!        != reset dynamical variable =!
!        !============================!
!
!        phQ0=phQ; phQ0=phQ; E0_e=E_e; p0_e=p_e; d0_e=d_e; w0_e=w_e;w0_h=w_h
!        E0_h=E_h; p0_h=p_h; d0_h=d_h
!        time2   = io_time()
!        write(stdout,"(A5,1X,I10,1X,F15.6,1X,A10)") 'Step=',istep,(time2-time1),"seconds"        
!      enddo
!      
!      !=====================!
!      != store information =!
!      !=====================!    
!
!      pes_elec(0,isnap,iaver)=E_e(isurface_e)
!      pes_hole(0,isnap,iaver)=E_h(isurface_h)
!      do ibasis=1,nbasis
!        csit_elec(ibasis,isnap)=csit_elec(ibasis,isnap)+abs(c_e(ibasis))**2
!        csit_hole(ibasis,isnap)=csit_hole(ibasis,isnap)+abs(c_h(ibasis))**2
!        wsit_elec(ibasis,isnap)=wsit_elec(ibasis,isnap)+abs(w_e(ibasis))**2
!        wsit_hole(ibasis,isnap)=wsit_hole(ibasis,isnap)+abs(w_h(ibasis))**2
!        psit_elec(ibasis,isnap)=psit_elec(ibasis,isnap)+p_e(ibasis,isurface_e)**2
!        psit_hole(ibasis,isnap)=psit_hole(ibasis,isnap)+p_h(ibasis,isurface_h)**2
!        pes_elec(ibasis,isnap,iaver)=e_e(ibasis)
!        pes_hole(ibasis,isnap,iaver)=e_h(ibasis)
!      enddo      
!      
!      do ifreem =1 , nfreem
!        !!平均核构型
!        xsit(ifreem,isnap)=xsit(ifreem,isnap)+phQ(ifreem)
!        !!平均动能
!        ksit(ifreem,isnap)=ksit(ifreem,isnap)+0.5d0*phP(ifreem)**2
      enddo
!           
    enddo
!    !!
  enddo
!  
!  csit_elec=csit_elec/naver
!  csit_hole=csit_hole/naver
!  wsit_elec=wsit_elec/naver
!  wsit_hole=wsit_hole/naver
!  psit_elec=psit_elec/naver
!  psit_hole=psit_hole/naver
!  xsit=xsit/naver*Au2ang*dsqrt(au2amu)
!  ksit=ksit/naver*Au2eV
!
!  !====================!
!  != save information =!
!  !====================!
!  call saveresult()
!  call cpu_time(t1)
!  write(6,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'
!  write(stdout,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'  
!  
!  close(stdout)
!  
!  stop
  
end program lvcsh


